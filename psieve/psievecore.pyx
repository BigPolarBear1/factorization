# cython: language_level=3, boundscheck=False, wraparound=False, cdivision=True, initializedcheck=False
"""
psievecore -- search  X^2 - a*y^2 = 4*n*K  over many multipliers K, in C.

a is the leftover of an SIQS relation (odd, squarefree, either sign, |a| < 2^62).
For each usable K the 2^w classes X = +-sqrt(4nK) mod each prime of a are sieved over a
short window of X that starts at the crossing X = sqrt(4nK) and moves away from it
(upwards for a > 0, downwards for a < 0), because that is where (X^2-4nK)/a is smallest
and squares are densest.

A multiplier K is used only if
  * every prime factor f of K has (a|f) = +1   (otherwise f divides both X and y: a scaled copy)
  * 4nK is a square mod every prime of a       (otherwise no X exists)

    core = PSieveCore(n, a, a_primes)
    hits = core.run(k_start, k_stop, window, max_hits)   ->  [(K, X, y), ...],  core.next_k
"""
from libc.stdlib cimport malloc, calloc, free
from libc.string cimport memset
from libc.stdint cimport uint64_t, int64_t, uint8_t, int32_t, uint32_t
from math import isqrt

cdef extern from *:
    int __builtin_ctzll(unsigned long long) nogil

cdef int jac(int64_t a, int64_t m) noexcept nogil:
    cdef int t = 1
    cdef int64_t tmp
    a %= m
    if a < 0:
        a += m
    while a:
        while (a & 1) == 0:
            a >>= 1
            if (m & 7) == 3 or (m & 7) == 5:
                t = -t
        tmp = a; a = m; m = tmp
        if (a & 3) == 3 and (m & 3) == 3:
            t = -t
        a %= m
    return t if m == 1 else 0

DEF MAXW = 12          # primes in a
DEF MAXSP = 40         # sieve primes
DEF SMALL = 65536      # trial division bound for the K filter

cdef class PSieveCore:
    cdef object n, a, absa, M, T4
    cdef int w, nsp, sgn
    cdef int64_t aq[MAXW]
    cdef int32_t* sq[MAXW]           # sq[j][t] = a square root of t mod aq[j], -1 if none
    cdef int sp[MAXSP]
    cdef int n4p[MAXSP]              # 4n mod p
    cdef int inv[MAXSP]              # 1/(a*64) mod p
    cdef uint64_t** pat[MAXSP]       # pat[s][t] -> p words
    cdef int32_t* smallp
    cdef uint8_t* smallbad
    cdef int nsmall
    cdef uint64_t ua
    cdef int a8
    cdef readonly long long next_k, intervals, multipliers, survivors
    cdef readonly object sieve_primes

    def __init__(self, n, a, a_primes, int nsieve=18):
        cdef int j, p, t, x, s, i, o, kk, sh, np_ = 0, s64
        cdef int64_t q
        cdef uint64_t* base
        cdef uint64_t* wv
        cdef uint8_t* chi
        cdef int ainv, stepinv, nb, r
        self.n = n; self.a = a; self.absa = abs(a); self.T4 = 4 * n
        self.sgn = 1 if a > 0 else -1
        self.w = len(a_primes)
        if self.w < 1 or self.w > MAXW or self.absa >= (1 << 62) or a % 2 == 0:
            raise ValueError("unsupported a")
        prod = 1
        for q in a_primes:
            prod *= q
        if prod != self.absa:
            raise ValueError("a must be squarefree and equal to the product of a_primes")
        self.ua = <uint64_t> self.absa
        self.a8 = self.absa % 8 if a > 0 else (8 - self.absa % 8) % 8
        self.M = []
        for j in range(self.w):
            q = a_primes[j]
            if q > 50_000_000:
                raise ValueError("prime of a too large for a root table")
            self.aq[j] = q
            self.sq[j] = <int32_t*> malloc(q * sizeof(int32_t))
            for t in range(q):
                self.sq[j][t] = -1
            for x in range(q):
                t = <int> ((<int64_t> x * x) % q)
                if self.sq[j][t] < 0:
                    self.sq[j][t] = x
            mq = self.absa // q
            self.M.append(mq * pow(mq, -1, q) % self.absa)
        # small primes for the K filter
        cdef uint8_t* comp = <uint8_t*> calloc(SMALL, 1)
        self.smallp = <int32_t*> malloc(7000 * sizeof(int32_t))
        self.smallbad = <uint8_t*> malloc(7000)
        self.nsmall = 0
        for i in range(2, SMALL):
            if not comp[i]:
                self.smallp[self.nsmall] = i
                if i == 2:
                    self.smallbad[self.nsmall] = 0 if self.a8 == 1 else 1
                else:
                    self.smallbad[self.nsmall] = 0 if jac(<int64_t> (a % i), i) == 1 else 1
                self.nsmall += 1
                o = i * 2
                while o < SMALL:
                    comp[o] = 1
                    o += i
        # sieve primes and their patterns
        self.sieve_primes = []
        for i in range(1, self.nsmall):
            p = self.smallp[i]
            if np_ == nsieve or np_ == MAXSP:
                break
            if n % p == 0 or a % p == 0:
                continue
            self.sp[np_] = p
            self.sieve_primes.append(p)
            self.n4p[np_] = (4 * n) % p
            ainv = pow(a % p, -1, p)
            self.inv[np_] = pow((a * 64) % p, -1, p)
            stepinv = ainv                       # residue x sits at position x/a
            chi = <uint8_t*> calloc(p, 1)        # 1 for squares and 0
            for x in range(p):
                chi[(x * x) % p] = 1
            self.pat[np_] = <uint64_t**> malloc(p * sizeof(uint64_t*))
            nb = p + 64
            base = <uint64_t*> malloc(((nb >> 6) + 2) * sizeof(uint64_t))
            s64 = 64 % p
            for t in range(p):
                memset(base, 0, ((nb >> 6) + 2) * sizeof(uint64_t))
                for x in range(p):
                    r = ((x * x - t) % p + p) % p
                    if chi[(r * ainv) % p]:
                        o = (x * stepinv) % p
                        while o < nb:
                            base[o >> 6] |= (<uint64_t> 1) << (o & 63)
                            o += p
                wv = <uint64_t*> malloc(p * sizeof(uint64_t))
                o = 0
                for j in range(p):
                    kk = o >> 6
                    sh = o & 63
                    wv[j] = (base[kk] >> sh) | ((base[kk + 1] << 1) << (63 - sh))
                    o += s64
                    if o >= p:
                        o -= p
                self.pat[np_][t] = wv
            free(base); free(chi)
            np_ += 1
        self.nsp = np_
        free(comp)
        self.next_k = 1
        self.intervals = 0; self.multipliers = 0; self.survivors = 0

    def __dealloc__(self):
        cdef int j, t
        for j in range(self.w):
            free(self.sq[j])
        for j in range(self.nsp):
            if self.pat[j] != NULL:
                for t in range(self.sp[j]):
                    free(self.pat[j][t])
                free(self.pat[j])
        free(self.smallp); free(self.smallbad)

    def run(self, long long k_start, long long k_stop, int window, int max_hits=1, long long max_multipliers=0, long long offset=0):
        """Search K in [k_start, k_stop). Stops early after max_hits hits (or max_multipliers usable K);
        self.next_k is where to resume. The window covers positions offset .. offset+window-1
        counted from the crossing, so a later call can extend an earlier one."""
        cdef long long lo, hi, K, CH = 1 << 16, m, rem
        cdef int nw = (window + 63) >> 6, i, j, s, f, cls, ncls = 1 << self.w, wd, bitpos, e
        cdef int w = self.w, nsp = self.nsp, p, good
        cdef uint8_t* bad = <uint8_t*> malloc(CH)
        cdef int64_t* remv = <int64_t*> malloc(CH * sizeof(int64_t))
        cdef uint64_t* iv = <uint64_t*> malloc(nw * sizeof(uint64_t))
        cdef uint64_t term[MAXW]
        cdef uint64_t nterm[MAXW]
        cdef uint64_t ua = self.ua, bt, d, crm, x
        cdef uint64_t* pp[MAXSP]
        cdef int crp[MAXSP]
        cdef int aoff[MAXSP]
        cdef object offX
        cdef int k0, idx, lim, seg
        cdef uint64_t* pw
        cdef uint64_t acc, lastmask
        cdef object T, cr, X, c, y, a = self.a, T4 = self.T4, absa = self.absa, M = self.M
        cdef long long nmult = 0
        hits = []
        offX = absa * offset
        for j in range(nsp):
            aoff[j] = <int> ((absa % self.sp[j]) * (offset % self.sp[j]) % self.sp[j])
        lo = k_start
        while lo < k_stop and len(hits) < max_hits:
            hi = lo + CH
            if hi > k_stop:
                hi = k_stop
            # ---- K filter for this chunk: strip small primes, mark K with a bad prime
            for m in range(hi - lo):
                bad[m] = 0
                remv[m] = lo + m
            for i in range(self.nsmall):
                f = self.smallp[i]
                if f > hi:
                    break
                m = ((lo + f - 1) // f) * f - lo
                if self.smallbad[i]:
                    while m < hi - lo:
                        bad[m] = 1
                        m += f
                else:
                    while m < hi - lo:
                        rem = remv[m]
                        while rem % f == 0:
                            rem //= f
                        remv[m] = rem
                        m += f
            for j in range(w):                       # K must be coprime to a
                m = ((lo + self.aq[j] - 1) // self.aq[j]) * self.aq[j] - lo
                while m < hi - lo:
                    bad[m] = 1
                    m += self.aq[j]
            for m in range(hi - lo):
                if len(hits) >= max_hits or (max_multipliers and nmult >= max_multipliers):
                    break
                if bad[m]:
                    continue
                K = lo + m
                rem = remv[m]
                if rem > 1:                          # one large prime factor left (or a product of two: skip)
                    if rem > <long long> SMALL * SMALL:
                        continue
                    if jac(<int64_t> (a % rem), rem) != 1:
                        continue
                T = T4 * K
                good = 1
                for j in range(w):
                    s = self.sq[j][<long long> (T % self.aq[j])]
                    if s <= 0:                       # non-residue, or K shares the prime (excluded above)
                        good = 0
                        break
                    term[j] = <uint64_t> (s * M[j] % absa)
                    nterm[j] = (ua - term[j]) % ua
                if not good:
                    continue
                nmult += 1
                cr = isqrt(T)
                crm = <uint64_t> (cr % absa)
                for j in range(nsp):
                    p = self.sp[j]
                    pp[j] = self.pat[j][(self.n4p[j] * (K % p)) % p]
                    crp[j] = <int> (cr % p)
                lim = window
                if self.sgn < 0:                     # definite side: X runs from the crossing down to 0
                    c = cr // absa + 1 - offset
                    if c <= 0:
                        continue
                    if c < lim:
                        lim = <int> c
                for cls in range(ncls):
                    bt = 0
                    for j in range(w):
                        bt += term[j] if (cls >> j) & 1 else nterm[j]
                        if bt >= ua:
                            bt -= ua
                    if self.sgn > 0:
                        d = (bt + ua - (crm + 1) % ua) % ua          # X0 = cr + 1 + d
                    else:
                        d = (crm + ua - bt) % ua                     # X0 = cr - d
                    acc = 0
                    for wd in range(nw):
                        iv[wd] = ~(<uint64_t> 0)
                    for j in range(nsp):
                        p = self.sp[j]
                        if self.sgn > 0:
                            idx = <int> ((crp[j] + 1 + d % p + aoff[j]) % p)
                        else:
                            idx = <int> ((crp[j] + 2 * p - <int> (d % p) - aoff[j]) % p)
                        k0 = (idx * self.inv[j]) % p
                        pw = pp[j]
                        wd = 0
                        while wd < nw:                       # straight runs between wrap-arounds, so the compiler can vectorise
                            seg = p - k0
                            if seg > nw - wd:
                                seg = nw - wd
                            for e in range(seg):
                                iv[wd + e] &= pw[k0 + e]
                            wd += seg
                            k0 = 0
                    self.intervals += 1
                    for wd in range(nw):
                        x = iv[wd]
                        while x:
                            bitpos = __builtin_ctzll(x)
                            x &= x - 1
                            e = (wd << 6) + bitpos
                            if e >= lim:
                                continue
                            self.survivors += 1
                            if self.sgn > 0:
                                X = cr + 1 + d + offX + absa * e
                            else:
                                X = cr - d - offX - absa * e
                                if X <= 0:
                                    continue
                            c = (X * X - T) // a
                            if c <= 0:
                                continue
                            y = isqrt(c)
                            if y * y == c:
                                hits.append((K, X, y))
            lo = lo + m if (len(hits) >= max_hits or (max_multipliers and nmult >= max_multipliers)) else hi
            if max_multipliers and nmult >= max_multipliers:
                break
        self.next_k = lo
        self.multipliers += nmult
        free(bad); free(remv); free(iv)
        return hits
