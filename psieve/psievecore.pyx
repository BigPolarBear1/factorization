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
    cdef int inv[3][MAXSP]           # 1/(a*step*64) mod p for step 1, 2, 4
    cdef uint64_t** pat[3][MAXSP]    # pat[step index][s][t] -> p words
    cdef int stp_of[16]              # per K mod 16: X must lie in one class mod 1, 2 or 4 (0 = no X can work)
    cdef int res_of[16]
    cdef uint64_t mask2[16][64]      # per K mod 16 and X0 mod 64: positions whose quotient is a square mod 64
    cdef int32_t* smallp
    cdef uint8_t* smallbad
    cdef int nsmall
    cdef uint64_t ua
    cdef int a8, two_mode
    cdef int64_t a_c
    cdef int64_t n4q[MAXW]           # 4n mod each prime of a
    cdef readonly long long next_k, intervals, multipliers, survivors
    cdef readonly object sieve_primes

    def __init__(self, n, a, a_primes, int nsieve=18):
        cdef int j, p, t, x, s, i, o, kk, sh, np_ = 0, s64, cnt
        cdef int64_t q
        cdef uint64_t* base
        cdef uint64_t* wv
        cdef uint8_t* chi
        cdef int ainv, stepinv, nb, r, si, k16, e, cls8
        cdef uint8_t sq64[64]
        cdef uint8_t al[64]
        cdef int a64inv
        cdef uint64_t mk
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
        self.a_c = <int64_t> a
        self.a8 = self.absa % 8 if a > 0 else (8 - self.absa % 8) % 8
        self.M = []
        for j in range(self.w):
            q = a_primes[j]
            if q > 50_000_000:
                raise ValueError("prime of a too large for a root table")
            self.aq[j] = q
            self.n4q[j] = (4 * n) % q
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
                    # parity of K:  a = 1 mod 8: any K ; a = 5 mod 8: K odd ; a = 3 mod 4: K odd or K = 2 mod 4
                    self.two_mode = 0 if self.a8 == 1 else (1 if self.a8 == 5 else 2)
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
            chi = <uint8_t*> calloc(p, 1)        # 1 for squares and 0
            for x in range(p):
                chi[(x * x) % p] = 1
            nb = p + 64
            base = <uint64_t*> malloc(((nb >> 6) + 2) * sizeof(uint64_t))
            s64 = 64 % p
            for si in range(3):
                self.inv[si][np_] = pow((a * (1 << si) * 64) % p, -1, p)
                stepinv = pow((a * (1 << si)) % p, -1, p)       # residue x sits at position x/(a*step)
                self.pat[si][np_] = <uint64_t**> malloc(p * sizeof(uint64_t*))
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
                    self.pat[si][np_][t] = wv
            free(base); free(chi)
            np_ += 1
        self.nsp = np_
        free(comp)
        # 2-adic conditions: (X^2-4nK)/a must be a square mod 64. Depends on K mod 16 and X mod 64.
        for x in range(64):
            sq64[x] = 0
        for x in range(64):
            sq64[(x * x) & 63] = 1
        a64inv = pow(a % 64, -1, 64)
        for k16 in range(16):
            t = (4 * (n % 16) * k16) % 64
            cnt = 0
            for x in range(64):
                al[x] = sq64[(((x * x - t) % 64 + 64) % 64 * a64inv) % 64]
                cnt += al[x]
            if cnt == 0:
                self.stp_of[k16] = 0
                self.res_of[k16] = 0
                continue
            si = 0
            for j in (2, 1):                          # all live X in one class mod 4, else mod 2
                cls8 = -1
                r = 1
                for x in range(64):
                    if al[x]:
                        if cls8 < 0:
                            cls8 = x & ((1 << j) - 1)
                        elif cls8 != (x & ((1 << j) - 1)):
                            r = 0
                if r:
                    si = j
                    break
            self.stp_of[k16] = 1 << si
            self.res_of[k16] = 0
            if si:
                for x in range(64):
                    if al[x]:
                        self.res_of[k16] = x & ((1 << si) - 1)
                        break
            for x in range(64):
                mk = 0
                for e in range(64):
                    if al[(x + (self.absa % 64) * (1 << si) * e * (1 if a > 0 else 63)) % 64]:
                        mk |= (<uint64_t> 1) << e
                self.mask2[k16][x] = mk
        self.next_k = 1
        self.intervals = 0; self.multipliers = 0; self.survivors = 0

    def __dealloc__(self):
        cdef int j, t, si
        for j in range(self.w):
            free(self.sq[j])
        for si in range(3):
            for j in range(self.nsp):
                if self.pat[si][j] != NULL:
                    for t in range(self.sp[j]):
                        free(self.pat[si][j][t])
                    free(self.pat[si][j])
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
        cdef int sroot[MAXW]
        cdef long long nh = 0
        cdef uint64_t nterm[MAXW]
        cdef uint64_t ua = self.ua, bt, d, crm, x
        cdef uint64_t* pp[MAXSP]
        cdef int crp[MAXSP]
        cdef int aoff[MAXSP]
        cdef object offX
        cdef int k0, idx, lim, seg, stp, si, res, cr64, cr4, guard
        cdef uint64_t* pw
        cdef uint64_t m2, stepu
        cdef object stepX
        cdef uint64_t acc, lastmask
        cdef object T, cr, X, c, y, a = self.a, T4 = self.T4, absa = self.absa, M = self.M
        cdef long long nmult = 0
        hits = []
        lo = k_start
        while lo < k_stop and nh < max_hits:
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
                if i == 0 and self.two_mode == 2:
                    while m < hi - lo:
                        if ((lo + m) & 3) == 0:
                            bad[m] = 1
                        else:
                            remv[m] >>= 1
                        m += 2
                elif self.smallbad[i]:
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
                if nh >= max_hits or (max_multipliers and nmult >= max_multipliers):
                    break
                if bad[m]:
                    continue
                K = lo + m
                stp = self.stp_of[K & 15]
                if stp == 0:                         # no X can make the quotient a 2-adic square
                    continue
                good = 1
                for j in range(w):                   # 4nK must be a square mod every prime of a (all in C)
                    s = self.sq[j][(self.n4q[j] * (K % self.aq[j])) % self.aq[j]]
                    if s <= 0:
                        good = 0
                        break
                    sroot[j] = s
                if not good:
                    continue
                rem = remv[m]
                if rem > 1:                          # one large prime factor left (or a product of two: skip)
                    if rem > <long long> SMALL * SMALL:
                        continue
                    if jac(self.a_c % rem, rem) != 1:
                        continue
                T = T4 * K
                for j in range(w):
                    term[j] = <uint64_t> (sroot[j] * M[j] % absa)
                    nterm[j] = (ua - term[j]) % ua
                nmult += 1
                cr = isqrt(T)
                crm = <uint64_t> (cr % absa)
                si = 0 if stp == 1 else (1 if stp == 2 else 2)
                res = self.res_of[K & 15]
                stepu = ua * stp
                stepX = absa * stp
                offX = stepX * offset
                cr64 = <int> (cr % 64)
                cr4 = cr64 & 3
                for j in range(nsp):
                    p = self.sp[j]
                    pp[j] = self.pat[si][j][(self.n4p[j] * (K % p)) % p]
                    crp[j] = <int> (cr % p)
                    aoff[j] = <int> (((stepu % p) * (offset % p)) % p)
                lim = window
                if self.sgn < 0:                     # definite side: X runs from the crossing down to 0
                    c = cr // stepX + 1 - offset
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
                        guard = 0
                        while stp > 1 and ((cr4 + 1 + <int> (d & 3)) & (stp - 1)) != res and guard < 4:
                            d += ua                                  # move to the X of this class with the right residue mod step
                            guard += 1
                        m2 = self.mask2[K & 15][(cr64 + 1 + <int> (d & 63) + <int> ((stepu & 63) * (offset & 63))) & 63]
                    else:
                        d = (crm + ua - bt) % ua                     # X0 = cr - d
                        guard = 0
                        while stp > 1 and ((cr4 + 4 - <int> (d & 3)) & (stp - 1)) != res and guard < 4:
                            d += ua
                            guard += 1
                        m2 = self.mask2[K & 15][(cr64 + 128 - <int> (d & 63) - <int> (((stepu & 63) * (offset & 63)) & 63)) & 63]
                    acc = 0
                    for wd in range(nw):
                        iv[wd] = m2
                    for j in range(nsp):
                        p = self.sp[j]
                        if self.sgn > 0:
                            idx = <int> ((crp[j] + 1 + d % p + aoff[j]) % p)
                        else:
                            idx = <int> ((crp[j] + 2 * p - <int> (d % p) - aoff[j]) % p)
                        k0 = (idx * self.inv[si][j]) % p
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
                                X = cr + 1 + d + offX + stepX * e
                            else:
                                X = cr - d - offX - stepX * e
                                if X <= 0:
                                    continue
                            c = (X * X - T) // a
                            if c <= 0:
                                continue
                            y = isqrt(c)
                            if y * y == c:
                                hits.append((K, X, y))
                                nh += 1
            lo = lo + m if (nh >= max_hits or (max_multipliers and nmult >= max_multipliers)) else hi
            if max_multipliers and nmult >= max_multipliers:
                break
        self.next_k = lo
        self.multipliers += nmult
        free(bad); free(remv); free(iv)
        return hits
