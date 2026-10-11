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
DEF NGP = 32           # sieve primes available (the first NGP odd primes)
DEF NWMAX = 128        # a pattern is stored with NWMAX extra words so a window can be read in one straight run

# Shared sieve patterns. For prime p, pattern[chi][t][bit i] says whether (i^2 - t) is 0 or has Legendre symbol
# +1 (chi = 0) / -1 (chi = 1). Position i of a window stands for X = X0 + (a*step)*i, and with u = a*step the
# condition "(X^2 - 4nK)/a is a square or 0 mod p" becomes exactly this with t = 4nK/u^2 and chi from (a|p),
# so one table serves every leftover.
cdef int _gp[NGP]
cdef uint64_t** _gpat[2][NGP]
cdef bint _gp_ready = False

cdef void _gp_build():
    global _gp_ready
    cdef int c = 0, x = 3, d, p, chi, t, i, r, o, kk, sh, s64, j, isp, nb
    cdef uint8_t* q
    cdef uint64_t* base
    cdef uint64_t* wv
    while c < NGP:
        isp = 1
        d = 3
        while d * d <= x:
            if x % d == 0:
                isp = 0
                break
            d += 2
        if isp:
            _gp[c] = x
            c += 1
        x += 2
    for c in range(NGP):
        p = _gp[c]
        q = <uint8_t*> calloc(p, 1)          # 1 = nonzero square
        for i in range(1, p):
            q[(i * i) % p] = 1
        nb = p + 64
        base = <uint64_t*> malloc(((nb >> 6) + 2) * sizeof(uint64_t))
        s64 = 64 % p
        for chi in range(2):
            _gpat[chi][c] = <uint64_t**> malloc(p * sizeof(uint64_t*))
            for t in range(p):
                memset(base, 0, ((nb >> 6) + 2) * sizeof(uint64_t))
                for i in range(p):
                    r = ((i * i - t) % p + p) % p
                    if r == 0 or (q[r] == 1) == (chi == 0):
                        o = i
                        while o < nb:
                            base[o >> 6] |= (<uint64_t> 1) << (o & 63)
                            o += p
                wv = <uint64_t*> malloc((p + NWMAX) * sizeof(uint64_t))
                o = 0
                for j in range(p):
                    kk = o >> 6
                    sh = o & 63
                    wv[j] = (base[kk] >> sh) | ((base[kk + 1] << 1) << (63 - sh))
                    o += s64
                    if o >= p:
                        o -= p
                for j in range(p, p + NWMAX):
                    wv[j] = wv[j - p]
                _gpat[chi][c][t] = wv
        free(q); free(base)
    _gp_ready = True
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
    cdef int usq[3][MAXSP]           # 1/(a*step)^2 mod p
    cdef int gidx[MAXSP]             # index of the sieve prime in the shared table
    cdef object grp                  # products of sieve primes below 2^62: one big-number reduction serves a whole group
    cdef int grp_of[MAXSP]
    cdef int ngrp
    cdef int gu[4][MAXSP]            # g*|a| mod p for g = 0..3
    cdef uint8_t* ktab[3][MAXSP]     # ktab[step][s][x] = x/(a*step*64) mod p
    cdef int gchi[MAXSP]             # 0 if a is a square mod p, 1 if not
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
        cdef int ainv, stepinv, nb, r, si, k16, e, cls8, u
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
        # sieve primes: taken from the shared table, skipping primes of n or a
        if not _gp_ready:
            _gp_build()
        self.sieve_primes = []
        for i in range(NGP):
            p = _gp[i]
            if np_ == nsieve or np_ == MAXSP:
                break
            if n % p == 0 or a % p == 0:
                continue
            self.sp[np_] = p
            self.gidx[np_] = i
            self.gchi[np_] = 0 if jac(<int64_t> (a % p), p) == 1 else 1
            self.sieve_primes.append(p)
            self.n4p[np_] = (4 * n) % p
            for si in range(3):
                u = (a * (1 << si)) % p
                self.inv[si][np_] = pow((u * 64) % p, -1, p)
                self.usq[si][np_] = pow((u * u) % p, -1, p)
                self.ktab[si][np_] = <uint8_t*> malloc(p)
                for x in range(p):
                    self.ktab[si][np_][x] = (x * self.inv[si][np_]) % p
            for x in range(4):
                self.gu[x][np_] = (x * (self.absa % p)) % p
            np_ += 1
        self.nsp = np_
        self.grp = []
        prodg = 1
        for j in range(np_):
            if prodg * self.sp[j] >= (1 << 62):
                self.grp.append(prodg)
                prodg = 1
            prodg *= self.sp[j]
            self.grp_of[j] = len(self.grp)
        self.grp.append(prodg)
        self.ngrp = len(self.grp)
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
                free(self.ktab[si][j])
        free(self.smallp); free(self.smallbad)

    def run(self, long long k_start, long long k_stop, int window, int max_hits=1, long long max_multipliers=0, long long offset=0, klist=None):
        """Search K in [k_start, k_stop). Stops early after max_hits hits (or max_multipliers usable K);
        self.next_k is where to resume. The window covers positions offset .. offset+window-1
        counted from the crossing, so a later call can extend an earlier one.
        With klist (an array('q') of multipliers already known to have only good prime factors) the range
        arguments are ignored and exactly those K are searched; next_k is left alone."""
        cdef long long lo, hi, K, CH = 1 << 16, m, rem
        cdef int nw = (window + 63) >> 6, i, j, s, f, cls, ncls = 1 << self.w, wd, bitpos, e
        cdef int w = self.w, nsp = self.nsp, p, good
        cdef uint8_t* bad = <uint8_t*> malloc(CH)
        cdef int64_t* remv = <int64_t*> malloc(CH * sizeof(int64_t))
        cdef uint64_t* iv = <uint64_t*> malloc(nw * sizeof(uint64_t))
        cdef uint64_t term[MAXW]
        cdef uint64_t tm2[MAXW]
        cdef uint64_t crg[8]
        cdef list grp = self.grp
        cdef int tmq[MAXW][MAXSP]
        cdef int dmq[MAXSP]
        cdef int baseq[MAXSP]
        cdef uint64_t dbase
        cdef int q, gray
        cdef bint up, wrap
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
        cdef long long[::1] kl
        cdef bint uselist = klist is not None
        hits = []
        if uselist:
            kl = klist
            k_start = 0
            k_stop = kl.shape[0]
        lo = k_start
        while lo < k_stop and nh < max_hits:
            hi = lo + CH
            if hi > k_stop:
                hi = k_stop
            # ---- K filter for this chunk: strip small primes, mark K with a bad prime
            for m in range(hi - lo if not uselist else 0):
                bad[m] = 0
                remv[m] = lo + m
            for i in range(self.nsmall if not uselist else 0):
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
            for j in range(w if not uselist else 0): # K must be coprime to a
                m = ((lo + self.aq[j] - 1) // self.aq[j]) * self.aq[j] - lo
                while m < hi - lo:
                    bad[m] = 1
                    m += self.aq[j]
            for m in range(hi - lo):
                if nh >= max_hits or (max_multipliers and nmult >= max_multipliers):
                    break
                if uselist:
                    K = kl[lo + m]
                else:
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
                rem = remv[m] if not uselist else 1
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
                cr64 = <int> (cr & 63)
                for j in range(self.ngrp):
                    crg[j] = <uint64_t> (cr % grp[j])
                cr4 = cr64 & 3
                for j in range(nsp):
                    p = self.sp[j]
                    pp[j] = _gpat[self.gchi[j]][self.gidx[j]][(((self.n4p[j] * (K % p)) % p) * self.usq[si][j]) % p]
                    crp[j] = <int> (crg[self.grp_of[j]] % p)
                    aoff[j] = <int> (((stepu % p) * (offset % p)) % p)
                lim = window
                if self.sgn < 0:                     # definite side: X runs from the crossing down to 0
                    c = cr // stepX + 1 - offset
                    if c <= 0:
                        continue
                    if c < lim:
                        lim = <int> c
                # Classes are walked in Gray-code order, so each one differs from the last in the sign of one root
                # and every residue is updated by an addition instead of being recomputed.
                bt = 0
                for j in range(w):
                    bt += nterm[j]
                    if bt >= ua:
                        bt -= ua
                    tm2[j] = (2 * term[j]) % ua
                if self.sgn > 0:
                    dbase = (bt + ua - (crm + 1) % ua) % ua           # X0 = cr + 1 + d
                else:
                    dbase = (crm + ua - bt) % ua                      # X0 = cr - d
                for q in range(nsp):
                    p = self.sp[q]
                    dmq[q] = <int> (dbase % p)
                    if self.sgn > 0:
                        baseq[q] = (crp[q] + 1 + aoff[q]) % p
                    else:
                        baseq[q] = (crp[q] + 2 * p - aoff[q]) % p
                    for j in range(w):
                        tmq[j][q] = <int> (tm2[j] % p)
                cls = 0
                for gray in range(ncls):
                    if gray:
                        j = __builtin_ctzll(<unsigned long long> gray)
                        cls ^= 1 << j
                        up = ((cls >> j) & 1) == (1 if self.sgn > 0 else 0)   # does d go up by 2*term[j] ?
                        if up:
                            dbase += tm2[j]
                            wrap = dbase >= ua
                            if wrap:
                                dbase -= ua
                            for q in range(nsp):
                                dmq[q] += tmq[j][q] - (self.gu[1][q] if wrap else 0)
                                if dmq[q] >= self.sp[q]:
                                    dmq[q] -= self.sp[q]
                                elif dmq[q] < 0:
                                    dmq[q] += self.sp[q]
                        else:
                            wrap = dbase >= tm2[j]                    # no borrow: the added |a| cancels again
                            dbase = dbase + ua - tm2[j]
                            if wrap:
                                dbase -= ua
                            for q in range(nsp):
                                dmq[q] += (0 if wrap else self.gu[1][q]) - tmq[j][q]
                                if dmq[q] >= self.sp[q]:
                                    dmq[q] -= self.sp[q]
                                elif dmq[q] < 0:
                                    dmq[q] += self.sp[q]
                    d = dbase
                    guard = 0
                    if self.sgn > 0:
                        while stp > 1 and ((cr4 + 1 + <int> (d & 3)) & (stp - 1)) != res and guard < 3:
                            d += ua                                  # move to the X of this class with the right residue mod step
                            guard += 1
                        m2 = self.mask2[K & 15][(cr64 + 1 + <int> (d & 63) + <int> ((stepu & 63) * (offset & 63))) & 63]
                    else:
                        while stp > 1 and ((cr4 + 4 - <int> (d & 3)) & (stp - 1)) != res and guard < 3:
                            d += ua
                            guard += 1
                        m2 = self.mask2[K & 15][(cr64 + 128 - <int> (d & 63) - <int> (((stepu & 63) * (offset & 63)) & 63)) & 63]
                    acc = 0
                    for wd in range(nw):
                        iv[wd] = m2
                    for j in range(nsp):
                        p = self.sp[j]
                        if self.sgn > 0:
                            idx = baseq[j] + dmq[j] + self.gu[guard][j]
                        else:
                            idx = baseq[j] + 2 * p - dmq[j] - self.gu[guard][j]
                        while idx >= p:
                            idx -= p
                        acc = 0
                        k0 = self.ktab[si][j][idx]
                        pw = pp[j] + k0
                        if nw <= NWMAX:                      # one straight run: the stored pattern is long enough
                            for wd in range(nw):
                                x = iv[wd] & pw[wd]
                                iv[wd] = x
                                acc |= x
                        else:
                            wd = 0
                            while wd < nw:
                                seg = NWMAX if nw - wd > NWMAX else nw - wd
                                for e in range(seg):
                                    x = iv[wd + e] & pw[e]
                                    iv[wd + e] = x
                                    acc |= x
                                wd += seg
                                k0 = (k0 + seg) % p
                                pw = pp[j] + k0
                        if acc == 0:                     # nothing left in this window: the remaining primes cannot matter
                            break
                    self.intervals += 1
                    if acc == 0:
                        continue
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
        if not uselist:
            self.next_k = lo
        self.multipliers += nmult
        free(bad); free(remv); free(iv)
        return hits


# ---------------------------------------------------------------------------------------------------------
# Multipliers in order of expected hits per unit of sieving.
# A squarefree K whose prime factors all split (a is a square mod each) carries 2^w ideals, so its chance of
# holding a solution is proportional to 2^w, while the crossing sqrt(4nK) grows like K^(1/2) and the hit density
# falls like K^(-1/4). value(K) = w*log(2) - log(K)/4. The search takes K in decreasing value.
from libc.math cimport log as clog
from libc.stdlib cimport realloc
import array as _array

cdef struct KBucket:
    long long* v
    Py_ssize_t n, cap

cdef inline void _kb_add(KBucket* b, long long K) noexcept nogil:
    if b.n == b.cap:
        b.cap = b.cap * 2 if b.cap else 256
        b.v = <long long*> realloc(b.v, b.cap * sizeof(long long))
    b.v[b.n] = K
    b.n += 1

cdef void _kdfs(int* good, double* lg, int ngood, int i, long long K, double v, double vhi, double vlo, double step,
                int a8, KBucket* bk, int nb) noexcept nogil:
    cdef int j, b
    cdef double v2
    cdef long long f
    if vlo <= v < vhi:
        b = <int> ((vhi - v) / step)
        if b < nb:
            _kb_add(&bk[b], K)
    if a8 == 1 or a8 == 3 or a8 == 7:        # the prime 2: splits (a = 1 mod 8) or ramifies (a = 3 mod 4)
        v2 = v + (0.5198603854199589 if a8 == 1 else -0.17328679513998632)
        if vlo <= v2 < vhi and K < (<long long> 1 << 61):
            b = <int> ((vhi - v2) / step)
            if b < nb:
                _kb_add(&bk[b], 2 * K)
    for j in range(i, ngood):
        f = good[j]
        v2 = v + 0.6931471805599453 - 0.25 * lg[j]
        if f > 16 and v2 + 0.52 < vlo:       # larger primes only lower the value further
            break
        if K > (<long long> 1 << 61) // f:
            break
        _kdfs(good, lg, ngood, j + 1, K * f, v2, vhi, vlo, step, a8, bk, nb)

def multipliers_by_value(good_primes, int a8, double vhi, double vlo, double step):
    """All squarefree products K of good_primes (ascending odd primes) with vlo <= value(K) < vhi,
    returned as a list of array('q'), one per value band of width step, best band first."""
    cdef int ngood = len(good_primes), j, nb = <int> ((vhi - vlo) / step) + 1
    cdef int* good = <int*> malloc(max(ngood, 1) * sizeof(int))
    cdef double* lg = <double*> malloc(max(ngood, 1) * sizeof(double))
    cdef KBucket* bk = <KBucket*> calloc(nb, sizeof(KBucket))
    for j in range(ngood):
        good[j] = good_primes[j]
        lg[j] = clog(<double> good[j])
    with nogil:
        _kdfs(good, lg, ngood, 0, 1, 0.0, vhi, vlo, step, a8, bk, nb)
    out = []
    for j in range(nb):
        arr = _array.array('q')
        if bk[j].n:
            arr.frombytes((<char*> bk[j].v)[:bk[j].n * sizeof(long long)])
        free(bk[j].v)
        out.append(arr)
    free(good); free(lg); free(bk)
    return out

def good_primes(a, int fmax):
    """Odd primes f < fmax with a a nonzero square mod f, ascending, as array('i')."""
    cdef int i, j
    cdef uint8_t* comp = <uint8_t*> calloc(fmax + 1, 1)
    cdef long long am
    out = _array.array('i')
    cdef list tmp = []
    i = 3
    while i < fmax:
        if not comp[i]:
            am = a % i
            if am != 0 and jac(am, i) == 1:
                tmp.append(i)
            if <long long> i * i < fmax:
                j = i * i
                while j < fmax:
                    comp[j] = 1
                    j += 2 * i
        i += 2
    free(comp)
    out.extend(tmp)
    return out

cdef class KOrder:
    """Hands out multipliers K for one leftover a in decreasing value, a chunk at a time.
    Only one value band is held in memory; park() drops even that (it is rebuilt on the next take)."""
    cdef object good, cur
    cdef int a8
    cdef double hi, cur_hi, step
    cdef Py_ssize_t pos
    cdef readonly long long used
    def __init__(self, a, int fmax=1_000_000, double step=0.05):
        cdef double top = 0.0
        self.good = good_primes(a, fmax)
        self.a8 = a % 8
        for f in self.good:
            if f < 16:
                top += 0.6931471805599453 - 0.25 * clog(<double> f)
        self.hi = top + (0.52 if self.a8 == 1 else 0.0) + 1e-9
        self.step = step
        self.cur = _array.array('q'); self.cur_hi = self.hi; self.pos = 0; self.used = 0
    cdef object _band(self, double vhi):
        out = _array.array('q')
        for part in multipliers_by_value(self.good, self.a8, vhi, vhi - self.step, self.step):
            out.extend(part)
        return out
    def take(self, Py_ssize_t count):
        """Next `count` multipliers (fewer only if a band ends), as array('q')."""
        if self.cur is None:                      # parked: rebuild the band we were in
            self.cur = self._band(self.cur_hi)
        while self.pos >= len(self.cur):
            self.cur_hi = self.hi
            self.cur = self._band(self.hi)
            self.hi -= self.step
            self.pos = 0
        out = self.cur[self.pos:self.pos + count]
        self.pos += len(out)
        self.used += len(out)
        return out
    def park(self):
        """Free the current band; the position in it is remembered."""
        if self.cur is not None and len(self.cur) > 4096:
            self.cur = None
