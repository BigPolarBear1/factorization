# cython: language_level=3, boundscheck=False, wraparound=False, cdivision=True, initializedcheck=False
"""
psievecore -- search  X^2 - a*y^2 = 4*n*K  over many multipliers K, in C.

a is the leftover of an SIQS relation (odd, squarefree, either sign, |a| < 2^96, prime factors below 5*10^7).
For each usable K the 2^w classes X = +-sqrt(4nK) mod each prime of a are sieved over a
short window of X that starts at the crossing X = sqrt(4nK) and moves away from it
(upwards for a > 0, downwards for a < 0), because that is where (X^2-4nK)/a is smallest
and squares are densest.

A multiplier K is used only if
  * every prime factor f of K has (a|f) = +1 (otherwise f divides both X and y: a scaled copy), or f is a prime of a taken once
  * 4nK is a square mod every prime of a       (otherwise no X exists)

    core  = PSieveCore(n, a, a_primes, 32)
    order = KOrder(a, 100000, 0.05, a_primes, n)          # usable K, best first
    hits  = core.run_fast(order.take(16384), 4096)        ->  [(K, X, y), ...]

run_fast() is the search used by QSv3_simd: everything per multiplier is machine arithmetic, |a| may be up to 96 bits.
run() is the older, slower routine; it is kept because it can also walk a plain range of K (|a| < 2^62 only).
"""
from libc.stdlib cimport malloc, calloc, free
from libc.string cimport memset
from libc.stdint cimport uint64_t, int64_t, uint8_t, int32_t, uint32_t
from math import isqrt

cdef extern from *:
    """
    typedef unsigned __int128 u128_t;
    typedef __int128 i128_t;
    """
    ctypedef unsigned long long u128_t
    ctypedef long long i128_t
    int __builtin_ctzll(unsigned long long) nogil
from libc.math cimport sqrtl

cdef inline int bmod(uint64_t x, uint64_t m, int p) noexcept nogil:
    # x mod p with a precomputed m = (2^64-1)//p
    cdef uint64_t r = x - <uint64_t> ((<u128_t> x * m) >> 64) * <uint64_t> p
    if r >= <uint64_t> p:
        r -= p
    return <int> r

cdef inline int bmod2(u128_t x, uint64_t m, int c64, int p) noexcept nogil:
    # the same for a 128-bit x (c64 = 2^64 mod p)
    if (x >> 64) == 0:
        return bmod(<uint64_t> x, m, p)
    return bmod(<uint64_t> (bmod(<uint64_t> (x >> 64), m, p) * c64 + bmod(<uint64_t> x, m, p)), m, p)

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

cdef inline uint64_t powmod(uint64_t a, uint64_t e, uint64_t m) noexcept nogil:
    cdef uint64_t r = 1
    a %= m
    while e:
        if e & 1:
            r = r * a % m
        a = a * a % m
        e >>= 1
    return r

cdef inline uint64_t powb(uint64_t a, uint64_t e, int p, uint64_t m) noexcept nogil:
    # a^e mod p for p < 2^31, with m = (2^64-1)//p (no division)
    cdef uint64_t r = 1
    while e:
        if e & 1:
            r = <uint64_t> bmod(r * a, m, p)
        a = <uint64_t> bmod(a * a, m, p)
        e >>= 1
    return r

cdef int64_t sqrtmod(uint64_t a, int p, uint64_t q, int s, uint64_t c, uint64_t bm) noexcept nogil:
    # Square root of a mod a prime p < 2^31, -1 if a is not a square. Tonelli-Shanks: p - 1 = q * 2^s,
    # c = g^q for a non-residue g, bm = (2^64-1)//p.
    cdef uint64_t r, t, t2, b
    cdef int m = s, i, k
    if s == 1:
        r = powb(a, (<uint64_t> p + 1) >> 2, p, bm)
        return <int64_t> r if <uint64_t> bmod(r * r, bm, p) == a else -1
    r = powb(a, (q + 1) >> 1, p, bm)
    t = powb(a, q, p, bm)
    while t != 1:
        i = 0
        t2 = t
        while t2 != 1:
            t2 = <uint64_t> bmod(t2 * t2, bm, p)
            i += 1
            if i == m:
                return -1
        b = c
        for k in range(m - i - 1):
            b = <uint64_t> bmod(b * b, bm, p)
        r = <uint64_t> bmod(r * b, bm, p)
        c = <uint64_t> bmod(b * b, bm, p)
        t = <uint64_t> bmod(t * c, bm, p)
        m = i
    return <int64_t> r

DEF MAXW = 12          # primes in a
DEF TABMAX = 8192      # primes of a up to here get a table of square roots (shared by all leftovers); larger ones are computed per multiplier
DEF NGP = 32           # sieve primes available (the first NGP odd primes)
DEF NWMAX = 128        # a pattern is stored with NWMAX extra words so a window can be read in one straight run

# Shared sieve patterns. For prime p, pattern[chi][t][bit i] says whether (i^2 - t) is 0 or has Legendre symbol
# +1 (chi = 0) / -1 (chi = 1). Position i of a window stands for X = X0 + (a*step)*i, and with u = a*step the
# condition "(X^2 - 4nK)/a is a square or 0 mod p" becomes exactly this with t = 4nK/u^2 and chi from (a|p),
# so one table serves every leftover.
cdef int _gp[NGP]
cdef uint8_t** _gmul[NGP]          # _gmul[i][m][x] = x*m mod p: every per-leftover table of residues is one of these
cdef uint8_t* _ginv[NGP]           # inverses mod p
cdef uint64_t _gbm[NGP]
cdef int _gc64[NGP]
cdef int32_t* _sqtab[TABMAX + 1]   # _sqtab[q][t] = a square root of t mod q, -1 if none (built on first use)
_m2cache = {}
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
        _gmul[c] = <uint8_t**> malloc(p * sizeof(uint8_t*))
        _ginv[c] = <uint8_t*> calloc(p, 1)
        for t in range(p):
            _gmul[c][t] = <uint8_t*> malloc(p)
            for i in range(p):
                _gmul[c][t][i] = (i * t) % p
                if (i * t) % p == 1:
                    _ginv[c][t] = i
        _gbm[c] = <uint64_t> 18446744073709551615 // <uint64_t> p
        _gc64[c] = <int> ((<uint64_t> 18446744073709551615 % <uint64_t> p + 1) % p)
    _gp_ready = True
DEF MAXSP = 40         # sieve primes
DEF SMALL = 65536      # trial division bound for the K filter

cdef class PSieveCore:
    cdef object n, a, absa, M, T4
    cdef int w, nsp, sgn
    cdef int64_t aq[MAXW]
    cdef int32_t* sq[MAXW]           # sq[j][t] = a square root of t mod aq[j], -1 if none (NULL for a prime above TABMAX)
    cdef uint64_t tsq[MAXW]          # Tonelli-Shanks constants for those primes: aq - 1 = tsq * 2^tss, tsc = g^tsq
    cdef uint64_t tsc[MAXW]
    cdef uint64_t tbm[MAXW]
    cdef int tss[MAXW]
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
    cdef int64_t anq[MAXW]           # -n*(a/q) mod q: decides a multiplier K that shares the prime q with a
    cdef readonly long long next_k, intervals, multipliers, survivors
    cdef uint64_t bm[MAXSP]          # (2^64-1)//p
    cdef int c64[MAXSP]              # 2^64 mod p
    cdef int kgc[3][4][MAXSP]        # ktab of g*|a| mod p
    cdef uint8_t* tix[3][MAXSP]      # tix[step][s][K mod p] = pattern row for this K
    cdef u128_t Mc2[MAXW]
    cdef u128_t ua2
    cdef bint big                    # |a| >= 2^62: only run_fast can be used
    cdef u128_t n4lo
    cdef long double n4f
    cdef public int t1
    cdef long long kfast             # K up to here: sqrt(4nK) is taken in machine arithmetic
    cdef readonly object sieve_primes

    def __init__(self, n, a, a_primes, int nsieve=18):
        cdef int j, p, t, x, s, i, o, kk, sh, np_ = 0, s64, cnt
        cdef int64_t q
        cdef uint64_t* base
        cdef uint64_t* wv
        cdef uint8_t* chi
        cdef int ainv, stepinv, nb, r, si, k16, e, cls8, u
        cdef int am, nm, absm, a64, sg63
        cdef uint8_t sq64[64]
        cdef uint8_t al[64]
        cdef int a64inv
        cdef uint64_t mk
        self.n = n; self.a = a; self.absa = abs(a); self.T4 = 4 * n
        self.sgn = 1 if a > 0 else -1
        self.w = len(a_primes)
        if self.w < 1 or self.w > MAXW or self.absa >= (1 << 96) or a % 2 == 0:
            raise ValueError("unsupported a")
        prod = 1
        for q in a_primes:
            prod *= q
        if prod != self.absa:
            raise ValueError("a must be squarefree and equal to the product of a_primes")
        self.big = self.big or self.absa >= (1 << 62)
        self.ua = <uint64_t> self.absa if not self.big else 0
        self.a_c = <int64_t> a if not self.big else 0
        self.a8 = self.absa % 8 if a > 0 else (8 - self.absa % 8) % 8
        self.M = []
        for j in range(self.w):
            q = a_primes[j]
            if q >= (1 << 31):
                raise ValueError("prime of a too large")
            self.aq[j] = q
            self.n4q[j] = (4 * n) % q
            self.anq[j] = (-n * (a // q)) % q
            if q > TABMAX:
                self.sq[j] = NULL
                self.big = True                       # run() needs the tables
                self.tss[j] = 0
                x = q - 1
                while x % 2 == 0:
                    x //= 2
                    self.tss[j] += 1
                self.tsq[j] = x
                t = 2
                while jac(t, q) != -1:
                    t += 1
                self.tsc[j] = pow(t, x, q)
                self.tbm[j] = <uint64_t> (((1 << 64) - 1) // q)
            else:
                if _sqtab[q] == NULL:
                    _sqtab[q] = <int32_t*> malloc(q * sizeof(int32_t))
                    for t in range(q):
                        _sqtab[q][t] = -1
                    for x in range(q):
                        t = <int> ((<int64_t> x * x) % q)
                        if _sqtab[q][t] < 0:
                            _sqtab[q][t] = x
                self.sq[j] = _sqtab[q]
            mq = self.absa // q
            self.M.append(mq * pow(mq, -1, q) % self.absa)
        self.smallp = NULL
        self.smallbad = NULL
        self.nsmall = 0
        self.two_mode = 0 if self.a8 == 1 else (1 if self.a8 == 5 else 2)
        # sieve primes: taken from the shared tables, skipping primes of n or a
        if not _gp_ready:
            _gp_build()
        self.sieve_primes = []
        for i in range(NGP):
            p = _gp[i]
            if np_ == nsieve or np_ == MAXSP:
                break
            am = a % p
            nm = n % p
            if nm == 0 or am == 0:
                continue
            self.sp[np_] = p
            self.gidx[np_] = i
            self.gchi[np_] = 0 if jac(am, p) == 1 else 1
            self.sieve_primes.append(p)
            self.n4p[np_] = (4 * nm) % p
            absm = am if self.sgn > 0 else p - am
            for si in range(3):
                u = (am << si) % p
                self.inv[si][np_] = _ginv[i][(u * 64) % p]
                self.usq[si][np_] = _ginv[i][(u * u) % p]
                self.ktab[si][np_] = _gmul[i][self.inv[si][np_]]
                self.tix[si][np_] = _gmul[i][(self.n4p[np_] * self.usq[si][np_]) % p]
            for x in range(4):
                self.gu[x][np_] = (x * absm) % p
                for si in range(3):
                    self.kgc[si][x][np_] = self.ktab[si][np_][self.gu[x][np_]]
            self.bm[np_] = _gbm[i]
            self.c64[np_] = _gc64[i]
            np_ += 1
        self.nsp = np_
        self.grp = None
        # 2-adic conditions: (X^2-4nK)/a must be a square mod 64. Depends on K mod 16 and X mod 64
        # (and only on a mod 64, its sign and n mod 16: computed once per combination).
        a64 = self.absa % 64
        sg63 = 1 if a > 0 else 63
        key = (a64, sg63, n % 16)
        if key in _m2cache:
            ent = _m2cache[key]
            for k16 in range(16):
                self.stp_of[k16] = ent[0][k16]
                self.res_of[k16] = ent[1][k16]
                for x in range(64):
                    self.mask2[k16][x] = ent[2][k16 * 64 + x]
        else:
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
                        if al[(x + a64 * (1 << si) * e * sg63) % 64]:
                            mk |= (<uint64_t> 1) << e
                    self.mask2[k16][x] = mk
            _m2cache[key] = ([self.stp_of[k16] for k16 in range(16)], [self.res_of[k16] for k16 in range(16)],
                             [self.mask2[k16][x] for k16 in range(16) for x in range(64)])
        self.next_k = 1
        self.intervals = 0; self.multipliers = 0; self.survivors = 0
        for j in range(self.w):
            mj = self.M[j]
            self.Mc2[j] = (<u128_t> (<uint64_t> (mj >> 64)) << 64) | <u128_t> (<uint64_t> (mj & 0xFFFFFFFFFFFFFFFF))
        self.ua2 = (<u128_t> (<uint64_t> (self.absa >> 64)) << 64) | <u128_t> (<uint64_t> (self.absa & 0xFFFFFFFFFFFFFFFF))
        T4 = self.T4
        self.n4lo = (<u128_t> (<uint64_t> ((T4 >> 64) & ((1 << 64) - 1))) << 64) | <u128_t> (<uint64_t> (T4 & ((1 << 64) - 1)))
        self.n4f = (<long double> (<uint64_t> ((T4 >> 128) & ((1 << 64) - 1))) * 18446744073709551616.0 + <long double> (<uint64_t> ((T4 >> 64) & ((1 << 64) - 1)))) * 18446744073709551616.0 + <long double> (<uint64_t> (T4 & ((1 << 64) - 1)))
        kf = (1 << 183) // T4 if T4 < (1 << 190) else 0
        self.kfast = kf if kf < (1 << 62) else (1 << 62)
        self.t1 = 16

    cdef inline int64_t rootq(self, int j, int64_t t) noexcept nogil:
        # a square root of t mod the j-th prime of a: positive if t is a nonzero square, else 0 or -1
        if self.sq[j] != NULL:
            return self.sq[j][t]
        if t == 0:
            return 0
        return sqrtmod(<uint64_t> t, <int> self.aq[j], self.tsq[j], self.tss[j], self.tsc[j], self.tbm[j])

    def __dealloc__(self):
        if self.smallp != NULL:
            free(self.smallp)
            free(self.smallbad)

    cdef _prep_run(self):
        # tables only the older run() needs
        cdef int i, o, j
        cdef uint8_t* comp
        if self.grp is not None:
            return
        comp = <uint8_t*> calloc(SMALL, 1)
        self.smallp = <int32_t*> malloc(7000 * sizeof(int32_t))
        self.smallbad = <uint8_t*> malloc(7000)
        self.nsmall = 0
        a = self.a
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
        free(comp)
        self.grp = []
        prodg = 1
        for j in range(self.nsp):
            if prodg * self.sp[j] >= (1 << 62):
                self.grp.append(prodg)
                prodg = 1
            prodg *= self.sp[j]
            self.grp_of[j] = len(self.grp)
        self.grp.append(prodg)
        self.ngrp = len(self.grp)

    def run(self, long long k_start, long long k_stop, int window, int max_hits=1, long long max_multipliers=0, long long offset=0, klist=None):
        """Search K in [k_start, k_stop). Stops early after max_hits hits (or max_multipliers usable K);
        self.next_k is where to resume. The window covers positions offset .. offset+window-1
        counted from the crossing, so a later call can extend an earlier one.
        With klist (an array('q') of multipliers already known to have only good prime factors) the range
        arguments are ignored and exactly those K are searched; next_k is left alone."""
        if self.big:
            if klist is None or window > 64 * NWMAX or offset:
                raise ValueError("|a| >= 2^62: only run_fast(klist, window <= %d) is available" % (64 * NWMAX))
            return self.run_fast(klist, window, max_hits)
        self._prep_run()
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
        cdef int act[MAXW]
        cdef int nact
        cdef long long kq
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
                nact = 0
                for j in range(w):                   # 4nK must be a square mod every prime of a (all in C)
                    kq = K % self.aq[j]
                    if kq == 0:
                        # K shares the prime q with a: X = 0 mod q is the only root (no sign choice for this prime),
                        # and (X^2-4nK)/a = -4n(K/q)/(a/q) mod q must itself be a square mod q
                        kq = (K // self.aq[j]) % self.aq[j]
                        if kq == 0 or self.sq[j][(self.anq[j] * kq) % self.aq[j]] <= 0:
                            good = 0
                            break
                        sroot[j] = 0
                        continue
                    s = self.sq[j][(self.n4q[j] * kq) % self.aq[j]]
                    if s <= 0:
                        good = 0
                        break
                    sroot[j] = s
                    act[nact] = j
                    nact += 1
                if not good:
                    continue
                ncls = 1 << nact
                rem = remv[m] if not uselist else 1
                if rem > 1:                          # one large prime factor left (or a product of two: skip)
                    if rem > <long long> SMALL * SMALL:
                        continue
                    if jac(self.a_c % rem, rem) != 1:
                        continue
                T = T4 * K
                for j in range(nact):                # only the primes with a sign choice
                    term[j] = <uint64_t> (sroot[act[j]] * M[act[j]] % absa)
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
                for j in range(nact):
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
                    for j in range(nact):
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

    def run_fast(self, klist, int window, int max_hits=1 << 30):
        """List mode of run() with every per-multiplier step in machine arithmetic (no big numbers until a
        position survives the whole sieve). Same hits as run(0, 0, window, max_hits, 0, 0, klist)."""
        cdef int nw = (window + 63) >> 6
        if nw > NWMAX or nw < 1:
            return self.run(0, 0, window, max_hits, 0, 0, klist)
        cdef long long[::1] kl = klist
        cdef Py_ssize_t nk = kl.shape[0], m
        cdef long long K, kq, nh = 0, nmult = 0
        cdef int w = self.w, nsp = self.nsp, p, good, i, j, s, q, gray, cls, ncls, nact, wd, bitpos, e
        cdef int stp, si, res, cr64, cr4, guard, lim, k0, kp, crp, tq, kt, kgw, v, sgn = self.sgn
        cdef uint64_t iv[NWMAX]
        cdef u128_t term[MAXW]
        cdef u128_t tm2[MAXW]
        cdef u128_t nterm[MAXW]
        cdef int sroot[MAXW]
        cdef int act[MAXW]
        cdef int dl[MAXW][4][MAXSP]
        cdef int kd[MAXSP]
        cdef int KB[MAXSP]
        cdef int spv[MAXSP]
        cdef uint64_t* pp[MAXSP]
        cdef int* dv
        cdef uint64_t* pw
        cdef int l1[MAXSP]
        cdef int l2[MAXSP]
        cdef int sp1[MAXSP]
        cdef int kg1[4][MAXSP]
        cdef int n1, n2, g, t1 = self.t1, g4, b
        cdef uint64_t* p0
        cdef uint64_t* p1
        cdef uint64_t* p2
        cdef uint64_t* p3
        cdef uint64_t* pw4[4]
        cdef uint64_t* pw8[8]
        cdef uint64_t* p4
        cdef uint64_t* p5
        cdef uint64_t* p6
        cdef uint64_t* p7
        cdef uint64_t x, m2, acc, crlo, crhi
        cdef u128_t ua = self.ua2, bt, d, crm, dbase, stepu      # |a| may exceed 64 bits: offsets inside a class are 128-bit
        cdef u128_t cr, Tlo, cq
        cdef i128_t r
        cdef long double Tf
        cdef bint up, wrap
        cdef object T, X, c, y, a = self.a, T4 = self.T4, absa = self.absa, crpy, stepX, dpy
        hits = []
        for q in range(nsp):
            spv[q] = self.sp[q]
        for m in range(nk):
            if nh >= max_hits:
                break
            K = kl[m]
            stp = self.stp_of[K & 15]
            if stp == 0:
                continue
            good = 1
            nact = 0
            for j in range(w):
                kq = K % self.aq[j]
                if kq == 0:
                    kq = (K // self.aq[j]) % self.aq[j]
                    if kq == 0 or self.rootq(j, (self.anq[j] * kq) % self.aq[j]) <= 0:
                        good = 0
                        break
                    sroot[j] = 0
                    continue
                s = <int> self.rootq(j, (self.n4q[j] * kq) % self.aq[j])
                if s <= 0:
                    good = 0
                    break
                sroot[j] = s
                act[nact] = j
                nact += 1
            if not good:
                continue
            ncls = 1 << nact
            nmult += 1
            # ---- cr = floor(sqrt(4nK))
            if K <= self.kfast:
                Tlo = self.n4lo * <u128_t> (<uint64_t> K)
                Tf = self.n4f * <long double> K
                cr = <u128_t> sqrtl(Tf)
                r = <i128_t> (Tlo - cr * cr)
                if r >= 0:
                    cr += <u128_t> (<long double> r / (2.0 * <long double> cr))
                else:
                    cr -= <u128_t> (<long double> (-r) / (2.0 * <long double> cr))
                r = <i128_t> (Tlo - cr * cr)
                while r < 0:
                    cr -= 1
                    r += <i128_t> (2 * cr + 1)
                while r > <i128_t> (2 * cr):
                    r -= <i128_t> (2 * cr + 1)
                    cr += 1
            else:
                crpy = isqrt(T4 * K)
                cr = (<u128_t> (<uint64_t> (crpy >> 64)) << 64) | <u128_t> (<uint64_t> (crpy & 0xFFFFFFFFFFFFFFFF))
            crlo = <uint64_t> cr
            crhi = <uint64_t> (cr >> 64)
            crm = cr % ua
            si = 0 if stp == 1 else (1 if stp == 2 else 2)
            res = self.res_of[K & 15]
            stepu = ua * <u128_t> (<uint64_t> stp)
            cr64 = <int> (crlo & 63)
            cr4 = cr64 & 3
            lim = window
            if sgn < 0:
                if cr < stepu * <u128_t> (<uint64_t> window):
                    cq = cr / stepu + 1
                    if <int> cq < lim:
                        lim = <int> cq
            bt = 0
            for j in range(nact):
                term[j] = (<u128_t> (<uint64_t> sroot[act[j]]) * self.Mc2[act[j]]) % ua
                nterm[j] = (ua - term[j]) % ua
                bt += nterm[j]
                if bt >= ua:
                    bt -= ua
                tm2[j] = (2 * term[j]) % ua
            if sgn > 0:
                dbase = (bt + ua - (crm + 1) % ua) % ua
            else:
                dbase = (crm + ua - bt) % ua
            # tier 1: the first t1 sieve primes that can remove anything for this K (a sieve prime dividing K passes
            # every position), tracked through the Gray walk. Tier 2: the rest, set up only for a class that survives.
            n1 = 0
            n2 = 0
            for q in range(nsp):
                p = spv[q]
                kp = bmod(<uint64_t> K, self.bm[q], p)
                if kp == 0:
                    continue
                if n1 == t1:
                    l2[n2] = q
                    n2 += 1
                    continue
                l1[n1] = q
                pp[n1] = _gpat[self.gchi[q]][self.gidx[q]][self.tix[si][q][kp]]
                crp = bmod(<uint64_t> (bmod(crhi, self.bm[q], p) * self.c64[q] + bmod(crlo, self.bm[q], p)), self.bm[q], p)
                if sgn > 0:
                    crp += 1
                    if crp == p:
                        crp = 0
                KB[n1] = self.ktab[si][q][crp]
                kd[n1] = self.ktab[si][q][bmod2(dbase, self.bm[q], self.c64[q], p)]
                sp1[n1] = p
                kgw = self.kgc[si][1][q]
                for g in range(4):
                    kg1[g][n1] = self.kgc[si][g][q]
                for j in range(nact):
                    kt = self.ktab[si][q][bmod2(tm2[j], self.bm[q], self.c64[q], p)]
                    dl[j][0][n1] = kt
                    v = kt - kgw
                    dl[j][1][n1] = v + p if v < 0 else v
                    dl[j][2][n1] = p - kt if kt else 0
                    v = kgw - kt
                    dl[j][3][n1] = v + p if v < 0 else v
                n1 += 1
            cls = 0
            for gray in range(ncls):
                if gray:
                    j = __builtin_ctzll(<unsigned long long> gray)
                    cls ^= 1 << j
                    up = ((cls >> j) & 1) == (1 if sgn > 0 else 0)
                    if up:
                        dbase += tm2[j]
                        wrap = dbase >= ua
                        if wrap:
                            dbase -= ua
                        dv = dl[j][1 if wrap else 0]
                    else:
                        wrap = dbase >= tm2[j]
                        dbase = dbase + ua - tm2[j]
                        if wrap:
                            dbase -= ua
                        dv = dl[j][2 if wrap else 3]
                    for q in range(n1):
                        v = kd[q] + dv[q]
                        kd[q] = v - sp1[q] if v >= sp1[q] else v
                d = dbase
                guard = 0
                if sgn > 0:
                    while stp > 1 and ((cr4 + 1 + <int> (<uint64_t> d & 3)) & (stp - 1)) != res and guard < 3:
                        d += ua
                        guard += 1
                    m2 = self.mask2[K & 15][(cr64 + 1 + <int> (<uint64_t> d & 63)) & 63]
                else:
                    while stp > 1 and ((cr4 + 4 - <int> (<uint64_t> d & 3)) & (stp - 1)) != res and guard < 3:
                        d += ua
                        guard += 1
                    m2 = self.mask2[K & 15][(cr64 + 128 - <int> (<uint64_t> d & 63)) & 63]
                acc = m2
                j = 0
                if n1 >= 8:                          # the first eight primes in one pass
                    for b in range(8):
                        p = sp1[b]
                        if sgn > 0:
                            k0 = KB[b] + kd[b] + kg1[guard][b]
                        else:
                            k0 = KB[b] + 2 * p - kd[b] - kg1[guard][b]
                        k0 = k0 - p if k0 >= p else k0
                        k0 = k0 - p if k0 >= p else k0
                        pw8[b] = pp[b] + k0
                    p0 = pw8[0]; p1 = pw8[1]; p2 = pw8[2]; p3 = pw8[3]; p4 = pw8[4]; p5 = pw8[5]; p6 = pw8[6]; p7 = pw8[7]
                    acc = 0
                    for wd in range(nw):
                        x = m2 & p0[wd] & p1[wd] & p2[wd] & p3[wd] & p4[wd] & p5[wd] & p6[wd] & p7[wd]
                        iv[wd] = x
                        acc |= x
                    j = 8
                else:
                    for wd in range(nw):
                        iv[wd] = m2
                if acc == 0:
                    j = n1
                while j < n1:                        # four primes per pass over the window
                    g4 = n1 - j
                    if g4 > 4:
                        g4 = 4
                    for b in range(g4):
                        p = sp1[j + b]
                        if sgn > 0:
                            k0 = KB[j + b] + kd[j + b] + kg1[guard][j + b]
                        else:
                            k0 = KB[j + b] + 2 * p - kd[j + b] - kg1[guard][j + b]
                        k0 = k0 - p if k0 >= p else k0
                        k0 = k0 - p if k0 >= p else k0
                        pw4[b] = pp[j + b] + k0
                    acc = 0
                    if g4 == 4:
                        p0 = pw4[0]; p1 = pw4[1]; p2 = pw4[2]; p3 = pw4[3]
                        for wd in range(nw):
                            x = iv[wd] & p0[wd] & p1[wd] & p2[wd] & p3[wd]
                            iv[wd] = x
                            acc |= x
                    else:
                        for b in range(g4):
                            pw = pw4[b]
                            acc = 0
                            for wd in range(nw):
                                x = iv[wd] & pw[wd]
                                iv[wd] = x
                                acc |= x
                    if acc == 0:
                        break
                    j += g4
                self.intervals += 1
                if acc == 0:
                    continue
                for j in range(n2):
                    q = l2[j]
                    p = spv[q]
                    crp = bmod(<uint64_t> (bmod(crhi, self.bm[q], p) * self.c64[q] + bmod(crlo, self.bm[q], p)), self.bm[q], p)
                    if sgn > 0:
                        k0 = crp + 1 + bmod2(d, self.bm[q], self.c64[q], p)
                    else:
                        k0 = crp + p - bmod2(d, self.bm[q], self.c64[q], p)
                    k0 = k0 - p if k0 >= p else k0
                    k0 = k0 - p if k0 >= p else k0
                    pw = _gpat[self.gchi[q]][self.gidx[q]][self.tix[si][q][bmod(<uint64_t> K, self.bm[q], p)]] + self.ktab[si][q][k0]
                    acc = 0
                    for wd in range(nw):
                        x = iv[wd] & pw[wd]
                        iv[wd] = x
                        acc |= x
                    if acc == 0:
                        break
                if acc == 0:
                    continue
                T = T4 * K
                crpy = (<object> crhi << 64) | <object> crlo
                dpy = (<object> (<uint64_t> (d >> 64)) << 64) | <object> (<uint64_t> d)
                stepX = absa * stp
                for wd in range(nw):
                    x = iv[wd]
                    while x:
                        bitpos = __builtin_ctzll(x)
                        x &= x - 1
                        e = (wd << 6) + bitpos
                        if e >= lim:
                            continue
                        self.survivors += 1
                        if sgn > 0:
                            X = crpy + 1 + dpy + stepX * e
                        else:
                            X = crpy - dpy - stepX * e
                            if X <= 0:
                                continue
                        c = (X * X - T) // a
                        if c <= 0:
                            continue
                        y = isqrt(c)
                        if y * y == c:
                            hits.append((K, X, y))
                            nh += 1
        self.multipliers += nmult
        return hits



# ---------------------------------------------------------------------------------------------------------
# Multipliers in order of expected hits per unit of sieving.
# A squarefree K whose prime factors all split (a is a square mod each) carries 2^w ideals, so its chance of
# holding a solution is proportional to 2^w, while the crossing sqrt(4nK) grows like K^(1/2) and the hit density
# falls like K^(-1/4). value(K) = w*log(2) - log(K)/4. The search takes K in decreasing value.
from libc.math cimport log as clog
from libc.stdlib cimport realloc
import array as _array

cdef int _opt_pow = 1, _opt_stp = 1
def set_opts(int p, int st):
    global _opt_pow, _opt_stp
    _opt_pow = p; _opt_stp = st

cdef struct KBucket:
    long long* v
    Py_ssize_t n, cap

cdef inline void _kb_add(KBucket* b, long long K) noexcept nogil:
    if b.n == b.cap:
        b.cap = b.cap * 2 if b.cap else 256
        b.v = <long long*> realloc(b.v, b.cap * sizeof(long long))
    b.v[b.n] = K
    b.n += 1

cdef struct KGen:
    int* good
    double* lg
    uint8_t* ram
    uint32_t* vec
    int* byc            # prime indices grouped by character, ascending inside a group
    int* cstart         # group c is byc[cstart[c] : cstart[c + 1]]
    int ngood, a8, nb
    double vhi, vlo, step, slack
    uint32_t vec2
    KBucket* bk

cdef inline void _kemit(KGen* g, long long K, double v, uint32_t need) noexcept nogil:
    # the prime 2 and the step in X it allows (a wider step is a longer window for the same work: +log(step)/2)
    #   a = 1 mod 8: 2 splits.  K odd or 2 mod 4 -> step 2, K = 0 mod 4 -> step 1
    #   a = 3 mod 4: 2 ramifies. K odd -> step 2, K = 2 mod 4 -> step 4
    #   a = 5 mod 8: K odd only, step 1
    # `need` is what the characters of the primes still to be multiplied in must add up to: K itself is usable
    # when need == 0, 2K (or 8K) when need == vec2.
    cdef double ve, vhi = g.vhi, vlo = g.vlo, step = g.step
    cdef long long Ke
    cdef int e, a8 = g.a8
    cdef KBucket* bk = g.bk
    if not _opt_stp:
        if vlo <= v < vhi and need == 0:
            _kb_add(&bk[<int> ((vhi - v) / step)], K)
        if a8 == 1 or a8 == 3 or a8 == 7:
            ve = v + (0.5198603854199589 if a8 == 1 else -0.17328679513998632)
            if vlo <= ve < vhi and K < (<long long> 1 << 60) and need == g.vec2:
                _kb_add(&bk[<int> ((vhi - ve) / step)], 2 * K)
    elif a8 == 1:
        ve = v + 0.34657359027997264
        if vlo <= ve < vhi and need == 0:
            _kb_add(&bk[<int> ((vhi - ve) / step)], K)
        Ke = K
        for e in range(1, 5):
            if Ke > (<long long> 1 << 60):
                break
            Ke *= 2
            ve = v + 0.6931471805599453 - 0.17328679513998632 * e + (0.34657359027997264 if e == 1 else 0.0)
            if vlo <= ve < vhi and need == (g.vec2 if e & 1 else 0):
                _kb_add(&bk[<int> ((vhi - ve) / step)], Ke)
    elif a8 == 5:
        if vlo <= v < vhi and need == 0:
            _kb_add(&bk[<int> ((vhi - v) / step)], K)
    else:
        ve = v + 0.34657359027997264
        if vlo <= ve < vhi and need == 0:
            _kb_add(&bk[<int> ((vhi - ve) / step)], K)
        ve = v + 0.6931471805599453 - 0.17328679513998632
        if vlo <= ve < vhi and K < (<long long> 1 << 60) and need == g.vec2:
            _kb_add(&bk[<int> ((vhi - ve) / step)], 2 * K)

cdef void _kdfs(KGen* g, int i, long long K, double v, uint32_t need) noexcept nogil:
    cdef int j, lo, hi, mid, c, pos, leaf = -1, ngood = g.ngood
    cdef double v2, ve, vlo = g.vlo, slack = g.slack
    cdef long long f, Ke
    cdef uint32_t nd, cc
    cdef double* lg = g.lg
    cdef int* good = g.good
    _kemit(g, K, v, need)
    for j in range(i, ngood):
        f = good[j]
        v2 = v + 0.6931471805599453 - 0.25 * lg[j]
        if f > 16 and v2 + slack < vlo:      # larger primes only lower the value further (slack = the most the prime 2 can add)
            break
        if K > (<long long> 1 << 60) // f:
            break
        # From some prime on, K*f can take no further factor (neither a larger prime nor f again). Those children
        # are leaves, and only the primes with the right character can be usable: they are read from the lists.
        if f > 16 and v2 - 0.25 * lg[j] + slack < vlo and (j + 1 == ngood or v2 + 0.6931471805599453 - 0.25 * lg[j + 1] + slack < vlo):
            leaf = j
            break
        Ke = K * f
        ve = v2
        nd = need ^ g.vec[j]
        while True:                          # f, f^2, f^3 ...: a higher power keeps the two ideals but costs size
            _kdfs(g, j + 1, Ke, ve, nd)
            if g.ram[j] or not _opt_pow or Ke > (<long long> 1 << 60) // f:      # a prime of a may divide K only once
                break
            Ke *= f
            nd ^= g.vec[j]
            ve -= 0.25 * lg[j]
            if ve + slack < vlo:
                break
    if leaf < 0:
        return
    for c in range(2):
        cc = need if c == 0 else need ^ g.vec2
        if c == 1 and g.vec2 == 0:
            break
        lo = g.cstart[cc]
        hi = g.cstart[cc + 1]
        pos = hi
        while lo < hi:                       # first prime of this character at or after `leaf`
            mid = (lo + hi) >> 1
            if g.byc[mid] >= leaf:
                hi = mid
            else:
                lo = mid + 1
        pos = lo
        hi = g.cstart[cc + 1]
        while pos < hi:
            j = g.byc[pos]
            f = good[j]
            v2 = v + 0.6931471805599453 - 0.25 * lg[j]
            if v2 + slack < vlo or K > (<long long> 1 << 60) // f:
                break
            _kemit(g, K * f, v2, need ^ g.vec[j])
            pos += 1

def multipliers_by_value(good_primes, int a8, double vhi, double vlo, double step, ramified=(), chars=None):
    """All products K of powers of good_primes (ascending odd primes) with vlo <= value(K) < vhi,
    returned as a list of array('q'), one per value band of width step, best band first.
    Primes listed in `ramified` (the primes of the leftover) are used at most once.
    chars = (vec per good prime, vec of 2, target): keep only K whose characters add up to the target, i.e.
    the K for which 4nK is a square at every prime of the leftover."""
    cdef int ngood = len(good_primes), j, nb = <int> ((vhi - vlo) / step) + 2, ncl = 1
    cdef KGen g
    cdef uint32_t tgt = 0, mx = 0
    cdef unsigned int[::1] cv
    g.good = <int*> malloc(max(ngood, 1) * sizeof(int))
    g.lg = <double*> malloc(max(ngood, 1) * sizeof(double))
    g.ram = <uint8_t*> calloc(max(ngood, 1), 1)
    g.vec = <uint32_t*> calloc(max(ngood, 1), sizeof(uint32_t))
    g.byc = <int*> malloc(max(ngood, 1) * sizeof(int))
    g.bk = <KBucket*> calloc(nb, sizeof(KBucket))
    g.vec2 = 0
    rs = set(ramified)
    for j in range(ngood):
        g.good[j] = good_primes[j]
        g.lg[j] = clog(<double> g.good[j])
        g.ram[j] = 1 if g.good[j] in rs else 0
    if chars is not None:
        cv = chars[0]
        for j in range(ngood):
            g.vec[j] = cv[j]
            if cv[j] > mx:
                mx = cv[j]
        g.vec2 = chars[1]; tgt = chars[2]
        if g.vec2 > mx:
            mx = g.vec2
        if tgt > mx:
            mx = tgt
    while <uint32_t> ncl <= mx:
        ncl *= 2
    g.cstart = <int*> calloc(ncl + 1, sizeof(int))
    for j in range(ngood):
        g.cstart[g.vec[j] + 1] += 1
    for j in range(ncl):
        g.cstart[j + 1] += g.cstart[j]
    cdef int* fill = <int*> malloc(ncl * sizeof(int))
    for j in range(ncl):
        fill[j] = g.cstart[j]
    for j in range(ngood):
        g.byc[fill[g.vec[j]]] = j
        fill[g.vec[j]] += 1
    free(fill)
    if _opt_stp:
        g.slack = 0.0 if a8 == 5 else (0.8665 if a8 == 1 else 0.52)
    else:
        g.slack = 0.52 if a8 == 1 else 0.0
    g.slack += 1e-6
    g.ngood = ngood; g.a8 = a8; g.nb = nb; g.vhi = vhi; g.vlo = vlo; g.step = step
    with nogil:
        _kdfs(&g, 0, 1, 0.0, tgt)
    out = []
    for j in range(nb):
        arr = _array.array('q')
        if g.bk[j].n:
            arr.frombytes((<char*> g.bk[j].v)[:g.bk[j].n * sizeof(long long)])
        free(g.bk[j].v)
        out.append(arr)
    free(g.good); free(g.lg); free(g.ram); free(g.bk); free(g.vec); free(g.byc); free(g.cstart)
    return out

def genus_chars(n, a, a_primes, good):
    """Characters for the K filter: bit j of vec[f] is set when multiplying K by the prime f flips the condition
    at the j-th prime q of a. For f = q itself (K shares the prime with a) the condition changes from
    (4nK|q) = 1 to (-n(a/q)(K/q)|q) = 1, and the bit records that change of target."""
    cdef int j, i, ng = len(good)
    cdef long long q, f
    cdef uint32_t t = 0, v2 = 0, b
    vec = _array.array('I', bytes(4 * ng))
    cdef unsigned int[::1] vv = vec
    qs = [int(x) for x in a_primes]
    for j in range(len(qs)):
        q = qs[j]
        tj = 1 if jac((n % q), q) == -1 else 0
        t |= tj << j
        if jac(2 % q, q) == -1:
            v2 |= 1 << j
        aj = 1 if jac(((-n * (a // q)) % q), q) == -1 else 0
        for i in range(ng):
            f = good[i]
            if f == q:
                b = tj ^ aj
            else:
                b = 1 if jac(f % q, q) == -1 else 0
            vv[i] |= b << j
    return vec, v2, t

cdef uint32_t _tdp[6600]
cdef uint64_t _tdm[6600]
cdef int _ntd = 0

cdef bint _isprime32(uint64_t v) noexcept nogil:
    # deterministic Miller-Rabin for v < 4759123141 (bases 2, 7, 61); v odd and > 61
    cdef uint64_t d = v - 1, x, b
    cdef int s = 0, i, j
    cdef bint comp
    while (d & 1) == 0:
        d >>= 1
        s += 1
    for i in range(3):
        b = 2 if i == 0 else (7 if i == 1 else 61)
        x = powmod(b, d, v)
        if x == 1 or x == v - 1:
            continue
        comp = True
        for j in range(s - 1):
            x = x * x % v
            if x == v - 1:
                comp = False
                break
        if comp:
            return False
    return True

cdef bint _prp2(uint64_t v) noexcept nogil:
    # base-2 Fermat test for odd v < 2^64 (only used to give up early on a large prime cofactor)
    cdef u128_t r = 1, b = 2
    cdef uint64_t e = v - 1
    while e:
        if e & 1:
            r = r * b % v
        b = b * b % v
        e >>= 1
    return r == 1

def odd_part(g):
    """g = sign * a * r^2 with a odd and squarefree: returns (primes of a ascending, r), or None when g is even to an
    odd power, or when what is left after dividing out the primes below 65536 is not 1 or a single prime below 2^31."""
    global _ntd
    cdef int i, j, e, isp
    cdef u128_t v
    cdef uint64_t v64, r = 1, q
    if _ntd == 0:
        i = 3
        while i < 65536:
            isp = 1
            j = 3
            while j * j <= i:
                if i % j == 0:
                    isp = 0
                    break
                j += 2
            if isp:
                _tdp[_ntd] = i
                _tdm[_ntd] = <uint64_t> 18446744073709551615 // <uint64_t> i
                _ntd += 1
            i += 2
    if g < 0:
        g = -g
    if g == 0 or (g >> 126):
        return None
    v = (<u128_t> (<uint64_t> (g >> 64)) << 64) | <u128_t> (<uint64_t> (g & 0xFFFFFFFFFFFFFFFF))
    e = 0
    while (<uint64_t> v & 1) == 0:
        v >>= 1
        e += 1
    if e & 1:
        return None
    r <<= e >> 1
    odd = []
    i = 0
    while i < _ntd and (v >> 64) != 0:
        q = _tdp[i]
        if v % q == 0:
            e = 0
            while v % q == 0:
                v = v / <u128_t> q                   # both operands 128-bit: with a 64-bit divisor Cython narrows v first
                e += 1
            if e & 1:
                odd.append(q)
            for j in range(e >> 1):
                r *= q
        i += 1
    if (v >> 64) != 0:
        return None
    v64 = <uint64_t> v
    while i < _ntd:
        q = _tdp[i]
        if q * q > v64:
            break
        if i == 1027 and v64 >= (<uint64_t> 1 << 31) and _prp2(v64):
            return None                              # what is left is one prime, and it is too large
        j = bmod(v64, _tdm[i], <int> q)              # may be off by one multiple of q for v64 near 2^64
        if j == 0 or j == <int> q:
            e = 0
            while v64 % q == 0:
                v64 //= q
                e += 1
            if e & 1:
                odd.append(q)
            for j in range(e >> 1):
                r *= q
        i += 1
    if v64 > 1:
        if v64 >= (<uint64_t> 1 << 31):
            return None
        odd.append(v64)
    return odd, r

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

cdef double _topw = 3.0

cdef class KOrder:
    """Hands out multipliers K for one leftover a in decreasing value, a chunk at a time.
    The tree of products is walked once per block of value bands; park() frees the block (it is rebuilt on demand)."""
    cdef object good, bands, ram, chars
    cdef int a8
    cdef double hi, block_hi, step, block
    cdef Py_ssize_t pos, band
    cdef readonly long long used
    def __init__(self, a, int fmax=1_000_000, double step=0.05, a_primes=None, n=None):
        cdef double top = 0.0
        self.good = good_primes(a, fmax)
        self.ram = tuple(int(q) for q in a_primes if q > 2) if a_primes else ()
        if a_primes:
            # a prime of a may divide K once: it halves the number of root classes instead of doubling the ideals,
            # so per unit of sieving it is worth exactly as much as a split prime
            merged = sorted(set(self.good) | set(int(q) for q in a_primes if q > 2 and q < (1 << 30)))
            self.good = _array.array('i', merged)
        self.chars = genus_chars(n, a, a_primes, self.good) if (n is not None and a_primes) else None
        self.a8 = a % 8
        for f in self.good:
            if f < 16:
                top += 0.6931471805599453 - 0.25 * clog(<double> f)
        self.hi = top + 0.87 + 1e-9
        self.step = step; self.block = 0.5
        self.bands = None; self.block_hi = self.hi; self.band = 0; self.pos = 0; self.used = 0
    def take(self, Py_ssize_t count):
        """Next `count` multipliers (fewer only if a band ends), as array('q')."""
        if self.bands is None:                    # first use, or parked: rebuild the block we were in
            self.bands = multipliers_by_value(self.good, self.a8, self.block_hi, self.block_hi - self.block, self.step, self.ram, self.chars)
        while self.band >= len(self.bands) or self.pos >= len(self.bands[self.band]):
            if self.band < len(self.bands):
                self.band += 1
                self.pos = 0
                continue
            total = sum(len(b) for b in self.bands)
            self.block_hi -= self.block
            if total > 1_000_000 and self.block > self.step * 1.5:
                self.block /= 2                   # keep a block to about a million multipliers
            self.bands = multipliers_by_value(self.good, self.a8, self.block_hi, self.block_hi - self.block, self.step, self.ram, self.chars)
            self.band = 0
            self.pos = 0
        out = self.bands[self.band][self.pos:self.pos + count]
        self.pos += len(out)
        self.used += len(out)
        return out
    def top(self, Py_ssize_t count):
        """The best `count` multipliers in one go (for a leftover that is looked at once and dropped)."""
        global _topw
        cdef double width = _topw - 0.5 if _topw > 1.5 else 1.0     # start a little above where the last leftover ended
        while True:
            bands = multipliers_by_value(self.good, self.a8, self.hi, self.hi - width, self.step, self.ram, self.chars)
            if sum(len(b) for b in bands) >= count or width > 14.0:
                break
            width += 0.5
        _topw = width
        out = _array.array('q')
        for b in bands:
            out.extend(b)
            if len(out) >= count:
                break
        return out[:count]
    def save(self):
        """Position in the ordering, to continue later with a new KOrder for the same leftover."""
        return (self.block_hi, self.block, self.band, self.pos, self.used)
    def load(self, st):
        self.block_hi, self.block, self.band, self.pos, self.used = st
        self.bands = None
    def park(self):
        """Free the current block; the position in it is remembered."""
        if self.bands is not None and sum(len(b) for b in self.bands) > 4096:
            self.bands = None
