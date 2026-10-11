# cython: language_level=3, boundscheck=False, wraparound=False, cdivision=True, initializedcheck=False
"""
amulsearch -- find a small multiplier m such that c = m * y^2, by Legendre signatures.

For r "signature" primes p_j, sig(x) has bit j set when x is a non-residue mod p_j.
sig(x*y) = sig(x) ^ sig(y) and sig(y^2) = 0, so c = m*y^2 forces sig(c) == sig(m).
All usable m (squarefree, |m| < limit, every prime factor in the factor base, both signs)
are stored once in a hash table keyed by sig(m). A candidate c then costs r table lookups
to get sig(c) and one hash probe to get every m that can work.

    s = AMulSearch(sig_primes, fbase, limit)
    s.scan(c0, c1, c2, count)  ->  [(i, m), ...] for c(i) = c2*i^2 + c1*i + c0, 0 <= i < count
    s.lookup(c)                ->  [m, ...] for one integer c

Matches are candidates: confirm with c % m == 0 and c // m a perfect square.
A c that is 0 mod a signature prime is skipped (its symbol there says nothing).
"""
from libc.stdlib cimport malloc, calloc, free, realloc
from libc.string cimport memset
from libc.stdint cimport uint64_t, uint8_t, int32_t

cdef class AMulSearch:
    cdef int r
    cdef int* p            # signature primes
    cdef uint8_t** chi     # chi[j][x] = 1 non-residue, 0 residue, 2 for x == 0
    cdef uint64_t* keys
    cdef int32_t* vals     # 0 = empty slot
    cdef uint64_t hmask
    cdef int hshift
    cdef readonly Py_ssize_t entries, limit
    cdef readonly object sig_primes

    def __cinit__(self, sig_primes, fbase, Py_ssize_t limit):
        cdef Py_ssize_t j, x, q, m, i, L = limit
        cdef int r = len(sig_primes), pj
        cdef uint64_t s, h
        if r < 1 or r > 63:
            raise ValueError("need 1..63 signature primes")
        self.r = r
        self.limit = L
        self.sig_primes = list(sig_primes)
        self.p = <int*> malloc(r * sizeof(int))
        self.chi = <uint8_t**> calloc(r, sizeof(uint8_t*))
        for j in range(r):
            pj = sig_primes[j]
            if pj < 3 or pj > 60000:
                raise ValueError("signature primes must be odd and below 60000")
            self.p[j] = pj
            self.chi[j] = <uint8_t*> malloc(pj)
            memset(self.chi[j], 1, pj)
            for x in range(1, pj):
                self.chi[j][(x * x) % pj] = 0
            self.chi[j][0] = 2
        # smallest prime factor sieve
        cdef int32_t* spf = <int32_t*> malloc(L * sizeof(int32_t))
        cdef uint64_t* sg = <uint64_t*> calloc(L, sizeof(uint64_t))
        cdef uint8_t* ok = <uint8_t*> calloc(L, 1)
        for x in range(L):
            spf[x] = <int32_t> x
        x = 2
        while x * x < L:
            if spf[x] == x:
                m = x * x
                while m < L:
                    if spf[m] == m:
                        spf[m] = <int32_t> x
                    m += x
            x += 1
        allowed = set(fbase) | {2}
        allowed -= set(sig_primes)
        cdef uint8_t* isok = <uint8_t*> calloc(L, 1)
        for q in allowed:
            if 1 < q < L:
                isok[q] = 1
        cdef uint64_t sneg = 0
        for j in range(r):
            if self.p[j] % 4 == 3:
                sneg |= (<uint64_t> 1) << j
        cdef Py_ssize_t cnt = 0
        ok[1] = 1
        cnt = 1
        for m in range(2, L):
            q = spf[m]
            if q == m:
                if isok[m]:
                    s = 0
                    for j in range(r):
                        if self.chi[j][m % self.p[j]] == 1:
                            s |= (<uint64_t> 1) << j
                    sg[m] = s
                    ok[m] = 1
                    cnt += 1
            else:
                x = m // q
                if ok[x] and isok[q] and x % q != 0:
                    sg[m] = sg[q] ^ sg[x]
                    ok[m] = 1
                    cnt += 1
        self.entries = 2 * cnt
        cdef int bits = 4
        while ((<Py_ssize_t> 1) << bits) < 4 * self.entries:
            bits += 1
        self.hshift = 64 - bits
        self.hmask = ((<uint64_t> 1) << bits) - 1
        self.keys = <uint64_t*> calloc((<size_t> 1) << bits, sizeof(uint64_t))
        self.vals = <int32_t*> calloc((<size_t> 1) << bits, sizeof(int32_t))
        if self.keys == NULL or self.vals == NULL:
            raise MemoryError()
        for m in range(1, L):
            if ok[m]:
                self._put(sg[m], <int32_t> m)
                self._put(sg[m] ^ sneg, -(<int32_t> m))
        free(spf); free(sg); free(ok); free(isok)

    def __dealloc__(self):
        cdef int j
        if self.chi != NULL:
            for j in range(self.r):
                free(self.chi[j])
        free(self.chi); free(self.p); free(self.keys); free(self.vals)

    cdef inline void _put(self, uint64_t s, int32_t m) noexcept nogil:
        cdef uint64_t h = (s * <uint64_t> 0x9E3779B97F4A7C15) >> self.hshift
        while self.vals[h] != 0:
            h = (h + 1) & self.hmask
        self.keys[h] = s
        self.vals[h] = m

    def scan(self, c0, c1, c2, Py_ssize_t count):
        """Candidates (i, m) with sig(c2*i^2 + c1*i + c0) == sig(m), 0 <= i < count.

        The symbol of c(i) mod p repeats with period p, so each prime gives one bit plane
        (a p-bit pattern tiled over the interval). The planes are then turned into one
        signature per position 64 positions at a time with a 64x64 bit transpose."""
        cdef Py_ssize_t nw = (count + 63) >> 6, w, x, nb, o, kk, pos, i
        cdef int j, pj, rr, d1, d2, r = self.r, sh, s64, jj, k
        cdef uint8_t* t
        cdef uint8_t v
        cdef uint64_t s, h, msk, tt, bad
        cdef uint64_t A[64]
        cdef uint64_t* pl = <uint64_t*> malloc((r + 1) * nw * sizeof(uint64_t))
        cdef uint64_t* mb = <uint64_t*> malloc((60000 // 64 + 4) * 2 * sizeof(uint64_t))
        cdef uint64_t* zb = mb + (60000 // 64 + 4)
        cdef int* r0 = <int*> malloc(r * 3 * sizeof(int))
        cdef uint64_t* row
        cdef uint64_t* badrow = pl + r * nw
        if pl == NULL or mb == NULL or r0 == NULL:
            free(pl); free(mb); free(r0)
            raise MemoryError()
        for j in range(r):
            pj = self.p[j]
            r0[3 * j] = c0 % pj
            r0[3 * j + 1] = (c1 + c2) % pj
            r0[3 * j + 2] = (2 * c2) % pj
        out = []
        cdef list hits_i = [], hits_m = []
        cdef Py_ssize_t ncand = 0, cap = 1024
        cdef Py_ssize_t* ci = <Py_ssize_t*> malloc(cap * sizeof(Py_ssize_t))
        cdef int32_t* cm = <int32_t*> malloc(cap * sizeof(int32_t))
        with nogil:
            memset(badrow, 0, nw * sizeof(uint64_t))
            for j in range(r):
                pj = self.p[j]
                t = self.chi[j]
                nb = pj + 64
                memset(mb, 0, ((nb >> 6) + 2) * sizeof(uint64_t))
                memset(zb, 0, ((nb >> 6) + 2) * sizeof(uint64_t))
                rr = r0[3 * j]; d1 = r0[3 * j + 1]; d2 = r0[3 * j + 2]
                bad = 0
                for x in range(pj):                      # one period of the symbol of c(i)
                    v = t[rr]
                    if v == 1:
                        o = x
                        while o < nb:
                            mb[o >> 6] |= (<uint64_t> 1) << (o & 63)
                            o += pj
                    elif v == 2:
                        bad = 1
                        o = x
                        while o < nb:
                            zb[o >> 6] |= (<uint64_t> 1) << (o & 63)
                            o += pj
                    rr += d1
                    if rr >= pj:
                        rr -= pj
                    d1 += d2
                    if d1 >= pj:
                        d1 -= pj
                row = pl + j * nw
                s64 = 64 % pj
                o = 0
                for w in range(nw):                      # word w starts at bit 64*w mod p
                    kk = o >> 6
                    sh = o & 63
                    row[w] = (mb[kk] >> sh) | ((mb[kk + 1] << 1) << (63 - sh))
                    if bad:
                        badrow[w] |= (zb[kk] >> sh) | ((zb[kk + 1] << 1) << (63 - sh))
                    o += s64
                    if o >= pj:
                        o -= pj
            for w in range(nw):
                for j in range(r):
                    A[j] = pl[j * nw + w]
                for j in range(r, 64):
                    A[j] = 0
                # 64x64 bit transpose: afterwards A[pos] bit j == plane j bit pos
                jj = 32
                msk = <uint64_t> 0x00000000FFFFFFFF
                while jj:
                    k = 0
                    while k < 64:
                        tt = ((A[k] >> jj) ^ A[k | jj]) & msk
                        A[k | jj] ^= tt
                        A[k] ^= tt << jj
                        k = ((k | jj) + 1) & ~jj
                    jj >>= 1
                    msk ^= msk << jj
                bad = badrow[w]
                for pos in range(64):
                    i = (w << 6) + pos
                    if i >= count:
                        break
                    if (bad >> pos) & 1:
                        continue
                    s = A[pos]
                    h = (s * <uint64_t> 0x9E3779B97F4A7C15) >> self.hshift
                    while self.vals[h] != 0:
                        if self.keys[h] == s:
                            if ncand == cap:
                                cap *= 2
                                ci = <Py_ssize_t*> realloc(ci, cap * sizeof(Py_ssize_t))
                                cm = <int32_t*> realloc(cm, cap * sizeof(int32_t))
                            ci[ncand] = i
                            cm[ncand] = self.vals[h]
                            ncand += 1
                        h = (h + 1) & self.hmask
        for x in range(ncand):
            out.append((ci[x], cm[x]))
        free(pl); free(mb); free(r0); free(ci); free(cm)
        return out

    def scan_slow(self, c0, c1, c2, Py_ssize_t count):
        """Candidates (i, m) with sig(c2*i^2 + c1*i + c0) == sig(m), 0 <= i < count."""
        cdef Py_ssize_t i
        cdef int j, pj, rr, d1, d2
        cdef uint8_t* t
        cdef uint8_t v
        cdef uint64_t bit, s, h
        cdef uint64_t* sg = <uint64_t*> calloc(count, sizeof(uint64_t))
        if sg == NULL:
            raise MemoryError()
        cdef int* r0 = <int*> malloc(self.r * 3 * sizeof(int))
        for j in range(self.r):
            pj = self.p[j]
            r0[3 * j] = c0 % pj
            r0[3 * j + 1] = (c1 + c2) % pj
            r0[3 * j + 2] = (2 * c2) % pj
        out = []
        cdef uint64_t BAD = (<uint64_t> 1) << 63
        with nogil:
            for j in range(self.r):
                pj = self.p[j]
                t = self.chi[j]
                bit = (<uint64_t> 1) << j
                rr = r0[3 * j]; d1 = r0[3 * j + 1]; d2 = r0[3 * j + 2]
                for i in range(count):
                    v = t[rr]
                    if v == 1:
                        sg[i] |= bit
                    elif v == 2:
                        sg[i] |= BAD
                    rr += d1
                    if rr >= pj:
                        rr -= pj
                    d1 += d2
                    if d1 >= pj:
                        d1 -= pj
        for i in range(count):
            s = sg[i]
            if s & BAD:
                continue
            h = (s * <uint64_t> 0x9E3779B97F4A7C15) >> self.hshift
            while self.vals[h] != 0:
                if self.keys[h] == s:
                    out.append((i, self.vals[h]))
                h = (h + 1) & self.hmask
        free(sg); free(r0)
        return out

    def lookup(self, c):
        """Every stored m with sig(m) == sig(c); [] when c is 0 mod a signature prime."""
        cdef uint64_t s = 0, h
        cdef int j, v
        for j in range(self.r):
            v = self.chi[j][c % self.p[j]]
            if v == 2:
                return []
            if v == 1:
                s |= (<uint64_t> 1) << j
        out = []
        h = (s * <uint64_t> 0x9E3779B97F4A7C15) >> self.hshift
        while self.vals[h] != 0:
            if self.keys[h] == s:
                out.append(self.vals[h])
            h = (h + 1) & self.hmask
        return out
