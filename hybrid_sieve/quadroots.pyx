# cython: language_level=3, boundscheck=False, wraparound=False, cdivision=True
import array

ctypedef unsigned long long u64


cdef inline u64 powmod(u64 a, u64 e, u64 m) noexcept nogil:
    cdef u64 r = 1
    a %= m
    while e:
        if e & 1:
            r = r * a % m
        a = a * a % m
        e >>= 1
    return r


cdef int jacobi(long long a, long long n) noexcept nogil:
    # n odd and positive
    cdef int t = 1
    cdef long long tmp
    a %= n
    while a:
        while (a & 1) == 0:
            a >>= 1
            if (n & 7) == 3 or (n & 7) == 5:
                t = -t
        tmp = a; a = n; n = tmp
        if (a & 3) == 3 and (n & 3) == 3:
            t = -t
        a %= n
    return t if n == 1 else 0


cdef u64 sqrtmod(u64 a, u64 p, u64 q, int s, u64 c) noexcept nogil:
    # Tonelli-Shanks. a is a residue, p - 1 = q * 2^s, c = g^q for a non-residue g.
    # Returns 0 if no root is found (a not a residue, or p not prime).
    cdef u64 r, t, t2, b
    cdef int m = s, i, k
    if s == 1:
        r = powmod(a, (p + 1) >> 2, p)
        return r if r * r % p == a else 0
    r = powmod(a, (q + 1) >> 1, p)
    t = powmod(a, q, p)
    while t != 1:
        i = 0
        t2 = t
        while t2 != 1:
            t2 = t2 * t2 % p
            i += 1
            if i == m:
                return 0
        b = c
        for k in range(m - i - 1):
            b = b * b % p
        r = r * b % p
        c = b * b % p
        t = t * c % p
        m = i
    return r if r * r % p == a else 0


cdef class QuadRoots:
    """Roots for every usable quad without a square root per (prime, quad).

    quads : the k in 1..k_max that are squarefree and have no prime factor
            above `bound`, ascending.
    roots(qi, out) : for quad = quads[qi], sets out[i] to an x with
            quad * x^2 = hmap[i][1]^2 / hmap[i][0]  (mod primes[i]),
            which is n when hmap[i] = (k, sqrt(n*k)); or -1 where no root exists.
    factors(qi), factors_of(k) : the prime factors of a kept quad.
    root(k, i) : the same single value for quad = k and primes[i]; also -1 when
            k is not in quads.
    Primes must be below 2^31. Memory: 4 bytes * len(primes) * pi(bound)."""
    cdef Py_ssize_t nfb, nq, k_max
    cdef int[::1] primes, g, base, ck, tbl, off, fac, qidx
    cdef readonly object quads, small

    def __init__(self, primes, hmap, Py_ssize_t bound, Py_ssize_t k_max):
        self._filter(bound, k_max)
        self._table(primes, hmap)

    cdef _filter(self, Py_ssize_t bound, Py_ssize_t k_max):
        # smallest-prime-factor sieve, then keep squarefree bound-smooth k
        cdef Py_ssize_t d, m, k, x, l, last, nq = 0, nf = 0, cnt, ns = 0
        cdef bint ok
        spf_a = array.array('i', range(k_max + 1))
        cdef int[::1] spf = spf_a
        d = 2
        while d * d <= k_max:
            if spf[d] == d:
                m = d * d
                while m <= k_max:
                    if spf[m] == m:
                        spf[m] = d
                    m += d
            d += 1
        pidx_a = array.array('i', [0]) * (k_max + 1)   # prime -> index in small
        cdef int[::1] pidx = pidx_a
        small = array.array('i')
        for l in range(2, min(bound, k_max) + 1):
            if spf[l] == l:
                pidx[l] = ns
                ns += 1
                small.append(l)
        self.small = small
        mark_a = array.array('b', [0]) * (k_max + 1)
        cdef signed char[::1] mark = mark_a
        for k in range(1, k_max + 1):
            x = k; last = 0; ok = True; cnt = 0
            while x > 1:
                l = spf[x]
                if l > bound or l == last:
                    ok = False
                    break
                last = l
                x //= l
                cnt += 1
            if ok:
                mark[k] = 1
                nq += 1
                nf += cnt
        quads = array.array('i', [0]) * nq
        off_a = array.array('i', [0]) * (nq + 1)
        fac_a = array.array('i', [0]) * max(nf, 1)
        qidx_a = array.array('i', [-1]) * (k_max + 1)   # k -> index in quads, -1 if dropped
        cdef int[::1] qv = quads, off = off_a, fac = fac_a, qidx = qidx_a
        nq = 0; nf = 0
        for k in range(1, k_max + 1):
            if mark[k]:
                qidx[k] = nq
                qv[nq] = k
                off[nq] = nf
                x = k
                while x > 1:
                    l = spf[x]
                    fac[nf] = pidx[l]
                    nf += 1
                    x //= l
                nq += 1
        off[nq] = nf
        self.nq = nq
        self.k_max = k_max
        self.qidx = qidx_a
        self.quads = quads
        self.off = off_a
        self.fac = fac_a

    cdef _table(self, primes, hmap):
        cdef Py_ssize_t nfb = len(primes), ns = len(self.small), i, f
        cdef u64 p, gg, q, c, a, r, k0, x0, l
        cdef int s, ch
        if len(hmap) != nfb:
            raise ValueError("hmap and primes differ in length")
        self.nfb = nfb
        self.primes = array.array('i', primes)
        self.g = array.array('i', [0]) * nfb
        self.base = array.array('i', [0]) * nfb
        self.ck = array.array('i', [0]) * nfb
        # tbl[f*nfb + i] = +1/sqrt(l) if l is a residue mod p, -1/sqrt(l*g) if not, 0 if p | l
        self.tbl = array.array('i', [0]) * max(ns * nfb, 1)
        cdef int[::1] small = self.small, pr = self.primes, gv = self.g
        cdef int[::1] base = self.base, ck = self.ck, tbl = self.tbl
        for i in range(nfb):
            if pr[i] < 2:
                raise ValueError("bad modulus %d" % pr[i])
            p = pr[i]
            k0 = hmap[i][0]
            x0 = hmap[i][1]
            if p == 2 or k0 % p == 0:
                continue                           # ck stays 0: always skipped
            gg = 2
            while gg < p and jacobi(gg, p) != -1:
                gg += 1
            if (p & 1) == 0 or gg == p:
                raise ValueError("%d is not an odd prime" % p)
            q = p - 1
            s = 0
            while (q & 1) == 0:
                q >>= 1
                s += 1
            c = powmod(gg, q, p)
            ch = jacobi(k0, p)
            a = k0 % p if ch == 1 else k0 % p * gg % p
            r = sqrtmod(a, p, q, s, c)
            if r == 0:
                raise ValueError("%d is not an odd prime" % p)
            gv[i] = <int>gg
            ck[i] = ch
            # x0^2 = n*k0, so x0 * sqrt(k0/quad) / k0 squares to n/quad
            base[i] = <int>(x0 % p * r % p * powmod(k0 % p, p - 2, p) % p)
            for f in range(ns):
                l = small[f] % p
                if l == 0:
                    continue
                ch = jacobi(l, p)
                a = l if ch == 1 else l * gg % p
                r = sqrtmod(a, p, q, s, c)
                if r == 0:
                    raise ValueError("%d is not an odd prime" % p)
                r = powmod(r, p - 2, p)
                if ch == 1:
                    tbl[f * nfb + i] = <int>r
                else:
                    tbl[f * nfb + i] = -<int>r

    cdef int _root(self, Py_ssize_t lo, Py_ssize_t hi, Py_ssize_t i) noexcept:
        cdef Py_ssize_t t, nfb = self.nfb
        cdef int v, nr = 0, c = self.ck[i]
        cdef u64 p, x
        if c == 0:
            return -1
        p = self.primes[i]
        x = self.base[i]
        for t in range(lo, hi):
            v = self.tbl[self.fac[t] * nfb + i]
            if v == 0:
                return -1
            if v < 0:
                v = -v
                nr += 1
            x = x * v % p
        # quad and k must have the same character: nr odd iff c == -1
        if (nr & 1) != (c < 0):
            return -1
        for t in range(nr >> 1):
            x = x * self.g[i] % p
        return <int>x

    def roots(self, Py_ssize_t qi, int[::1] out):
        cdef Py_ssize_t i, lo, hi
        if qi < 0 or qi >= self.nq:
            raise IndexError("quad index out of range")
        if out.shape[0] < self.nfb:
            raise ValueError("out is shorter than the factor base")
        lo = self.off[qi]
        hi = self.off[qi + 1]
        for i in range(self.nfb):
            out[i] = self._root(lo, hi, i)

    def factors(self, Py_ssize_t qi):
        """Prime factors of quads[qi], ascending (empty for 1)."""
        cdef Py_ssize_t t
        if qi < 0 or qi >= self.nq:
            raise IndexError("quad index out of range")
        return [self.small[self.fac[t]] for t in range(self.off[qi], self.off[qi + 1])]

    def factors_of(self, Py_ssize_t k):
        """Prime factors of k, or None if k is not in quads."""
        if k < 1 or k > self.k_max or self.qidx[k] < 0:
            return None
        return self.factors(self.qidx[k])

    def root(self, Py_ssize_t k, Py_ssize_t i):
        cdef Py_ssize_t qi
        if i < 0 or i >= self.nfb:
            raise IndexError("prime index out of range")
        if k < 1 or k > self.k_max:
            return -1
        qi = self.qidx[k]
        if qi < 0:
            return -1
        return self._root(self.off[qi], self.off[qi + 1], i)
