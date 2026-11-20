###DISCLAIMER: I have generated this code with claude because it is well documented linear algebra logic and it saves me time.


# cython: language_level=3, boundscheck=False, wraparound=False, cdivision=True, initializedcheck=False
"""
almostsq -- products of B-smooth numbers with few odd-exponent primes.

GF(2) elimination and the single/pair search run in C on packed 64-bit words.

    import almostsq

    # from the numbers themselves (factored with sympy)
    for odd_primes, used in almostsq.almost_squares(smooths, keep=10): ...

    # from parity vectors you already have: masks[i] bit j set  <=>  prime j
    # has odd exponent in relation i; m = number of primes
    for prime_mask, history_mask in almostsq.from_masks(masks, m, keep=10): ...

history_mask bit i set  <=>  relation i is part of the product.
Assumes a little-endian host (x86, ARM).
"""
from libc.stdlib cimport malloc, calloc, free
from libc.string cimport memcpy
from libc.stdint cimport uint64_t
from math import isqrt, prod

cdef extern int __builtin_popcountll(unsigned long long) nogil


cdef inline void _insert(int* cw, Py_ssize_t* ci, Py_ssize_t* cj,
                         Py_ssize_t* count, Py_ssize_t keep,
                         int w, Py_ssize_t i, Py_ssize_t j) noexcept nogil:
    """Keep the `keep` lightest candidates, sorted, earlier ones first on ties."""
    cdef Py_ssize_t p
    if count[0] == keep:
        if w >= cw[keep - 1]:
            return
        p = keep - 1
    else:
        p = count[0]
        count[0] += 1
    while p > 0 and cw[p - 1] > w:
        cw[p] = cw[p - 1]; ci[p] = ci[p - 1]; cj[p] = cj[p - 1]
        p -= 1
    cw[p] = w; ci[p] = i; cj[p] = j


cdef class Reduced:
    """Relation matrix [M | I] over GF(2), packed; reduce() then search()."""
    cdef uint64_t* P          # n x wp : prime parities
    cdef uint64_t* H          # n x wh : history (which relations were combined)
    cdef readonly Py_ssize_t n, m
    cdef Py_ssize_t wp, wh

    def __cinit__(self, masks, Py_ssize_t m):
        cdef Py_ssize_t i, n = len(masks)
        cdef bytes b
        if m < 0:
            raise ValueError("m must be >= 0")
        self.n = n
        self.m = m
        self.wp = max(1, (m + 63) // 64)
        self.wh = max(1, (n + 63) // 64)
        self.P = <uint64_t*> calloc(max(n, 1) * self.wp, 8)
        self.H = <uint64_t*> calloc(max(n, 1) * self.wh, 8)
        if self.P == NULL or self.H == NULL:
            raise MemoryError()
        for i in range(n):
            mask = masks[i]
            if mask < 0 or mask >> m:
                raise ValueError(f"mask {i} has bits outside 0..{m - 1}")
            b = mask.to_bytes(self.wp * 8, "little")
            memcpy(self.P + i * self.wp, <const char*> b, self.wp * 8)
            self.H[i * self.wh + i // 64] = (<uint64_t> 1) << (i % 64)

    def __dealloc__(self):
        free(self.P)
        free(self.H)

    cpdef reduce(self):
        """Full reduction: every pivot column ends up in exactly one row."""
        cdef Py_ssize_t n = self.n, wp = self.wp, wh = self.wh
        cdef Py_ssize_t cur, i, k, pw
        cdef uint64_t bit
        cdef uint64_t* rp
        cdef uint64_t* rh
        cdef uint64_t* tp
        cdef uint64_t* th
        cdef bint in_p
        with nogil:
            for cur in range(n):
                rp = self.P + cur * wp
                rh = self.H + cur * wh
                in_p = False
                pw = 0
                bit = 0
                for k in range(wp):
                    if rp[k]:
                        pw = k; bit = rp[k] & (~rp[k] + 1); in_p = True
                        break
                if not in_p:                # prime part cancelled: pivot on
                    for k in range(wh):     # history (never all-zero)
                        if rh[k]:
                            pw = k; bit = rh[k] & (~rh[k] + 1)
                            break
                for i in range(n):
                    if i == cur:
                        continue
                    tp = self.P + i * wp
                    th = self.H + i * wh
                    if in_p:
                        if tp[pw] & bit:
                            for k in range(pw, wp):
                                tp[k] ^= rp[k]
                            for k in range(wh):
                                th[k] ^= rh[k]
                    elif th[pw] & bit:
                        for k in range(pw, wh):
                            th[k] ^= rh[k]

    cdef object _pm(self, Py_ssize_t i):
        return int.from_bytes((<char*> (self.P + i * self.wp))[:self.wp * 8], "little")

    cdef object _hm(self, Py_ssize_t i):
        return int.from_bytes((<char*> (self.H + i * self.wh))[:self.wh * 8], "little")

    def rows(self):
        """All rows as (prime_mask, history_mask) Python ints."""
        return [(self._pm(i), self._hm(i)) for i in range(self.n)]

    def search(self, Py_ssize_t keep=10):
        """The `keep` lightest single rows and pairs of rows, lightest first.

        Pairs involving an exact-square row are skipped (they only repeat the
        other row's odd primes).
        """
        cdef Py_ssize_t n = self.n, wp = self.wp
        cdef Py_ssize_t i, j, k, count = 0
        cdef int c, bound, d
        cdef uint64_t* a
        cdef uint64_t* b
        if keep <= 0 or n == 0:
            return []
        cdef int* w = <int*> malloc(n * sizeof(int))
        cdef int* cw = <int*> malloc(keep * sizeof(int))
        cdef Py_ssize_t* ci = <Py_ssize_t*> malloc(keep * sizeof(Py_ssize_t))
        cdef Py_ssize_t* cj = <Py_ssize_t*> malloc(keep * sizeof(Py_ssize_t))
        if w == NULL or cw == NULL or ci == NULL or cj == NULL:
            free(w); free(cw); free(ci); free(cj)
            raise MemoryError()
        try:
            with nogil:
                for i in range(n):
                    a = self.P + i * wp
                    c = 0
                    for k in range(wp):
                        c += __builtin_popcountll(a[k])
                    w[i] = c
                    _insert(cw, ci, cj, &count, keep, c, i, -1)
                for i in range(n):
                    if w[i] == 0:
                        continue
                    a = self.P + i * wp
                    for j in range(i + 1, n):
                        if w[j] == 0:
                            continue
                        bound = cw[keep - 1] if count == keep else 0x7fffffff
                        d = w[i] - w[j]
                        if d < 0:
                            d = -d
                        if d >= bound:          # weight(a^b) >= |w_a - w_b|
                            continue
                        b = self.P + j * wp
                        c = 0
                        for k in range(wp):
                            c += __builtin_popcountll(a[k] ^ b[k])
                            if c >= bound:
                                break
                        if c < bound:
                            _insert(cw, ci, cj, &count, keep, c, i, j)
            out = []
            for k in range(count):
                pm = self._pm(ci[k]); hm = self._hm(ci[k])
                if cj[k] >= 0:
                    pm ^= self._pm(cj[k]); hm ^= self._hm(cj[k])
                out.append((pm, hm))
            return out
        finally:
            free(w); free(cw); free(ci); free(cj)


def reduce_with_history(masks, Py_ssize_t m):
    """Reduced rows of [M | I] as a list of (prime_mask, history_mask)."""
    cdef Reduced r = Reduced(masks, m)
    r.reduce()
    return r.rows()


def from_masks(masks, Py_ssize_t m, Py_ssize_t keep=10):
    """Lightest combinations for relations given as parity bitmasks."""
    cdef Reduced r = Reduced(masks, m)
    r.reduce()
    return r.search(keep)


def parity_masks(smooths):
    """Factor with sympy. Returns (masks, primes); -1 counts as a prime."""
    from sympy import factorint
    facs = [factorint(x) for x in smooths]
    primes = sorted({p for f in facs for p in f})
    col = {p: j for j, p in enumerate(primes)}
    masks = []
    for f in facs:
        v = 0
        for p, e in f.items():
            if e & 1:
                v |= 1 << col[p]
        masks.append(v)
    return masks, primes


def bits(mask):
    """Indices of the set bits of a Python int."""
    out = []
    i = 0
    while mask:
        if mask & 1:
            out.append(i)
        mask >>= 1
        i += 1
    return out


def almost_squares(smooths, Py_ssize_t keep=10):
    """Return [(odd_primes, used_smooths), ...], fewest odd primes first.

    Zeros are rejected; 1s and duplicates are dropped. Every result is checked:
    product(used) / product(odd_primes) must be a perfect square.
    """
    if 0 in smooths:
        raise ValueError("0 is not a valid smooth number")
    smooths = [x for x in dict.fromkeys(smooths) if x != 1]
    if not smooths:
        return []
    masks, primes = parity_masks(smooths)
    out = []
    for pm, hm in from_masks(masks, len(primes), keep):
        used = [smooths[i] for i in bits(hm)]
        odd = [primes[j] for j in bits(pm)]
        q, r = divmod(prod(used), prod(odd))
        if r or q < 0 or isqrt(q) ** 2 != q:
            raise RuntimeError(f"self-check failed for {used}")
        out.append((odd, used))
    return out
