###DISCLAIMER: I have generated this code with claude because it is well documented sieving logic and it saves me time.

# cython: language_level=3, boundscheck=False, wraparound=False, cdivision=True
import array

cdef extern from *:
    int __builtin_ctzll(unsigned long long) nogil   # GCC/Clang

ctypedef unsigned long long u64


def new_interval(Py_ssize_t n):
    """n positions, all set to 1, packed 64 per word."""
    cdef Py_ssize_t nw = (n + 63) >> 6
    a = array.array('Q', [0xFFFFFFFFFFFFFFFF]) * nw
    if n & 63:
        a[nw - 1] = ((<u64>1) << (n & 63)) - 1      # clear bits past n
    return a


def make_pattern(Py_ssize_t p, sols, step=1):
    """p words; bit i is set iff (step*i % p) is in sols.
    p must be odd and coprime to step."""
    cdef Py_ssize_t nbits = p + 64, r, i, j, k, o = 0, s64 = 64 % p
    cdef Py_ssize_t inv = pow(step, -1, p)
    cdef int s
    base = array.array('Q', [0]) * ((nbits >> 6) + 2)
    pat = array.array('Q', [0]) * p
    cdef u64[::1] m = base, w = pat
    for r in sols:                                  # one period plus 64 bits
        i = r % p
        if i < 0:
            i += p
        i = (i * inv) % p                           # residue r -> position r/step
        while i < nbits:
            m[i >> 6] |= (<u64>1) << (i & 63)
            i += p
    for j in range(p):                              # word j starts at bit 64*j mod p
        k = o >> 6
        s = o & 63
        w[j] = (m[k] >> s) | ((m[k + 1] << 1) << (63 - s))
        o += s64
        if o >= p:
            o -= p
    return pat


def sieve(u64[::1] iv, u64[::1] pat, root=0, step=1):
    """Clear every position i where (root + step*i) % p is not a solution.
    pat must come from make_pattern(p, sols, step) with the same step.
    root and step may be arbitrarily large Python ints."""
    cdef Py_ssize_t p = pat.shape[0], nw = iv.shape[0], j = 0, i, c
    # start word k:  64*k = root/step (mod p)
    cdef Py_ssize_t k = (root * pow(step, -1, p) * pow(64, -1, p)) % p
    while j < nw:
        c = p - k
        if c > nw - j:
            c = nw - j
        for i in range(c):
            iv[j + i] &= pat[k + i]
        j += c
        k = 0


def survivors(u64[::1] iv):
    """Positions still set to 1, ascending."""
    cdef Py_ssize_t j
    cdef u64 w
    out = []
    for j in range(iv.shape[0]):
        w = iv[j]
        while w:
            out.append((j << 6) + __builtin_ctzll(w))
            w &= w - 1
    return out
