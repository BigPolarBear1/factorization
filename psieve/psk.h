/* psk.h -- the inner search of psieve for one leftover a:  X^2 - a*y^2 = 4*n*K  over a list of multipliers K.
 *
 * For each K the classes X = +-sqrt(4nK) mod each prime of a are walked in Gray-code order. Each class is a
 * window of positions X = Xstart + (a*step)*i, i = 0..window-1, and a position survives when (X^2-4nK)/a is a
 * square mod 64 and mod every sieve prime. Survivors are returned; the caller checks them exactly.
 *
 * Everything that is "mod a sieve prime" is done for all sieve primes at once in short fixed-length loops over
 * 16-bit lanes (the compiler turns them into vector instructions):
 *   Xstart = -sum_j rho_j*(|a|/q_j) + Q*|a|     (rho_j = root_j/(|a|/q_j) mod q_j, Q a known integer)
 * so the pattern index of a class, Xstart/(a*step*64) mod p, needs only rho_j mod p and Q mod p, and a Gray
 * step changes it by +-2*rho_j*(|a|/q_j) plus a multiple of |a|.
 *
 * Square roots of 4nK mod the primes q_j of a are multiplied together from the square roots of the prime factors
 * of K (one table per leftover), so no modular square root is taken per multiplier.
 */
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <immintrin.h>
#include <x86intrin.h>

#define PSK_NP 32
#define PSK_MAXW 12
#define PSK_NW 64
#ifndef PSK_PROF
#define PSK_PROF 0          /* 1: count cycles per phase into PskCtx.tsc (PSieveCore.prof()) */
#endif

#if defined(__AVX512F__) && defined(__AVX512BW__) && defined(__AVX512VL__)
#define PSK_WIDE 1
#else
#define PSK_WIDE 0
#endif
#if defined(__AVX512F__)
#define PSK_BW 8
#define PSK_BLK __m512i
#define PSK_LOAD(p) _mm512_loadu_si512((const void *)(p))
#define PSK_STORE(p, x) _mm512_store_si512((void *)(p), x)
#define PSK_AND(a, b) _mm512_and_si512(a, b)
#define PSK_SET1(x) _mm512_set1_epi64((long long)(x))
#define PSK_NZ(x) _mm512_test_epi64_mask(x, x)
#elif defined(__AVX2__)
#define PSK_BW 4
#define PSK_BLK __m256i
#define PSK_LOAD(p) _mm256_loadu_si256((const void *)(p))
#define PSK_STORE(p, x) _mm256_store_si256((void *)(p), x)
#define PSK_AND(a, b) _mm256_and_si256(a, b)
#define PSK_SET1(x) _mm256_set1_epi64x((long long)(x))
#define PSK_NZ(x) (0xF ^ _mm256_movemask_pd(_mm256_castsi256_pd(_mm256_cmpeq_epi64(x, _mm256_setzero_si256()))))
#else
#define PSK_BW 1
#define PSK_BLK uint64_t
#define PSK_LOAD(p) (*(const uint64_t *)(p))
#define PSK_STORE(p, x) (*(uint64_t *)(p) = (x))
#define PSK_AND(a, b) ((a) & (b))
#define PSK_SET1(x) ((uint64_t)(x))
#define PSK_NZ(x) ((x) != 0)
#endif

typedef unsigned __int128 psk_u128;
typedef __int128 psk_i128;
typedef uint16_t psk_v[PSK_NP];          /* inside PskCtx: not aligned (the struct lives in a Python object) */
typedef uint16_t psk_va[PSK_NP] __attribute__((aligned(64)));

typedef struct {
    int w, nsp, sgn;
    psk_u128 ua;                              /* |a| */
    long double uaf;
    uint64_t aq[PSK_MAXW];                    /* primes q_j of a */
    uint64_t aqbm[PSK_MAXW];                  /* (2^64-1)/q */
    uint64_t cj[PSK_MAXW];                    /* (|a|/q)^-1 mod q */
    psk_u128 Aj[PSK_MAXW];                    /* |a|/q */
    psk_v p, m16, c16, c32, c48, one;         /* sieve primes, floor(65536/p), 2^16, 2^32, 2^48 mod p, 1 */
    psk_v beta[3], beta16[3];                 /* |a|/(a*step*64) mod p, and 65536 times it */
    psk_v nb[3][3];                           /* p - e*beta mod p, e = 0, 1, 2 (in 1..p) */
    psk_v alpha[3][PSK_MAXW];                 /* beta/q_j mod p */
    psk_v alpha16[3][PSK_MAXW];               /* 65536*beta/q_j mod p */
    psk_v tcoef[3];                           /* 4n/(a*step)^2 mod p: pattern row of K is K*tcoef mod p */
    psk_v kg[3][4];                           /* shift of the pattern index when Xstart moves by g*|a| */
    uint64_t **pat[PSK_NP];                   /* pattern rows of each sieve prime */
    int stp_of[16], res_of[16];
    const uint64_t *mask2;                    /* [16][64] */
    psk_u128 n4lo;
    long double n4f;
    long long kfast;                          /* K up to here: sqrt(4nK) from one long double square root */
    long long kwide;                          /* K up to here: the same with a 256-bit remainder (4nK < 2^248) */
    uint64_t n4hi;                            /* bits 128.. of 4n */
    /* square roots for coded multipliers */
    const uint32_t *sg;                       /* sg[f*w + j]: sqrt(f) or sqrt(f*g_j) mod q_j */
    const uint64_t *nrp;                      /* 4 bits per j: 1 where the g_j form was used */
    const uint16_t *ramj;                     /* j+1 if the prime is q_j, else 0 */
    uint32_t s4n[PSK_MAXW];                   /* (sqrt(4n) or sqrt(4n*g_j)) * cj */
    uint64_t nr4n;
    uint32_t sg2[PSK_MAXW];
    uint64_t nr2;
    uint32_t gip[PSK_MAXW][8];                /* g_j^-t mod q_j */
    long long intervals, multipliers, skipped;
    unsigned long long tsc[4];
} PskCtx;

typedef struct {
    int32_t kidx, wd;
    uint64_t bits, dlo, dhi, crlo, crhi;
    int32_t stp, pad;
} PskHit;

static uint64_t psk_ones[PSK_NW + 160 + 8];
static int psk_ones_ready = 0;

static inline uint64_t psk_bmod(uint64_t x, uint64_t m, uint64_t p) {
    uint64_t r = x - (uint64_t)(((psk_u128)x * m) >> 64) * p;
    return r >= p ? r - p : r;
}

/* x mod p, lane by lane (x < 65536) */
#define PSK_RED(c, x, r) do { int i_; for (i_ = 0; i_ < PSK_NP; i_++) { \
        uint16_t x_ = (x); \
        uint16_t q_ = (uint16_t)(((uint32_t)x_ * (c)->m16[i_]) >> 16); \
        uint16_t t_ = (uint16_t)(x_ - q_ * (c)->p[i_]); \
        uint16_t u_ = (uint16_t)(t_ - (c)->p[i_]); \
        (r)[i_] = t_ < u_ ? t_ : u_; } } while (0)

/* x * coef mod p for all sieve primes, x below 2^32; coef16 = 65536*coef mod p */
static inline void psk_redmul(const PskCtx *c, uint32_t x, const uint16_t *restrict coef, const uint16_t *restrict coef16, uint16_t *restrict r) {
    const uint16_t hi = (uint16_t)(x >> 16), lo = (uint16_t)x;
    psk_va a, b;
    int i;
    PSK_RED(c, hi, a);
    PSK_RED(c, lo, b);
    for (i = 0; i < PSK_NP; i++) a[i] = (uint16_t)(a[i] * coef16[i] + b[i] * coef[i]);      /* < 2*137^2 */
    PSK_RED(c, a[i_], r);
}

/* codes: codes[m] lists the odd prime factors of K = kl[m] (with multiplicity) as 8-bit indices + 1, lowest byte
 * first; 0 = no code, and then roots[] is read: roots[m*w+j] = a square root of 4nK mod q_j (> 0), 0 if q_j divides
 * K, and roots[m*w] = -1 if K is not usable. Either pointer may be NULL if it is never needed.
 * Returns the number of survivors written; *done says how many multipliers were finished. */
static long psk_run(PskCtx *c, const long long *kl, const uint64_t *codes, const int32_t *roots, long nk, int window,
                    PskHit *out, long maxout, long *done) {
    const int w = c->w, nsp = c->nsp, sgn = c->sgn, nw = (window + 63) >> 6;
    const psk_u128 ua = c->ua;
    const int ua4 = (int)((uint64_t)ua & 3);
    long nout = 0, m;
    uint64_t iv[PSK_NW + 8] __attribute__((aligned(64)));
    psk_va V[PSK_MAXW][4];
    psk_va idx, k0v, kp, rowt, hv, tau, tA, tB;
    psk_va perm;
    psk_u128 tm2[PSK_MAXW];
    uint64_t rho[PSK_MAXW];
    int f2[PSK_MAXW], act[PSK_MAXW];
    int l1[PSK_NP + 1];
    const uint64_t *rowp[16] __attribute__((aligned(64)));
    int i;

    if (!psk_ones_ready) {
        memset(psk_ones, 0xff, sizeof(psk_ones));
        psk_ones_ready = 1;
    }
    for (m = 0; m < nk; m++) {
        const long long K = kl[m];
        int stp = c->stp_of[K & 15], si, res, nact = 0, j, jj, n1, ncls, cls, gray, cr64, cr4, kap1 = 0, kap2 = 0;
        psk_u128 Tlo, cr, bt, dbase, c1, crm;
        psk_i128 r;
        long double Tf, qf, crf;
        uint64_t Q1, Qtot, code;
#if PSK_PROF
        unsigned long long t_a = __rdtsc(), t_b, t_c, t_d;
#endif
        if (nout > maxout / 2) break;
        if (stp == 0 || K > c->kwide) { c->skipped++; continue; }
        code = codes ? codes[m] : 0;
        if (code) {
            /* roots from the prime factors of K */
            uint64_t s[PSK_MAXW];
            uint64_t cnt = c->nr4n;
            uint32_t dead = 0;
            int e2 = __builtin_ctzll((unsigned long long)K);
            for (j = 0; j < w; j++) s[j] = c->s4n[j];
            while (e2--) {
                for (j = 0; j < w; j++) s[j] = psk_bmod(s[j] * c->sg2[j], c->aqbm[j], c->aq[j]);
                cnt += c->nr2;
            }
            while (code) {
                const int f = (int)(code & 255) - 1;
                const uint32_t *sf = c->sg + (size_t)f * w;
                code >>= 8;
                for (j = 0; j < w; j++) s[j] = psk_bmod(s[j] * sf[j], c->aqbm[j], c->aq[j]);
                cnt += c->nrp[f];
                if (c->ramj[f]) dead |= 1u << (c->ramj[f] - 1);
            }
            if (cnt & 0x1111111111111111ULL) {
                /* an odd count at a live prime: 4nK is not a square there */
                const uint64_t odd = cnt & 0x1111111111111111ULL;
                int bad = 0;
                for (j = 0; j < w; j++) if (((odd >> (4 * j)) & 1) && !((dead >> j) & 1)) bad = 1;
                if (bad) { c->skipped++; continue; }
            }
            for (j = 0; j < w; j++) {
                if ((dead >> j) & 1) continue;
                rho[nact] = psk_bmod(s[j] * c->gip[j][(cnt >> (4 * j + 1)) & 7], c->aqbm[j], c->aq[j]);
                act[nact++] = j;
            }
        } else {
            const int32_t *sr = roots + m * w;
            if (sr[0] < 0) { c->skipped++; continue; }
            for (j = 0; j < w; j++) if (sr[j] > 0) {
                rho[nact] = psk_bmod((uint64_t)sr[j] * c->cj[j], c->aqbm[j], c->aq[j]);
                act[nact++] = j;
            }
        }
#if PSK_PROF
        t_b = __rdtsc();
#endif
        /* cr = floor(sqrt(4nK)) */
        Tlo = c->n4lo * (psk_u128)(uint64_t)K;
        Tf = c->n4f * (long double)K;
        crf = sqrtl(Tf);
        if (crf < 9.0e18L) cr = (psk_u128)(uint64_t)crf;
        else {
            const uint64_t ch = (uint64_t)(crf * (1.0L / 18446744073709551616.0L));
            cr = ((psk_u128)ch << 64) | (psk_u128)(uint64_t)(crf - (long double)ch * 18446744073709551616.0L);
        }
        if (K <= c->kfast) {
            /* one Newton correction; the remainder fits 128 bits (signed) */
            long double rf, df;
            r = (psk_i128)(Tlo - cr * cr);
            rf = (long double)(int64_t)(r >> 64) * 18446744073709551616.0L + (long double)(uint64_t)r;
            df = rf / (2.0L * crf);
            if (df >= 0) cr += (psk_u128)(uint64_t)df; else cr -= (psk_u128)(uint64_t)(-df);
        } else {
            /* the same with the remainder T - cr^2 taken in 256 bits */
            const uint64_t k = (uint64_t)K, c0 = (uint64_t)cr, c1_ = (uint64_t)(cr >> 64);
            const psk_u128 t0 = (psk_u128)(uint64_t)c->n4lo * k;
            const psk_u128 t1 = (psk_u128)(uint64_t)(c->n4lo >> 64) * k + (t0 >> 64);
            const psk_u128 Thi = (psk_u128)c->n4hi * k + (t1 >> 64);
            const psk_u128 Tl = ((psk_u128)(uint64_t)t1 << 64) | (uint64_t)t0;
            const psk_u128 s00 = (psk_u128)c0 * c0, s01 = (psk_u128)c0 * c1_, s11 = (psk_u128)c1_ * c1_;
            /* cr^2 = s00 + 2*s01*2^64 + s11*2^128 */
            psk_u128 Pl = s00, Ph = s11, add = s01 << 65, Rl, Rh;
            long double rf, df;
            Pl += add; Ph += (Pl < add) + (s01 >> 63);
            Rl = Tl - Pl; Rh = Thi - Ph - (Tl < Pl);
            rf = ((long double)(int64_t)(Rh >> 64) * 18446744073709551616.0L + (long double)(uint64_t)Rh) * 340282366920938463463374607431768211456.0L
               + ((long double)(uint64_t)(Rl >> 64) * 18446744073709551616.0L + (long double)(uint64_t)Rl);
            df = rf / (2.0L * crf);
            if (df >= 0) cr += (psk_u128)(uint64_t)df; else cr -= (psk_u128)(uint64_t)(-df);
        }
        r = (psk_i128)(Tlo - cr * cr);
        while (r < 0) { cr -= 1; r += (psk_i128)(2 * cr + 1); }
        while (r > (psk_i128)(2 * cr)) { r -= (psk_i128)(2 * cr + 1); cr += 1; }
        /* class offsets */
        bt = 0;
        for (jj = 0; jj < nact; jj++) {
            psk_u128 term = (psk_u128)rho[jj] * c->Aj[act[jj]], T2;
            bt += ua - term;
            if (bt >= ua) { bt -= ua; kap1++; }
            T2 = 2 * term;
            f2[jj] = T2 >= ua;
            tm2[jj] = f2[jj] ? T2 - ua : T2;
        }
        c1 = sgn > 0 ? cr + 1 : cr;
        qf = ((long double)(uint64_t)(c1 >> 64) * 18446744073709551616.0L + (long double)(uint64_t)c1) / c->uaf;
        if (qf >= 4.0e9L) { c->skipped++; continue; }
        Q1 = (uint64_t)qf;
        crm = c1 - (psk_u128)Q1 * ua;
        while ((psk_i128)crm < 0) { crm += ua; Q1--; }
        while (crm >= ua) { crm -= ua; Q1++; }
        if (sgn > 0) {
            if (bt < crm) { dbase = bt + ua - crm; kap2 = 1; } else dbase = bt - crm;
            Qtot = (uint64_t)nact + Q1 + (uint64_t)kap2 - (uint64_t)kap1;
        } else {
            if (crm < bt) { dbase = crm + ua - bt; kap2 = 1; } else dbase = crm - bt;
            Qtot = (uint64_t)nact + Q1 - (uint64_t)kap1 - (uint64_t)kap2;
        }
        if (Qtot >> 32) { c->skipped++; continue; }
        si = stp == 1 ? 0 : (stp == 2 ? 1 : 2);
        res = c->res_of[K & 15];
        cr64 = (int)((uint64_t)cr & 63);
        cr4 = cr64 & 3;
#if PSK_PROF
        t_c = __rdtsc();
#endif
        /* everything mod the sieve primes */
        {
            const uint16_t khi = (uint16_t)((uint64_t)K >> 16), klo = (uint16_t)K;
            psk_va a_, b_;
            if ((uint64_t)K >> 32) {
                psk_redmul(c, (uint32_t)((uint64_t)K >> 32), c->c32, c->c48, a_);
                psk_redmul(c, (uint32_t)K, c->one, c->c16, b_);
                for (i = 0; i < PSK_NP; i++) a_[i] = (uint16_t)(a_[i] + b_[i]);
            } else {
                PSK_RED(c, khi, a_);
                PSK_RED(c, klo, b_);
                for (i = 0; i < PSK_NP; i++) a_[i] = (uint16_t)(a_[i] * c->c16[i] + b_[i]);
            }
            PSK_RED(c, a_[i_], kp);
            for (i = 0; i < PSK_NP; i++) a_[i] = (uint16_t)(kp[i] * c->tcoef[si][i]);
            PSK_RED(c, a_[i_], rowt);
        }
        psk_redmul(c, (uint32_t)Qtot, c->beta[si], c->beta16[si], idx);
        for (jj = 0; jj < nact; jj++) {
            const uint16_t *nbA = c->nb[si][f2[jj]], *nbB = c->nb[si][f2[jj] + 1];
            uint16_t *restrict v0 = V[jj][sgn > 0 ? 0 : 2], *restrict v1 = V[jj][sgn > 0 ? 1 : 3];
            uint16_t *restrict v2 = V[jj][sgn > 0 ? 2 : 0], *restrict v3 = V[jj][sgn > 0 ? 3 : 1];
            psk_redmul(c, (uint32_t)rho[jj], c->alpha[si][act[jj]], c->alpha16[si][act[jj]], hv);
            for (i = 0; i < PSK_NP; i++) {
                const uint16_t p = c->p[i];
                uint16_t x = (uint16_t)(idx[i] + p - hv[i]), y = (uint16_t)(x - p);
                uint16_t t2 = (uint16_t)(2 * hv[i]), u2 = (uint16_t)(t2 - p);
                idx[i] = x < y ? x : y;
                tau[i] = t2 < u2 ? t2 : u2;
            }
            for (i = 0; i < PSK_NP; i++) {
                const uint16_t p = c->p[i];
                uint16_t xa = (uint16_t)(tau[i] + nbA[i]), ya = (uint16_t)(xa - p);
                uint16_t xb = (uint16_t)(tau[i] + nbB[i]), yb = (uint16_t)(xb - p);
                tA[i] = xa < ya ? xa : ya;
                tB[i] = xb < yb ? xb : yb;
            }
            for (i = 0; i < PSK_NP; i++) {
                const uint16_t p = c->p[i];
                uint16_t na = (uint16_t)(p - tA[i]), nna = (uint16_t)(na - p);
                uint16_t nbb = (uint16_t)(p - tB[i]), nnb = (uint16_t)(nbb - p);
                v0[i] = tA[i]; v1[i] = tB[i];
                v2[i] = na < nna ? na : nna;
                v3[i] = nbb < nnb ? nbb : nnb;
            }
        }
        n1 = 0;
        for (i = 0; i < nsp; i++) {                        /* a sieve prime dividing K passes every position */
            l1[n1] = i;
            n1 += kp[i] != 0;
        }
        for (i = 0; i < 16; i++) {
            if (i < n1) rowp[i] = c->pat[l1[i]][rowt[l1[i]]];
            else { rowp[i] = psk_ones; l1[i] = PSK_NP - 1; }
        }
        for (i = 0; i < PSK_NP; i++) perm[i] = (uint16_t)(i < 16 ? l1[i] : 0);
        ncls = 1 << nact;
        cls = 0;
        c->multipliers++;
        c->intervals += ncls;
#if PSK_PROF
        t_d = __rdtsc(); c->tsc[0] += t_b - t_a; c->tsc[1] += t_c - t_b; c->tsc[2] += t_d - t_c;
#endif
        {
#if PSK_WIDE
            const __m256i pv0 = _mm256_loadu_si256((const void *)c->p), pv1 = _mm256_loadu_si256((const void *)(c->p + 16));
            const __m256i permv = _mm256_load_si256((const void *)perm);
            const __m512i rp0 = _mm512_loadu_si512((const void *)rowp), rp1 = _mm512_loadu_si512((const void *)(rowp + 8));
            __m256i ix0 = _mm256_load_si256((const void *)idx), ix1 = _mm256_load_si256((const void *)(idx + 16));
#endif
        for (gray = 0; gray < ncls; gray++) {
            psk_u128 d;
            int g, wd;
            uint64_t m2, nz = 0;
            const uint16_t *kgv;
            const uint64_t *pw[16] __attribute__((aligned(64)));
            if (gray) {
                int v, up, wu, wd_;
                psk_u128 du, dd;
                jj = __builtin_ctz((unsigned)gray);
                cls ^= 1 << jj;
                up = (((cls >> jj) & 1) != 0) == (sgn > 0);
                du = dbase + tm2[jj];
                wu = du >= ua;
                du -= wu ? ua : 0;
                wd_ = dbase >= tm2[jj];
                dd = dbase - tm2[jj] + (wd_ ? 0 : ua);
                dbase = up ? du : dd;
                v = up ? wu : 3 - wd_;
                {
                    const uint16_t *dv = V[jj][v];
#if PSK_WIDE
                    ix0 = _mm256_add_epi16(ix0, _mm256_load_si256((const void *)dv));
                    ix1 = _mm256_add_epi16(ix1, _mm256_load_si256((const void *)(dv + 16)));
                    ix0 = _mm256_min_epu16(ix0, _mm256_sub_epi16(ix0, pv0));
                    ix1 = _mm256_min_epu16(ix1, _mm256_sub_epi16(ix1, pv1));
#else
                    for (i = 0; i < PSK_NP; i++) {
                        uint16_t x = (uint16_t)(idx[i] + dv[i]);
                        uint16_t y = (uint16_t)(x - c->p[i]);
                        idx[i] = x < y ? x : y;
                    }
#endif
                }
            }
            d = dbase;
            if (sgn > 0) {
                g = ((res - cr4 - 1 - (int)((uint64_t)d & 3)) * ua4) & (stp - 1);
                d += (psk_u128)(uint64_t)g * ua;
                m2 = c->mask2[(K & 15) * 64 + ((cr64 + 1 + (int)((uint64_t)d & 63)) & 63)];
            } else {
                g = ((cr4 - (int)((uint64_t)d & 3) - res) * ua4) & (stp - 1);
                d += (psk_u128)(uint64_t)g * ua;
                m2 = c->mask2[(K & 15) * 64 + ((cr64 + 128 - (int)((uint64_t)d & 63)) & 63)];
            }
            kgv = c->kg[si][g];
#if PSK_WIDE
            {
                __m256i k0a = _mm256_add_epi16(ix0, _mm256_loadu_si256((const void *)kgv));
                __m256i k0b = _mm256_add_epi16(ix1, _mm256_loadu_si256((const void *)(kgv + 16)));
                __m256i kc;
                k0a = _mm256_min_epu16(k0a, _mm256_sub_epi16(k0a, pv0));
                k0b = _mm256_min_epu16(k0b, _mm256_sub_epi16(k0b, pv1));
                _mm256_store_si256((void *)k0v, k0a);
                _mm256_store_si256((void *)(k0v + 16), k0b);
                kc = _mm256_permutex2var_epi16(k0a, permv, k0b);
                _mm512_store_si512((void *)pw, _mm512_add_epi64(rp0, _mm512_slli_epi64(_mm512_cvtepu16_epi64(_mm256_castsi256_si128(kc)), 3)));
                _mm512_store_si512((void *)(pw + 8), _mm512_add_epi64(rp1, _mm512_slli_epi64(_mm512_cvtepu16_epi64(_mm256_extracti128_si256(kc, 1)), 3)));
            }
#else
            for (i = 0; i < PSK_NP; i++) {
                uint16_t x = (uint16_t)(idx[i] + kgv[i]);
                uint16_t y = (uint16_t)(x - c->p[i]);
                k0v[i] = x < y ? x : y;
            }
            for (i = 0; i < 16; i++) pw[i] = rowp[i] + k0v[l1[i]];
#endif
            {
                const PSK_BLK mv = PSK_SET1(m2);
                for (wd = 0; wd < nw; wd += PSK_BW) {
                    PSK_BLK x0 = PSK_AND(mv, PSK_LOAD(pw[0] + wd)), x1 = PSK_AND(PSK_LOAD(pw[1] + wd), PSK_LOAD(pw[2] + wd));
                    PSK_BLK x2 = PSK_AND(PSK_LOAD(pw[3] + wd), PSK_LOAD(pw[4] + wd)), x3 = PSK_AND(PSK_LOAD(pw[5] + wd), PSK_LOAD(pw[6] + wd));
                    PSK_BLK x4 = PSK_AND(PSK_LOAD(pw[7] + wd), PSK_LOAD(pw[8] + wd)), x5 = PSK_AND(PSK_LOAD(pw[9] + wd), PSK_LOAD(pw[10] + wd));
                    PSK_BLK x6 = PSK_AND(PSK_LOAD(pw[11] + wd), PSK_LOAD(pw[12] + wd)), x7 = PSK_AND(PSK_LOAD(pw[13] + wd), PSK_LOAD(pw[14] + wd));
                    x0 = PSK_AND(x0, PSK_LOAD(pw[15] + wd));
                    x0 = PSK_AND(PSK_AND(PSK_AND(x0, x1), PSK_AND(x2, x3)), PSK_AND(PSK_AND(x4, x5), PSK_AND(x6, x7)));
                    PSK_STORE(iv + wd, x0);
                    nz |= (uint64_t)PSK_NZ(x0) << wd;
                }
                if (nw & (PSK_BW - 1)) nz &= (((uint64_t)1 << nw) - 1);
            }
            while (nz) {                                  /* rare: a word that passed 16 primes */
                uint64_t x;
                int bb;
                wd = __builtin_ctzll(nz);
                nz &= nz - 1;
                x = iv[wd];
                for (bb = 16; bb < n1 && x; bb++) {
                    const int q = l1[bb];
                    x &= c->pat[q][rowt[q]][k0v[q] + wd];
                }
                if (x && nout < maxout) {
                    PskHit *o = out + nout++;
                    o->kidx = (int32_t)m; o->wd = wd; o->bits = x;
                    o->dlo = (uint64_t)d; o->dhi = (uint64_t)(d >> 64);
                    o->crlo = (uint64_t)cr; o->crhi = (uint64_t)(cr >> 64);
                    o->stp = stp; o->pad = 0;
                }
            }
        }
        }
#if PSK_PROF
        c->tsc[3] += __rdtsc() - t_d;
#endif
    }
    *done = m;
    return nout;
}

/* ---- which primes below 65536 divide v (v < 2^126)? All primes are tested at once in vector loops:
 * v = sum l_i * 4096^i, so v mod p = sum l_i * (4096^i mod p) (below 2^32), and p | S  <=>  S * p^-1 mod 2^32 <= (2^32-1)/p. */
#define PSK_TDMAX 6560
static uint32_t psk_tdp[PSK_TDMAX], psk_tdinv[PSK_TDMAX], psk_tdlim[PSK_TDMAX];
static uint32_t psk_tdc[11][PSK_TDMAX];
static int psk_tdn = 0;

static void psk_td_init(void) {
    int i, j, n = 0;
    for (i = 3; i < 65536; i += 2) {
        int isp = 1;
        for (j = 3; j * j <= i; j += 2) if (i % j == 0) { isp = 0; break; }
        if (!isp) continue;
        {
            uint32_t p = (uint32_t)i, x = p;
            uint64_t cc = 1;
            for (j = 0; j < 5; j++) x *= 2 - p * x;
            psk_tdp[n] = p; psk_tdinv[n] = x; psk_tdlim[n] = 0xFFFFFFFFu / p;
            for (j = 0; j < 11; j++) { psk_tdc[j][n] = (uint32_t)cc; cc = (cc * 4096) % p; }
            n++;
        }
    }
    psk_tdn = n;
}

/* indices (into psk_tdp) of the primes dividing v, ascending; returns how many (at most maxout are written) */
static int psk_divscan(psk_u128 v, uint16_t *out, int maxout) {
    uint32_t l[11];
    uint32_t acc[512] __attribute__((aligned(64)));
    int nl = 0, k0, k, i, nout = 0;
    if (psk_tdn == 0) psk_td_init();
    while (v) { l[nl++] = (uint32_t)((uint64_t)v & 4095); v >>= 12; }
    if (nl == 0) return 0;
    for (k0 = 0; k0 < psk_tdn; k0 += 512) {
        const int kn = psk_tdn - k0 < 512 ? psk_tdn - k0 : 512;
        uint32_t any = 0;
        for (k = 0; k < kn; k++) acc[k] = l[0];
        for (i = 1; i < nl; i++) {
            const uint32_t li = l[i];
            const uint32_t *restrict ci = psk_tdc[i] + k0;
            for (k = 0; k < kn; k++) acc[k] += li * ci[k];
        }
        for (k = 0; k < kn; k++) {
            acc[k] = (uint32_t)(acc[k] * psk_tdinv[k0 + k]) <= psk_tdlim[k0 + k];
            any |= acc[k];
        }
        if (any)
            for (k = 0; k < kn; k++)
                if (acc[k] && nout < maxout) out[nout++] = (uint16_t)(k0 + k);
    }
    return nout;
}

/* ---- a nontrivial factor of an odd composite c < 2^62 (Pollard rho, Brent's variant, Montgomery arithmetic).
 * Returns 0 if none is found within the iteration limit. */
static inline uint64_t psk_mm(uint64_t a, uint64_t b, uint64_t c, uint64_t ninv) {
    const psk_u128 t = (psk_u128)a * b;
    const uint64_t m = (uint64_t)t * ninv;
    const uint64_t u = (uint64_t)((t + (psk_u128)m * c) >> 64);
    return u >= c ? u - c : u;
}

static uint64_t psk_gcd64(uint64_t a, uint64_t b) {
    while (b) { const uint64_t t = a % b; a = b; b = t; }
    return a;
}

static uint64_t psk_rho(uint64_t c, long maxit) {
    uint64_t ninv = c, one, add, x, y, q, g = 1, ys = 0;
    long it = 0, r = 1, k, i;
    int tries;
    for (i = 0; i < 6; i++) ninv *= 2 - c * ninv;          /* c^-1 mod 2^64 */
    ninv = (uint64_t)0 - ninv;
    one = (uint64_t)((((psk_u128)1) << 64) % c);
    for (tries = 1; tries <= 3; tries++) {
        add = one;
        for (i = 1; i < tries; i++) { add += one; if (add >= c) add -= c; }
        y = one + one; if (y >= c) y -= c;
        q = one; g = 1; r = 1; x = y;
        while (g == 1 && it < maxit) {
            x = y;
            for (i = 0; i < r; i++) { y = psk_mm(y, y, c, ninv) + add; if (y >= c) y -= c; }
            k = 0;
            while (k < r && g == 1) {
                const long lim = r - k < 64 ? r - k : 64;
                ys = y;
                for (i = 0; i < lim; i++) {
                    y = psk_mm(y, y, c, ninv) + add; if (y >= c) y -= c;
                    q = psk_mm(q, x > y ? x - y : y - x, c, ninv);
                }
                g = psk_gcd64(q, c);
                k += lim; it += lim;
            }
            r *= 2;
        }
        if (g == c) {                                       /* the batch overshot: redo it one step at a time */
            g = 1;
            for (i = 0; i < 80 && g == 1; i++) {
                ys = psk_mm(ys, ys, c, ninv) + add; if (ys >= c) ys -= c;
                g = psk_gcd64(x > ys ? x - ys : ys - x, c);
            }
        }
        if (g != 1 && g != c && g != 0) return g;
        if (it >= maxit) break;
    }
    return 0;
}
