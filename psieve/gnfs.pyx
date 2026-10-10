# cython: language_level=3
# Degree-d number field sieve for the SIQS project: it produces b-smooth relations for the SIQS matrix.
#
# The polynomial search (Kleinjung 2006), the polynomial arithmetic, the root finding and the lifting square root
# are taken from the project "General-number-field-sieve-Python" (its README lists the papers they implement), with
# the log-file calls replaced by prints. Changed from upstream: irreducibility() (degree 2), root_sieve() and root_sieve2() (same result,
# faster), and the square root, which is rewritten below (_Fq, _lift_sqrt). The upstream sieve, large primes, matrix and Lanczos/Wiedemann are not used:
# the sieve, the relation columns and the elimination at the bottom of this file are specific to the hybrid.

import math
import random
from timeit import default_timer
from datetime import datetime
import numpy as np
cimport cython
from libc.math cimport frexp, fabs
from numpy.polynomial.legendre import leggauss

def _log(line):
    print("[i]NFS poly: "+line.strip())


##### from upstream utils.py ############################################################

def invmod(a, m):
    (r, u, R, U) = (a, 1, m, 0)
    while R:
        q = r//R
        (r, u, R, U) = (R, U, r - q *R, u - q*U)
    return u%m

def compute_legendre_character(a, n):
    a = a%n
    t = 1
    while a:
        while not a&1:
            a = a>>1
            if n%8 == 3 or n%8 == 5: t = -t
        a, n = n, a
        if a%4 == n%4 and n%4 == 3: t = -t
        a = a%n
    if n == 1: return t
    return 0

def compute_sqrt_mod_p(n, p):
    n %= p
    if n == 1 : return 1
    P = p-1
    z = int(random.randint(2, P))
    while compute_legendre_character(z, p) != -1:
        z = int(random.randint(2, P))
    r = 0
    while not P&1:
        P >>= 1
        r += 1
    s = P
    generator = pow(z, s, p)
    lbd = pow(n, s, p)
    omega = pow(n, (s+1)>>1, p)

    while True:
        if not lbd: return 0
        if lbd == 1: return omega
        for m in range(1, r):
            if pow(lbd, 1<<m, p)==1: break

        tmp = pow(2, r-m-1, p-1)
        lbd = lbd*pow(generator, tmp<<1, p)%p
        omega = omega*pow(generator, tmp, p)%p

def is_prime(n):
    if n == 2 or n == 3:
        return True
    if not n&1 or n==1:
        return False
    r, s = 0, n - 1
    while not s&1:
        r += 1
        s >>= 1
    for _ in range(0, 50):
        a = random.randrange(2, n - 1)
        x = pow(a, s, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(0, r - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                break
        else:
            return False
    return True

def dickman(x, table):
    k,res = math.ceil(x),0
    delta = k-x
    if k-1 > len(table): table = get_dickman_table(k)
    tmp = 1
    for i in range(len(table[k-2])):
        res += table[k-2][i]*tmp
        tmp *= delta

    return res, table

def get_dickman_table(k):
    coeffs = [[1-math.log(2)]+[1/(i*(1<<i)) for i in range(1, 30)]]
    for i in range(3, k+1):
        new = [0]*30
        for u in range(1, 30):
            c = 0
            for j in range(u): c += coeffs[-1][j]/(u*pow(i, u-j))
            new[u] = c

        c = 0
        for j in range(1, len(new)): c += new[j-1]/(j+1)
        new[0] = c/(i-1)
        coeffs.append(new)

    return coeffs

def fac(n):
    res = 1
    for i in range(2, n+1): res *= i
    return res

def binom(k, n):
    res = 1

    for i in range(n-k+1, n+1): res *= i
    for i in range(2, k+1): res //= i

    return res

def central(k,x): return pow(x, k/2-1)/(math.exp(x/2)*pow(2, k/2)*math.gamma(k/2))

def non_central(k, l ,x):
    if x <= 0: return 0
    res = 0
    for i in range(100): res += central(k+2*i, x)*pow(l/2, i)/fac(i)
    
    return res*math.exp(-l/2)

def format_duration(delta):
    hours = delta.seconds//3600
    minutes = delta.seconds//60 - hours*60
    seconds = delta.seconds%60
    return str(delta.days)+" days, "+str(hours)+" hours, "+str(minutes)+" minutes, "+str(seconds)+" seconds"


##### from upstream polynomial_functions.py ############################################################

def get_derivative(f):
    res = [0]*(len(f)-1)
    for i in range(len(f)-1):
        res[i] = (len(f)-1-i)*f[i]
    return res

def power(poly, f, p, exp):
    if exp == 1: return poly

    tmp = power(poly, f, p, exp>>1)
    tmp = poly_prod(tmp, tmp)

    if exp&1:
        tmp = poly_prod(tmp, poly)
        return div_poly_mod(tmp, f, p)
    
    else: return div_poly_mod(tmp, f, p)

def shift(poly,k):
    res = [i for i in poly]
    for i in range(len(poly)-1):
        for j in range(i+1, len(poly)):
            res[j] += binom(j-i, len(poly)-1-i)*pow(k, j-i)*poly[i]
    return res

def irreducibility(f, p):
    g = [1,0]
    for i in range((len(f)-1)//2):
        g = power(g, f, p, p)
        tmp2 = [u for u in g]
        if len(tmp2) == 1: tmp2 = [-1, tmp2[0]]
        else: tmp2[-2] -= 1
        tmp = gcd_mod(f, tmp2, p)
        if len(tmp) > 1: return False
    return True

def find_roots_poly(f, p):
    tmp_f = [i%p for i in f]
    for k in range(len(f)):
        if tmp_f[k]:
            tmp_f = tmp_f[k:]
            break

    r = []
    tmp = [1,0]
    g = [1]
    tmp_p = p
    while tmp_p>1:
        if tmp_p&1:
            g = div_poly_mod(poly_prod(g, tmp), tmp_f, p)
        tmp = div_poly_mod(poly_prod(tmp, tmp), tmp_f, p)
        tmp_p >>= 1

    g = div_poly_mod(poly_prod(g, tmp), tmp_f, p)
    if len(g) == 1: g = [-1, g[0]]
    else: g[-2] -= 1
    g = gcd_mod(f, g, p)
    if g[-1] == 0:
        r.append(0)
        del g[-1]
    return r + roots(g, p)

def roots(g, p):
    if len(g) == 1: return []
    if len(g) == 2: return [-g[1]*invmod(g[0], p)%p]
    if len(g) == 3:
        tmp = (g[1]*g[1]-4*g[0]*g[2])%p
        if tmp == 0: return [-g[1]*invmod(2*g[0], p)%p]
        if compute_legendre_character(tmp, p) == -1: return []
        tmp = compute_sqrt_mod_p(tmp, p)*invmod(2*g[0], p)%p
        return [(-g[1]*invmod(g[0]<<1, p)+tmp)%p, (-g[1]*invmod(g[0]<<1, p)-tmp)%p]
    
    h = [1]
    while len(h) == 1 or h == g:
        a = random.randint(0, p-1)
        h = power([1, a], g, p, (p-1)>>1)
        for k in range(len(h)):
            if h[k]:
                h = h[k:]
                break
        h[-1] -= 1
        h = gcd_mod(h, g, p)
    r = roots(h, p)
    h = quotient_poly_mod(g, h, p)
    return r+roots(h, p)

def get_complex_roots(f):
    rho = 0
    for i in range(1, len(f)): rho += abs(f[i])
    rho /= abs(f[0])
    rho = max(1,rho)
    d = len(f)-1
    roots = []
    for i in range(d): roots.append(rho*pow(math.cos(2*math.pi/d)+math.sin(2*math.pi/d)*1j, i))
    next_roots = [None]*d

    for _ in range(1000):
        for i in range(d):
            prod = f[0]
            for k in range(d):
                if k != i: prod *= (roots[i]-roots[k])
            next_roots[i] = roots[i]-evaluate(f, roots[i])/prod
        roots = [i for i in next_roots]
    return roots

def poly_prod(a, b):
    res = [0]*(max(len(a), len(b))+min(len(a), len(b))-1)

    for i in range(len(a)):
        for j in range(len(b)):
            res[i+j] += a[i]*b[j]

    return res

def div_poly(a,b):
    remainder = [i for i in a]
    difference = len(a)-len(b)+1
    leading_coeff = b[0]

    for j in range(difference):
        quotient = -remainder[j]//leading_coeff
        for k in range(len(b)):
            remainder[j+k] += quotient*b[k]
            
    for k in range(len(remainder)):
        if remainder[k]: return remainder[k:]
        
    return [0]

def div_poly_mod(a, tmp_b, p):
    remainder = [i%p for i in a]
    b = [i%p for i in tmp_b]
    
    #print(remainder, b)
    while not b[0]: del b[0]
    
    difference = len(a)-len(b)+1
    coeff = invmod(-b[0], p)
    for j in range(difference):
        if remainder[j]:
            quotient = remainder[j]*coeff%p
            remainder[j] = 0
            for k in range(1,len(b)): remainder[j+k] = (remainder[j+k]+quotient*b[k]%p)%p
            
    for k in range(len(remainder)):
        if remainder[k]: return remainder[k:]
        
    return [0]

def quotient_poly_mod(a, b, p):
    remainder = [i%p for i in a]
    b = [i%p for i in b]
    
    while not b[0]: del b[0]
    
    difference = len(a)-len(b)+1
    coeff = invmod(-b[0], p)
    res = [0]*difference

    for j in range(difference):
        quotient = remainder[j]*coeff%p
        res[j] = -quotient
        for k in range(len(b)):
            remainder[j+k] = (remainder[j+k]+quotient*b[k]%p)%p
            
    for k in range(len(res)):
        if res[k]: return res[k:]
        
    return [0]

def gcd_mod(f, poly, p):
    while poly != [0]*len(poly):
        (f, poly) = (poly, div_poly_mod(f, poly, p))
    return f

def eval_mod(f, x, n):
    res = 0

    for i in range(len(f)-1):
        res += f[i]
        res *= x
        res %= n

    res += f[-1]
    res %= n

    return res

def evaluate(f, x):
    res = 0

    for i in range(len(f)-1):
        res += f[i]
        res *= x

    res += f[-1]

    return res

def eval_F(x, y, f, d):
    tmp = 0
    tmp2 = 1

    for k in range(d):
        tmp += f[k]*tmp2
        tmp *= x
        tmp2 *= y

    tmp += f[d]*tmp2

    return tmp


##### from upstream utils_polynomial_selection.py ############################################################

def minimise_Lnom(f, s, B, m0, m1):
    k, res = 1, get_Lnorm(poly_prod(f,f),s,B)
    while k > 0:
        flag = False
        F = shift(f, -k)
        tmp = get_Lnorm(poly_prod(F, F), s, B)
        if tmp < res:
            res, f, k, flag, m0 = tmp, [i for i in F], k<<1, True, m0+k*m1

        F = shift(f, k)
        tmp = get_Lnorm(poly_prod(F, F), s, B)
        if tmp < res:
            res, f, k, flag, m0 = tmp, [i for i in F], k<<1, True, m0-k*m1
        if not flag: k >>= 1

    return f, m0

def evaluate_polynomial_quality(f_x, B, m0, m1, primes, LOG_PATH):
    _, s = get_sieve_region(f_x, B)

    if len(f_x) > 2:

        if (s*abs(f_x[-3])-abs(f_x[-2]))//m0 > 0 and (s*s*abs(f_x[-3])-abs(f_x[-1]))//m0 > 0:
            f_x = root_sieve2(f_x, [m1,-m0], primes[:150], (s*abs(f_x[-3])-abs(f_x[-2]))//m0, (s*s*abs(f_x[-3])-abs(f_x[-1]))//m0)
            _, s = get_sieve_region(f_x, B)
            f_x, m0 = minimise_Lnom(f_x, s, B, m0, m1)

        elif (s*abs(f_x[-2])-abs(f_x[-1]))//m0 > 0:
            f_x = root_sieve(f_x, [m1,-m0], primes[:150], (s*abs(f_x[-2])-abs(f_x[-1]))//m0)
            _, s = get_sieve_region(f_x, B)
            f_x, m0 = minimise_Lnom(f_x, s, B, m0, m1)

    M, s = get_sieve_region(f_x, B)
    table = get_dickman_table(20)

    alpha = alpha_score(f_x, primes[:350])
    E_score = get_Escore(f_x, [m1, -m0], alpha[0], B, M, round(B/math.sqrt(s)), table)
    E_score_2 = get_Epscore(f_x, [m1, -m0], alpha, B, M, round(B/math.sqrt(s)), table)
    L_norm = get_Lnorm(poly_prod(f_x, f_x), s, B)

    _log("f(x) = "+str(f_x)+"    g(x) = "+str([m1, -m0]))
    _log("alpha = "+str(alpha[0]))
    _log("E-score = "+str(E_score))
    _log("E-score 2 = "+str(E_score_2))
    _log("L²-norm = "+str(L_norm))
    _log("skew = "+str(s)+"\n")
    
    return f_x, m0, M

def get_sieve_region(f, B):
    F = poly_prod(f, f)
    d = len(f)-1
    ratios = [math.log(1e-7+abs(f[i+1]/(1+abs(f[i])))) for i in range(d)]
    s = 2*int(math.exp(sum(ratios)/d)) # skew factor
    k = 1 # shift
    best_norm = None

    while k > 0:
        updated = False

        if s-k > 0:
            norm = get_Lnorm(F, s-k, B)
            if best_norm == None or norm < best_norm:
                s = s-k
                k <<= 1
                best_norm = norm
                updated = True

        if s+k < B:
            norm = get_Lnorm(F, s+k, B)
            if best_norm == None or norm < best_norm:
                s = s+k
                k <<= 1
                best_norm = norm
                updated = True

        if not updated: k >>= 1

    x = round(B*math.sqrt(s))

    return x, s

def get_Lnorm(F, s, B):
    sqrt = math.sqrt(s)
    n = len(F)
    
    base_X, base_Y = B*sqrt/2, 1+B/sqrt
    current_X, current_Y = pow(base_X, n), base_Y
    base_X, base_Y = base_X*base_X, base_Y*base_Y

    res = 0.0

    for i in range(0, n, 2):
        res += (2*current_X)*F[i]*(current_Y-1)/((i+1)*(n-i))
        current_X /= base_X
        current_Y *= base_Y

    return math.log(abs(res))/2

def get_Escore(f, g, alpha, B, x_limit, y_limit, table):
    n = 32
    log_B = math.log(B)
    
    x, w = leggauss(n)  # nodes & weights on [-1, 1]
    a, b = 0, math.pi
    xm = 0.5 * (b - a) * x + 0.5 * (b + a)
    wm = 0.5 * (b - a) * w

    res = 0
    for i in range(n):
        angle = xm[i]
        X = x_limit*math.cos(angle)
        Y = (y_limit-1)*math.sin(angle)+1
        res1, table = dickman((math.log(abs(eval_F(X, Y, f, len(f)-1)))-alpha)/log_B, table)
        res2, table = dickman((math.log(abs(eval_F(X, Y, g, len(g)-1))))/log_B, table)
        res += wm[i]*res1*res2

    return res

def get_Epscore(f, g, alpha, B, x_limit, y_limit, table):
    k,l,c = alpha[1],alpha[2],alpha[3]

    n = 16
    log_B = math.log(B)

    # Gauss–Legendre for mu in [c, c + 5c]
    mu_nodes, mu_weights = leggauss(n)
    a_mu, b_mu = c, 6*c
    mu_x = 0.5*(b_mu-a_mu)*mu_nodes + 0.5*(b_mu+a_mu)
    mu_w = 0.5*(b_mu-a_mu)*mu_weights

    # Gauss–Legendre for theta in [0, π]
    th_nodes, th_weights = leggauss(n)
    a_th, b_th = 0, math.pi
    th_x = 0.5*(b_th-a_th)*th_nodes + 0.5*(b_th+a_th)
    th_w = 0.5*(b_th-a_th)*th_weights

    res = 0
    for i in range(n):
        Y = non_central(k, l, mu_x[i])
        for j in range(n):
            angle = th_x[j]
            X_eval = x_limit*math.cos(angle)
            Y_eval = (y_limit-1)*math.sin(angle)+1
            res1, table = dickman((math.log(abs(eval_F(X_eval, Y_eval, f, len(f)-1)))-(mu_x[i]-c))/log_B, table)
            res2, table = dickman((math.log(abs(eval_F(X_eval, Y_eval, g, len(g)-1))))/log_B, table)
            res += mu_w[i]*th_w[j]*res1*res2*Y

    return res

def alpha_score(f, primes):
    E, F = 0, 0
    evals = [0]*primes[-1]
    tmp = [0]*len(f)
    for j in range(len(f)): tmp[j] = evaluate(f, j)
    for q in range(1, len(f)):
        for k in range(len(f)-1, q-1, -1): tmp[k] -= tmp[k-1]

    eval = tmp[0]
    for k in range(primes[-1]):
        evals[k] = eval
        eval += tmp[1]
        for q in range(1, len(f)-1): tmp[q] += tmp[q+1]

    f_prime = get_derivative(f)
    upto = len(f)+10
    baseline_term = 0
    for p in primes:
        log_p = math.log(p)
        baseline_term += log_p/(p-1)
        Ep,Fp = 0,0
        ramified_roots = []
        ramified = False
        for r in range(p):
            if not evals[r]%p:
                if not eval_mod(f_prime, r, p):
                    ramified = True
                    ramified_roots.append(r)
                    Ep += 1/p
                    Fp += 1/p
                else:
                    Ep += 1/(p-1)
                    Fp += (p+1)/((p-1)**2)

        if ramified:
            tmp2 = p
            for i in range(2, upto):
                new = []
                for r in ramified_roots:
                    if not eval_mod(f, r, tmp2*p):
                        for k in range(p): new.append(r+k*tmp2)
                        Ep += 1/tmp2
                        Fp += (2*i-1)/tmp2
                tmp2 *= p
                ramified_roots = new.copy()

            Ep += len(ramified_roots)/(pow(p, upto-2)*(p-1))
            Fp += len(ramified_roots)*(upto*upto+2*upto/(p-1)+(p+1)/((p-1)**2))/tmp2

        if not f[0]%p:
            frev = f[::-1]
            fdrev = get_derivative(frev)
            if eval_mod(fdrev,0,p):
                Ep += 1/(p-1)
                Fp += (p+1)/((p-1)**2)
            else:
                Ep += 1/p
                Fp += 1/p
                ramified_roots = [0]
                tmp2 = p
                for i in range(2, upto):
                    new = []
                    for r in ramified_roots:
                        if not eval_mod(frev, r, tmp2*p):
                            for k in range(p): new.append(r+k*tmp2)
                            Ep += 1/tmp2
                            Fp += (2*i-1)/tmp2
                    tmp2 *= p
                    ramified_roots = new.copy()

                Ep += len(ramified_roots)/(pow(p, upto-2)*(p-1))
                Fp += len(ramified_roots)*(upto*upto+2*upto/(p-1)+(p+1)/((p-1)**2))/tmp2
                
        tmpE, tmpF = p*Ep/(p+1), p*Fp/(p+1)
        E += tmpE*log_p
        F += (tmpF-tmpE*tmpE)*log_p*log_p
        
    k = 2*E-F/2
    lbd = F/2-E
    return E-baseline_term, k, lbd, baseline_term

def local_opt(f, g, B):
    m0, m1 = -g[1], g[0]

    d = len(f)-1
    if d <= 5:
        poly = []
        for i in range(3): poly.append(f[i]*binom(2-i, d-i))
        poly = poly_prod(poly, poly)
        poly2 = get_derivative(poly)
        zeros = get_complex_roots(poly2)
        min,mink = f[2]**2,0
        for r in zeros:
            if r.imag == 0:
                if evaluate(poly, math.ceil(r.real)) < min:
                    min, mink = evaluate(poly, math.ceil(r.real)), math.ceil(r.real)
                if evaluate(poly, math.floor(r.real)) < min:
                    min, mink = evaluate(poly, math.floor(r.real)), math.floor(r.real)
    else:
        poly = []
        for i in range(4): poly.append(f[i]*binom(3-i, d-i))
        zero = get_complex_roots(poly)
        min, mink = abs(f[3]), 0
        for r in zero:
            if r.imag == 0:
                if abs(evaluate(poly, math.ceil(r.real))) < min:
                    min, mink = abs(evaluate(poly, math.ceil(r.real))), math.ceil(r.real)
                if abs(evaluate(poly, math.floor(r.real))) < min:
                    min, mink = abs(evaluate(poly, math.floor(r.real))), math.floor(r.real)

    f = shift(f, mink)
    m0 -= mink*m1

    k, u, v, iteration = 1, 1, 1, 0
    _, s = get_sieve_region(f, B)
    min_n = get_Lnorm(poly_prod(f, f), s, B)
    while iteration < 500 and (k > 0 or u > 0 or v > 0):
        F = shift(f, -k)
        tmp = get_Lnorm(poly_prod(F, F), s, B)
        flag = False
        if tmp < min_n: f, min_n, k, flag, m0 = F.copy(), tmp, k<<1, True, m0+k*m1

        F = shift(f, k)
        tmp = get_Lnorm(poly_prod(F, F), s, B)
        if tmp < min_n: f, min_n, k, flag, m0 = F.copy(), tmp, k<<1, True, m0-k*m1
        if flag and not u: u = 1
        if flag and not v: v = 1
        elif not flag: k >>= 1

        if u:
            flag = False
            F = f.copy()
            F[-3] += u*m1
            F[-2] -= u*m0
            tmp = get_Lnorm(poly_prod(F, F), s, B)
            if tmp < min_n: f, min_n, u, flag = F.copy(), tmp, u<<1, True

            F = f.copy()
            F[-3] -= u*m1
            F[-2] += u*m0
            tmp = get_Lnorm(poly_prod(F, F), s, B)
            if tmp < min_n: f, min_n, u, flag = F.copy(), tmp, u<<1, True
            if flag and not k: k = 1
            if flag and not v: v = 1
            elif not flag: u >>= 1

        if v:
            flag = False
            F = f.copy()
            F[-2] += v*m1
            F[-1] -= v*m0
            tmp = get_Lnorm(poly_prod(F, F), s, B)
            if tmp < min_n: f, min_n, v, flag = F.copy(), tmp, v<<1, True

            F = f.copy()
            F[-2] -= v*m1
            F[-1] += v*m0
            tmp = get_Lnorm(poly_prod(F, F), s, B)
            if tmp < min_n: f, min_n, v, flag = F.copy(), tmp, v<<1, True
            if flag and not k: k = 1
            if flag and not u: u = 1
            elif not flag: v >>= 1

        iteration += 1
        _, s = get_sieve_region(f, B)

    return f, m0

def prime_combinations_with_indices(Q, l, B):
    n = len(Q)
    indices = [0] * l  # reuse buffer

    def backtrack(start, depth, product):
        if depth == l:
            yield tuple(indices)
            return
        for i in range(start, n - (l - depth) + 1):
            p = Q[i]
            if product * p >= B:
                break  # sorted Q means all further i will be too large
            indices[depth] = i
            yield from backtrack(i + 1, depth + 1, product * p)

    yield from backtrack(0, 0, 1)

def compute_e(m0, root_used, NB_ROOTS, prod, a_d, n, d):
    e = []
    tmp_m = m_mu(m0, root_used, [0]*NB_ROOTS, NB_ROOTS)
    base = poly(tmp_m, prod, a_d, n, d)[1]%prod
    for i in range(NB_ROOTS):
        line = [0]*d
        for j in range(d):
            if not i:
                tmp_m = m_mu(m0, root_used, [j]+[0]*(NB_ROOTS-1), NB_ROOTS)
                line[j] = poly(tmp_m, prod, a_d, n, d)[1]%prod
            elif j:
                tmp_m = m_mu(m0, root_used, [0]*i+[j]+[0]*(NB_ROOTS-i-1), NB_ROOTS)
                line[j] = (poly(tmp_m, prod, a_d, n, d)[1]-base)%prod
        e.append(line)

    return e

def compute_f(n, a_d, m0, d, prod, root_used, NB_ROOTS, e):
    f0 = (n-a_d*pow(m0, d))/(prod*prod*pow(m0, d-1))
    f = []

    for i in range(NB_ROOTS):
        line = [0]*d
        for j in range(d):
            line[j] = -(a_d*d*root_used[i][j]/pow(prod, 2)+e[i][j]/prod)
        f.append(line)

    return f, f0

def create_first_array(NB_ROOTS, f0, f, d):
    vec = [0]*(NB_ROOTS>>1)
    array1 = []

    while vec[-1] < d:
        U = (f0 + sum([f[j][vec[j]] for j in range(NB_ROOTS>>1)]))%1

        if not len(array1) or U > array1[-1][0]: array1.append([U, [i for i in vec]])
        else:
            tmp_a = 0
            tmp_b = len(array1)-1
            tmp = (tmp_a+tmp_b)>>1
            while tmp_a <= tmp_b:
                if array1[tmp][0] > U: tmp_b = tmp-1
                else: tmp_a = tmp+1
                tmp = (tmp_a+tmp_b)>>1
            array1.insert(tmp_a, [U, [i for i in vec]])
        vec[0] += 1
        for j in range(len(vec)-1):
            if vec[j] == d:
                vec[j] = 0
                vec[j+1] += 1
            else: break

    return array1

def create_second_array(NB_ROOTS, len_vec, d, f):
    vect = [0]*(NB_ROOTS-len_vec)
    array2 = []
    while vect[-1] < d:
        U = -sum([f[len_vec+j][vect[j]] for j in range(len(vect))])%1

        if not len(array2) or U > array2[-1][0]: array2.append([U, [i for i in vect]])
        else:
            tmp_a = 0
            tmp_b = len(array2)-1
            tmp = (tmp_a+tmp_b)>>1
            while tmp_a <= tmp_b:
                if array2[tmp][0] > U: tmp_b = tmp-1
                else: tmp_a = tmp+1
                tmp = (tmp_a+tmp_b)>>1
            array2.insert(tmp_a, [U, [i for i in vect]])
        vect[0] += 1
        for j in range(len(vect)-1):
            if vect[j] == d:
                vect[j] = 0
                vect[j+1] += 1
            else: break

    return array2

def m_mu(m0,roots,vec,l):
    return m0 + sum([roots[i][vec[i]] for i in range(l)])

def poly(m0, m1, a_d, n, d):
    c = [a_d]
    r = [n]
    for i in range(d-1, -1, -1):
        r.append((r[-1]-c[-1]*pow(m0, i+1))//m1)
        delta = -r[-1]*m1*invmod(m1,pow(m0, i))%(m1*pow(m0, i))
        c.append((r[-1]+delta)//pow(m0, i))
    for i in range(1, len(c)):
        if c[i] > m0//2:
            c[i] -= m0
            c[i-1] += m1
        elif c[i] < -m0//2:
            c[i] += m0
            c[i-1] -= m1
    return c


##### from upstream mono_cpu_polynomial_selection.py ############################################################

def Kleinjung_poly_search(n, primes, NB_ROOTS, PRIME_BOUND, MULTIPLIER, M, d, NB_POLY_COARSE_EVAL, NB_POLY_PRECISE_EVAL, LOG_PATH):
    t1 = datetime.now()
    P = []
    polys = []
    for p in primes:
        if p > PRIME_BOUND: break
        if p%d == 1: P.append(p)
    a_d = MULTIPLIER
    if d >= 4: admax = round(pow(pow(M, 2*d-2)/n, 1/(d-3)))
    else: admax = M

    cpt = 0
    avg = 0
    while a_d < admax and cpt < NB_POLY_COARSE_EVAL:
        tmp = a_d
        for p in primes:
            while not tmp%p: tmp//= p

        if tmp > 1: # If a_d is not primes[-1] smooth
            a_d += MULTIPLIER
            continue

        mw = math.ceil(pow(n/a_d, 1/d))
        ad1max = round(M*M/mw)
        if d > 2: ad2max = pow(pow(M, 2*d-6)/pow(mw, d-4), 1/(d-2))
        else: ad2max = M

        Q = []
        roots = []
        for p in P:
            if not a_d%p: continue

            f = [a_d]+[0]*d
            f[-1] = (-n)%p # Construct polynomial a_d*x^d - n (mod r)
            root = find_roots_poly(f, p)
            if len(root) > 0:
                Q.append(p)
                roots.append(root)

        if len(roots) >= NB_ROOTS:

            combinations = prime_combinations_with_indices(Q, NB_ROOTS, ad1max)

            for set in combinations:
                Q_used = []
                prod = 1
                for i in range(NB_ROOTS):
                    Q_used.append(Q[set[i]])
                    prod *= Q[set[i]]

                root_used = [roots[set[i]] for i in range(NB_ROOTS)]
                for i in range(NB_ROOTS): # Do some CRT
                    x = prod//Q_used[i]
                    tmp2 = x*invmod(x, Q_used[i])
                    for j in range(d): root_used[i][j] = root_used[i][j]*tmp2%prod

                m0 = mw+(-mw)%prod
                e = compute_e(m0, root_used, NB_ROOTS, prod, a_d, n, d)
                f, f0 = compute_f(n, a_d, m0, d, prod, root_used, NB_ROOTS, e)

                epsilon = ad2max/m0
                array1 = create_first_array(NB_ROOTS, f0, f, d)
                len_vec = NB_ROOTS>>1
                array2 = create_second_array(NB_ROOTS, len_vec, d, f)
                
                min = 0
                for j in range(len(array2)):
                    while min < len(array1) and array2[j][0]-epsilon > array1[min][0]: min += 1
                    if min == len(array1): break
                    z = min
                    while z < len(array1) and abs(array2[j][0]-array1[z][0]) < epsilon:
                        tmp = [poly(m_mu(m0, root_used, array1[z][1]+array2[j][1], NB_ROOTS), prod, a_d, n, d),
                               m_mu(m0, root_used, array1[z][1]+array2[j][1], NB_ROOTS),
                               prod]
                        cpt += 1
                        tmp[0], tmp[1] = local_opt(tmp[0], [prod,-tmp[1]], primes[-1])
                        _, s = get_sieve_region(tmp[0], primes[-1])
                        L_score = get_Lnorm(poly_prod(tmp[0], tmp[0]), s, primes[-1])
                        avg += L_score
                        if not len(polys):
                            tmp.append(L_score)
                            polys.append(tmp)
                        else:
                            if L_score < polys[-1][3]:
                                tmp.append(L_score)
                                a = 0
                                b = len(polys)-1
                                tmpu = (a+b)>>1
                                while a <= b:
                                    if polys[tmpu][3] < L_score: a = tmpu+1
                                    else: b = tmpu-1
                                    tmpu = (a+b)>>1
                                polys.insert(a, tmp)
                            elif len(polys) < NB_POLY_PRECISE_EVAL:
                                tmp.append(L_score)
                                polys.append(tmp)
                            if len(polys) > NB_POLY_PRECISE_EVAL: del polys[-1]

                        if cpt >= NB_POLY_COARSE_EVAL:

                            t2 = datetime.now()
                            _log("Polynomial search done in "+format_duration(t2-t1)+".\n")
                            _log(str(cpt)+" polynomials created, "+str(len(polys))+" kept for ranking")
                            _log("Average L2 score = "+str(avg/cpt))
                            _log("Ranking polynomials")

                            return select_best_poly_candidate(polys, primes)

                        z += 1
        a_d += MULTIPLIER

def select_best_poly_candidate(polys, primes):
    best_poly = None
    best_E = None
    table = get_dickman_table(10)

    for i in range(len(polys)):
        x_limit,s = get_sieve_region(polys[i][0], primes[-1])
        alpha = alpha_score(polys[i][0], primes[:300])
        polys[i].append(alpha)
        E_score = get_Epscore(polys[i][0], [polys[i][2], -polys[i][1]], alpha, primes[-1], x_limit, round(primes[-1]/math.sqrt(s)), table)
        if best_E == None or E_score > best_E:
            best_poly = polys[i]
            best_E = E_score

    return best_poly


######################################################################################################################
# Hybrid driver: degree-d number field sieve whose relations go into one matrix with the SIQS relations.
#
#   f(x)  degree d, c_d = leading coefficient, f(m0/m1) = 0 mod n        (Kleinjung polynomial search, above)
#   algebraic side: a - b*alpha, norm F(a,b) = sum f_i a^(d-i) b^i        smooth over the algebraic factor base
#   rational side:  G(a,b) = a*m1 - b*m0                                   smooth over the rational factor base
#
# nfs_sieve() sieves and hands out every relation as a sparse matrix row: the ideals with odd exponent on the algebraic
# side, the quadratic characters, and the primes with odd exponent on the rational side. A relation may keep one large
# prime per side; a large prime is just one more column. Nothing is eliminated here.
# A set S of rows whose algebraic columns and characters cancel has  g'(w)^2 * prod(c_d*a - b*w)  a square in Z[w]
# (w = c_d*alpha, root of the monic g). dep_square() takes its square root and maps it to Z/n: X with
#       X^2 = prod G(a,b)   (mod n)
# The rational primes of S that are left odd are cancelled by SIQS relations (X^2 = smooth mod n) in the same matrix.
######################################################################################################################

NFS_NB_ROOTS=3              # Kleinjung: number of primes in m1 (l in the paper)
NFS_PRIME_BOUND=300         # Kleinjung: largest prime allowed in m1
NFS_MULTIPLIER=1            # Kleinjung: the leading coefficient is a multiple of this
NFS_POLY_COARSE=100         # polynomials generated
NFS_POLY_PRECISE=50         # polynomials kept for the precise ranking
NFS_CHARS=64                # quadratic characters
NFS_SLACK=20                # sieve: verify positions whose logs fall short of the full size by at most this many bits
NFS_SMALL=30                # sieve: primes up to this are not sieved (they are covered by the slack)
NFS_LP_BITS=0               # large primes: a relation may keep a prime up to 2^NFS_LP_BITS outside the factor base (0 = off)
NFS_ALG=1<<60               # column numbers: a rational prime p is p (the sign is 1); the ideal (p, alpha - r) is NFS_ALG + (p<<28) + r, r = p at infinity
NFS_BATCH=2048              # sieve candidates whose smoothness is tested in one go
NFS_LINE_EXP=100            # order of the lines: b is taken as b * (b/phi(b))^(NFS_LINE_EXP/100), so lines with few pairs coprime to b come later
NFS_BLOCK=1<<18             # sieve: a line is processed in blocks of this many positions

_STATE={}
_T={}                       # seconds per stage of the current call

def _tick(k,t0):
    _T[k]=_T.get(k,0.0)+default_timer()-t0

try:                        # the square roots multiply integers of 10^4..10^6 bits; GMP is 10-50x faster there
    from gmpy2 import mpz as _big, gcd as _gcd
except ImportError:
    _big=int
    _gcd=math.gcd

def _is_smooth(v,fbprod):
    v=abs(v)
    if v<2:
        return v==1
    g=math.gcd(v,fbprod%v)
    while g>1:
        v//=g
        g=math.gcd(v,g)
    return v==1

def _cofactor(v,fbprod):
    # what is left of |v| after removing every prime of the factor base (fbprod = their product, as _big)
    v=_big(abs(v))
    if v<2:
        return int(v)
    g=_gcd(v,fbprod%v)
    while g>1:
        v//=g
        g=_gcd(v,g)
    return int(v)

def _cofactors(vals,fbprod):
    # _cofactor for many values at once: the product of a batch is formed pairwise, the factor base product is reduced
    # modulo it once and then down the same tree, so each value gets fbprod mod itself for a fraction of the cost of
    # reducing the whole product by every value
    out=[]
    for s0 in range(0,len(vals),NFS_BATCH):
        lv=[_big(abs(v)) for v in vals[s0:s0+NFS_BATCH]]
        tree=[lv]
        while len(tree[-1])>1:
            t=tree[-1]
            nxt=[t[i]*t[i+1] for i in range(0,len(t)-1,2)]
            if len(t)&1:
                nxt.append(t[-1])
            tree.append(nxt)
        rem=[fbprod%tree[-1][0]]
        for lvl in range(len(tree)-2,-1,-1):
            t=tree[lvl]
            rem=[rem[i>>1]%t[i] for i in range(len(t))]
        for v,r in zip(lv,rem):
            if v>1:
                g=_gcd(v,r)
                while g>1:
                    v//=g
                    g=_gcd(v,g)
            out.append(int(v))
    return out

def _factor_over(v,cand):
    # exponents of v over the primes in cand (ascending); returns ({p: e}, cofactor). A negative v gets {-1: 1}.
    ex={}
    if v<0:
        ex[-1]=1
        v=-v
    for p in cand:
        if v%p==0:
            e=0
            while v%p==0:
                v//=p
                e+=1
            ex[p]=e
    return ex,v

NFS_ROOTSIEVE_CELLS=1<<22     # root sieve: largest number of (u,v) cells looked at

def root_sieve(f, g, primes, U):
    # Same result as the upstream root_sieve; the innermost loop is one strided numpy add.
    U = max(1, min(U, NFS_ROOTSIEVE_CELLS>>1))
    array = np.zeros(U<<1)
    for p in primes:
        k = 1
        P = p
        while P < primes[-1]:
            c = math.log(p)/(pow(p, k-1)*(p+1))
            for x in range(P):
                eval1, eval2 = eval_mod(f, x, P), eval_mod(g, x, P)
                for u in range(P):
                    if (eval1+u*eval2)%P == 0:
                        array[(u+U)%P::P] += c
            k += 1
            P *= p
    u = int(np.argmax(array))-U
    f[-1] += u*g[-1]
    f[-2] += u*g[-2]
    return f

def root_sieve2(f, g, primes, U, V):
    # Same result as the upstream root_sieve2 (same additions in the same order for every cell). For a fixed x each
    # row i of the array belongs to one u = i-U mod P, hence one v. The modular inverse is taken once per x, and the
    # cells are filled with numpy: row by row when the rows are few, otherwise all rows at once per column step.
    if (U<<1)*(V<<1) > NFS_ROOTSIEVE_CELLS:               # the upstream ranges grow without limit with the size of n
        U = max(1, min(U, 32))
        V = max(1, min(V, NFS_ROOTSIEVE_CELLS//(4*U)))
    array = np.zeros((U<<1, V<<1))
    rows_i = np.arange(U<<1)
    for p in primes:
        k = 1
        P = p
        while P < primes[-1]:
            c = math.log(p)/(pow(p, k-1)*(p+1))
            um = (rows_i-U)%P
            by_row = (U<<1)*P < (V<<1)
            for x in range(P):
                eval1, eval2 = eval_mod(f, x, P), eval_mod(g, x, P)
                if eval2%p:
                    t = eval1*invmod(eval2, P)%P
                    cols = (V-(um*x+t))%P                  # first column j = v+V mod P of each row, v = -(u*x+t)
                    if by_row:
                        for i in range(U<<1):
                            array[i, int(cols[i])::P] += c
                    else:
                        while True:
                            m = cols < (V<<1)
                            if not m.any():
                                break
                            array[rows_i[m], cols[m]] += c
                            cols = cols+P
            k += 1
            P *= p
    maxu, maxv = divmod(int(np.argmax(array)), V<<1)
    u, v = maxu-U, maxv-V
    f[-1] += v*g[-1]
    f[-2] += v*g[-2]+u*g[-1]
    f[-3] += u*g[-2]
    return f

cdef enum:
    DMAX=10                 # largest supported degree

cdef void _fq_mul(const long long* a,const long long* b,const long long* gl,int d,long long p,long long* out) noexcept nogil:
    # out = a*b in F_p[x]/(g), coefficients low first, g monic with low coefficients gl; p < 2^31
    cdef long long t[2*DMAX]
    cdef long long c
    cdef int i,j
    for i in range(2*d-1):
        t[i]=0
    for i in range(d):
        c=a[i]
        if c:
            for j in range(d):
                t[i+j]=(t[i+j]+c*b[j])%p
    for i in range(2*d-2,d-1,-1):
        c=t[i]
        if c:
            for j in range(d):
                t[i-d+j]=(t[i-d+j]+(p-c)*gl[j])%p
    for i in range(d):
        out[i]=t[i]

cdef class _Fq:
    # The field F_p[x]/(g) for a prime p where g is irreducible; gives 1/sqrt of an element (start of the p-adic lift).
    cdef long long gl[DMAX]
    cdef long long zeta[DMAX]
    cdef long long p
    cdef int d,r
    cdef object s
    def __init__(self,g,p):
        cdef int i
        cdef long long z[DMAX]
        cdef long long t[DMAX]
        self.d=len(g)-1
        self.p=p
        if self.d>DMAX or p>=(1<<31):
            raise ValueError("degree or inert prime too large for the F_q code")
        for i in range(self.d):
            self.gl[i]=g[self.d-i]%p
        q=pow(p,self.d)
        self.s=q-1
        self.r=0
        while not self.s&1:
            self.s>>=1
            self.r+=1
        rnd=random.Random(12345)
        while True:                                 # a non-residue z, then zeta = z^s generates the 2-part
            for i in range(self.d):
                z[i]=rnd.randrange(p)
            self.powe(z,(q-1)>>1,t)
            t[0]=(t[0]+1)%p
            if self.is_zero(t):                     # z^((q-1)/2) = -1
                break
        self.powe(z,self.s,self.zeta)
    cdef int powe(self,const long long* a,object e,long long* out) except -1:
        cdef long long acc[DMAX]
        cdef long long base[DMAX]
        cdef int i
        for i in range(self.d):
            acc[i]=0
            base[i]=a[i]
        acc[0]=1
        while e:
            if e&1:
                _fq_mul(acc,base,self.gl,self.d,self.p,acc)
            e>>=1
            if e:
                _fq_mul(base,base,self.gl,self.d,self.p,base)
        for i in range(self.d):
            out[i]=acc[i]
        return 0
    cdef bint is_zero(self,const long long* a) noexcept:
        cdef int i
        for i in range(self.d):
            if a[i]:
                return False
        return True
    cdef bint is_one(self,const long long* a) noexcept:
        cdef int i
        if a[0]!=1:
            return False
        for i in range(1,self.d):
            if a[i]:
                return False
        return True
    def inv_sqrt(self,list poly):
        # poly: coefficients high first (any integers). Returns 1/sqrt(poly) mod (p, g) high first, or None if not a nonzero square.
        cdef long long nn[DMAX]
        cdef long long w[DMAX]
        cdef long long om[DMAX]
        cdef long long lb[DMAX]
        cdef long long tt[DMAX]
        cdef long long c[DMAX]
        cdef long long bb[DMAX]
        cdef long long x[DMAX]
        cdef int i,k,m,d=self.d
        cdef long long p=self.p
        for i in range(d):
            nn[i]=0
        k=len(poly)
        if k>d:
            return None
        for i in range(k):
            nn[i]=poly[k-1-i]%p
        if self.is_zero(nn):
            return None
        self.powe(nn,(self.s-1)>>1,w)               # w = n^((s-1)/2)
        _fq_mul(w,nn,self.gl,d,p,om)                # omega = n^((s+1)/2)
        _fq_mul(om,w,self.gl,d,p,lb)                # lambda = n^s
        for i in range(d):
            tt[i]=0
            c[i]=self.zeta[i]
        tt[0]=1                                     # tt = product of the corrections applied to omega
        m=self.r
        while not self.is_one(lb):                  # Tonelli-Shanks in F_q
            for i in range(d):
                x[i]=lb[i]
            k=0
            while not self.is_one(x):
                _fq_mul(x,x,self.gl,d,p,x)
                k+=1
                if k>=m:
                    return None                     # not a square
            for i in range(d):
                bb[i]=c[i]
            for i in range(m-k-1):
                _fq_mul(bb,bb,self.gl,d,p,bb)
            _fq_mul(om,bb,self.gl,d,p,om)
            _fq_mul(tt,bb,self.gl,d,p,tt)
            _fq_mul(bb,bb,self.gl,d,p,c)
            _fq_mul(lb,c,self.gl,d,p,lb)
            m=k
        # omega^2 = n and w^2*tt^2 = 1/n, so 1/omega = omega*w^2*tt^2
        _fq_mul(w,tt,self.gl,d,p,x)
        _fq_mul(x,x,self.gl,d,p,x)
        _fq_mul(x,om,self.gl,d,p,x)
        return [x[d-1-i] for i in range(d)]

@cython.boundscheck(False)
@cython.wraparound(False)
def _mulred(list x,list y,list g):
    # x*y reduced by the monic g, exact integer coefficients, high coefficient first (leading zeros stripped, like
    # div_poly(poly_prod(x,y),g)). One typed routine instead of two generic ones: this is called thousands of times
    # for every square root.
    cdef Py_ssize_t lx=len(x),ly=len(y),d=len(g)-1,i,j,n=lx+ly-1
    cdef list t=[0]*n
    cdef object c
    for i in range(lx):
        c=x[i]
        if c:
            for j in range(ly):
                t[i+j]=t[i+j]+c*y[j]
    for i in range(n-d):
        c=t[i]
        if c:
            for j in range(1,d+1):
                t[i+j]=t[i+j]-c*g[j]
    i=n-d if n>d else 0
    while i<n-1 and not t[i]:
        i+=1
    return t[i:]

def _mulmod(a,b,g,m):
    # a*b in Z[w]/(g) with coefficients mod m. g is monic with small coefficients, so the polynomial reduction is done
    # on the exact product (big * small) and each coefficient is reduced mod m only once: at these sizes one reduction
    # mod m costs about three multiplications.
    return [c%m for c in _mulred(a,b,g)]

def _lift_sqrt(root,gamma,g,p,bound):
    # Lift root = 1/sqrt(gamma) mod (p, g) and return sqrt(gamma) mod p^e with centred coefficients, p^e > 2*bound.
    # The inverse square root x is only lifted to half the target precision h = ceil(e/2) (Newton, the precision goes
    # 1, ..., h with each step at most doubling it and the last landing exactly on h). The square root itself then
    # comes from one correction step: y0 = gamma*x mod p^h, and y = y0 + p^h * (x * (gamma - y0^2)/p^h / 2) mod p^e.
    # The full-precision work is one squaring of a half-size y0 and one half-size product, instead of four full ones.
    d=len(g)-1
    e=int((bound.bit_length()+2)/math.log2(p))+2        # p^e > 2*bound
    h=(e+1)>>1
    chain=[h]
    while chain[-1]>1:
        chain.append((chain[-1]+1)>>1)
    # gamma reduced for every precision of the chain, from the top down (each from the one above it)
    ph=pow(p,h)
    red=[[c%ph for c in gamma]]
    mods=[ph]
    for k in chain[1:]:
        m=pow(p,k)
        red.append([c%m for c in red[-1]])
        mods.append(m)
    for lvl in range(len(chain)-2,-1,-1):
        modulo=mods[lvl]
        tmp=_mulmod(_mulmod(root,root,g,modulo),red[lvl],g,modulo)
        tmp=[-i%modulo for i in tmp]
        tmp=[0]*(d-len(tmp))+tmp
        tmp[-1]=(tmp[-1]+3)%modulo
        inv2=(modulo+1)>>1
        root=[inv2*i%modulo for i in _mulmod(root,tmp,g,modulo)]
    y0=_mulmod(root,red[0],g,ph)                        # root = 1/sqrt(gamma) mod p^h
    y0=[0]*(d-len(y0))+y0
    if e>h:
        pl=pow(p,e-h)
        sq=_mulred(y0,y0,g)                             # exact: gamma - y0^2 is divisible by p^h
        sq=[0]*(d-len(sq))+sq
        gg=[0]*(d-len(gamma))+list(gamma)
        t=[((gg[i]-sq[i])//ph)%pl for i in range(d)]
        corr=_mulmod(root,t,g,pl)
        inv2=(pl+1)>>1
        corr=[inv2*i%pl for i in corr]
        corr=[0]*(d-len(corr))+corr
        modulo=ph*pl
        y0=[(y0[i]+ph*corr[i])%modulo for i in range(d)]
    else:
        modulo=ph
    half=modulo>>1
    return [i-modulo if i>half else i for i in y0]

cdef inline int _ilog2(double x) noexcept nogil:
    # floor(log2 x) for x >= 1, 0 below
    cdef int e
    if x<1.0:
        return 0
    frexp(x,&e)
    return e-1

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef list _sieve_seg(unsigned short[::1] sg,unsigned short[::1] sa,Py_ssize_t seg_len,Py_ssize_t block,
                     long long[::1] pr,long long[::1] st_r,unsigned short[::1] lg_r,
                     long long[::1] pa,long long[::1] st_a,unsigned short[::1] lg_a,
                     double a0,double da,double b0,double db,double m1,double m0,double[::1] f,int d,
                     int init_g,int slack_g,int slack_a,int slack_lp,int slack_lp2):
    # Sieve one segment of a line, block by block (a block fits the cache). Position j is the pair
    # (a,b) = (a0 + j*da, b0 + j*db). sg / sa collect round(2*log2 p) for the primes dividing the rational value / the
    # norm; each prime's progression starts at st (counted from the start of the segment, st >= seg_len: never) and the
    # arrays st_r / st_a are advanced in place from block to block. A position is reported when both sums reach
    # 2*floor(log2 |value|) minus the slack of its side (slack_g rational, slack_a norm; in half bits). slack_lp is the
    # extra room for one large prime, granted to one side at a time; slack_lp2 is the (smaller) room each side gets when
    # both keep a large prime.
    # The exact test costs a polynomial evaluation, so positions are first screened 32 at a time against a threshold
    # that is a lower bound for the whole group (0 if a value changes sign inside it).
    cdef Py_ssize_t i,j,k,j0,j1,jm,base,blen,n_r=pr.shape[0],n_a=pa.shape[0]
    cdef long long p
    cdef unsigned short lg
    cdef double a,b,G,F,bp,Ga,Gb,F0,F1,F2,x
    cdef int eg,ef,t,tg,ta,which,dg,dn
    cdef list out=[]
    base=0
    while base<seg_len:
        blen=seg_len-base
        if blen>block:
            blen=block
        with nogil:
            for j in range(blen):
                sg[j]=init_g
                sa[j]=0
            for k in range(n_r):
                p=pr[k]
                lg=lg_r[k]
                i=st_r[k]-base
                while i<blen:
                    sg[i]+=lg
                    i+=p
                st_r[k]=i+base
            for k in range(n_a):
                p=pa[k]
                lg=lg_a[k]
                i=st_a[k]-base
                while i<blen:
                    sa[i]+=lg
                    i+=p
                st_a[k]=i+base
        j0=0
        while j0<blen:
            j1=j0+32
            if j1>blen:
                j1=blen
            # group thresholds
            Ga=(a0+(base+j0)*da)*m1-(b0+(base+j0)*db)*m0
            Gb=(a0+(base+j1-1)*da)*m1-(b0+(base+j1-1)*db)*m0
            if (Ga<0)!=(Gb<0):
                tg=-slack_g-slack_lp
            else:
                x=fabs(Ga)
                if fabs(Gb)<x:
                    x=fabs(Gb)
                tg=2*_ilog2(x)-slack_g-slack_lp
            jm=(j0+j1)>>1
            for which in range(3):
                j=j0 if which==0 else (jm if which==1 else j1-1)
                a=a0+(base+j)*da
                b=b0+(base+j)*db
                F=f[0]
                bp=1.0
                for t in range(1,d+1):
                    bp*=b
                    F=F*a+f[t]*bp
                if which==0:
                    F0=F
                elif which==1:
                    F1=F
                else:
                    F2=F
            if (F0<0)!=(F1<0) or (F1<0)!=(F2<0):
                ta=-slack_a-slack_lp
            else:
                x=fabs(F0)
                if fabs(F1)<x:
                    x=fabs(F1)
                if fabs(F2)<x:
                    x=fabs(F2)
                ta=2*(_ilog2(x)-2)-slack_a-slack_lp          # two bits of room for the variation inside the group
            for j in range(j0,j1):
                if <int>sg[j]<tg or <int>sa[j]<ta:
                    continue
                # exact test for this position
                a=a0+(base+j)*da
                b=b0+(base+j)*db
                eg=_ilog2(fabs(a*m1-b*m0))
                dg=2*eg-<int>sg[j]-slack_g              # how far the rational side is from its threshold (<= 0: reached)
                if dg>slack_lp:
                    continue
                F=f[0]
                bp=1.0
                for t in range(1,d+1):
                    bp*=b
                    F=F*a+f[t]*bp
                ef=_ilog2(fabs(F))
                dn=2*ef-<int>sa[j]-slack_a
                if dn<=0 or (dn<=slack_lp and dg<=0) or (dn<=slack_lp2 and dg<=slack_lp2):
                    out.append(base+j)
            j0=j1
        base+=blen
    return out


def _setup(n,d,alg_primes,poly=None):
    # one-time: polynomial, roots, characters, inert prime. poly = (f, m0, m1, characters) of an earlier set-up skips the
    # search: every process that sieves for the same matrix must use the same polynomial and the same characters
    t0=default_timer()
    primes=[2]+list(alg_primes)
    B=primes[-1]+1
    if poly is None:
        res=Kleinjung_poly_search(n,primes,NFS_NB_ROOTS,NFS_PRIME_BOUND,NFS_MULTIPLIER,int(pow(n,1/(d+1))),d,NFS_POLY_COARSE,NFS_POLY_PRECISE,None)
        if res is None:
            print("[i]NFS(d="+str(d)+"): the polynomial search found nothing for this degree and size")
            return None
        f_x,m0,m1=res[0],res[1],res[2]
        f_x,m0,M=evaluate_polynomial_quality(f_x,B,m0,m1,primes,None)
    else:
        f_x,m0,m1=list(poly[0]),poly[1],poly[2]
    assert eval_F(m0,m1,f_x,d)%n==0
    M,skew=get_sieve_region(f_x,B)
    leading=f_x[0]
    f_prime=get_derivative(f_x)
    g=[1,f_x[1]]
    for i in range(2,len(f_x)):
        g.append(f_x[i]*pow(leading,i-1))
    g_prime=get_derivative(g)
    st={'n':n,'d':d,'f':f_x,'m0':m0,'m1':m1,'M':M,'skew':skew,'leading':leading,'g':g,'primes':primes,
        'g_prime_sq':div_poly(poly_prod(g_prime,g_prime),g),
        'g_prime_eval':pow(leading,d-2,n)*eval_F(m0,m1,f_prime,d-1)%n if d>=2 else 1}
    R={}
    for p in primes:
        if p<=10*d or leading%p==0:
            # exhaustive: the upstream fast_roots() returns no roots when p divides the leading coefficient
            R[p]=[r for r in range(p) if eval_mod(f_x,r,p)==0]
        else:
            R[p]=sorted(set(find_roots_poly(f_x,p)))
    st['R']=R
    assert _factor_over(abs(leading),primes)[1]==1
    # (prime, root) pairs as arrays: p divides F(a,b) with b != 0 mod p exactly when a = b*r mod p for one of them
    st['ap']=np.array([p for p in primes for r in R[p]],dtype=np.int64)
    st['ar']=np.array([r for p in primes for r in R[p]],dtype=np.int64)
    st['apd']=st['ap'].astype(np.float64)               # the same in doubles, for _div_alg
    st['api']=1.0/st['apd']
    st['ard']=st['ar'].astype(np.float64)
    st['dbuf']=np.zeros(len(st['ap'])+1,dtype=np.int64)
    st['pp']=np.zeros(0,dtype=np.int64)                # the rational primes used so far (filled by nfs_sieve) and m0, m1 mod each
    st['div_lead']=[p for p in primes if leading%p==0]
    cq=[] if poly is None else list(poly[3])
    q=primes[-1]
    while poly is None and len(cq)<NFS_CHARS:
        q+=2 if q>2 else 1
        if not is_prime(q) or leading%q==0 or n%q==0:
            continue
        for r in find_roots_poly(f_x,q):
            if eval_mod(f_prime,r,q):
                cq.append((q,r))
    st['cq']=cq
    inert=0
    cand=[p for p in primes[::-1] if p>2]
    q=primes[-1]
    while inert==0 and len(cand)<len(primes)+2000:
        for p in cand:
            if p<(1<<31) and m1%p and leading%p and n%p and irreducibility(g,p):
                inert=p
                break
        cand=[]
        while len(cand)<200:
            q+=2
            if is_prime(q):
                cand.append(q)
        if q>50*primes[-1]:
            break
    if inert==0:
        print("[i]NFS(d="+str(d)+"): no prime found where f stays irreducible, so the lifting square root cannot run for this polynomial")
        return None
    st['inert']=inert
    st['fq']=_Fq(g,inert)
    # complex embeddings, for the size of the square roots: |coefficient j of beta| <= sum_i |L[j][i]| * |sigma_i(beta)|
    om=np.roots(np.array(g,dtype=np.float64))
    st['om']=om
    gp=np.polyval(np.array(get_derivative(g),dtype=np.float64),om)
    st['lg_gp']=np.log2(np.abs(gp))
    Lm=np.zeros((d,d))
    for i in range(d):
        num=np.poly(np.delete(om,i))[::-1]/gp[i]           # Lagrange polynomial of root i, low coefficient first
        Lm[:,i]=np.log2(np.abs(num)+1e-300)
    st['lg_L']=Lm
    f_norm=0
    tmp=1
    for x in g:
        f_norm+=x*x*tmp
        tmp*=leading
    st['f_norm']=int(math.sqrt(f_norm))+1
    st.update({'seen':set(),'free_done':set(),'line':0,'width':{},'rat_sieved':0,'rat_all':set([2]),'nrows':0,'counts':[0,0,0,0,0]})
    st['cbit']={qs:k for k,qs in enumerate(cq)}        # character -> bit; the two bits above them: sign of the norm, one per relation
    print("[i]NFS(d="+str(d)+"): f = "+str(f_x)+"   g = "+str(m1)+"*x - "+str(m0)+"   skew "+str(skew)+", sieve half-width "+str(M))
    print("[i]NFS(d="+str(d)+"): "+str(len(primes))+" algebraic primes with "+str(sum(len(R[p]) for p in primes))+" degree-1 ideals, "+str(len(cq))+" characters, inert prime "+str(inert)+"  (set-up %.1fs)"%(default_timer()-t0))
    return st

cdef extern from "math.h":
    double rint(double) nogil

@cython.boundscheck(False)
@cython.wraparound(False)
def _div_alg(double a,double b,double[::1] p,double[::1] pinv,double[::1] r,long long[::1] out):
    # which (p, r) have a = b*r mod p: the test in doubles (every product stays below 2^53, so it is exact)
    cdef Py_ssize_t k,m=0,n=p.shape[0]
    cdef double x
    with nogil:
        for k in range(n):
            x=a-b*r[k]
            if x-rint(x*pinv[k])*p[k]==0.0:
                out[m]=k
                m+=1
    return m

@cython.boundscheck(False)
@cython.wraparound(False)
def _div_rat(double a,double b,double[::1] p,double[::1] pinv,double[::1] m1p,double[::1] m0p,long long[::1] out):
    # which p divide a*m1 - b*m0, the same way: a is reduced mod p first so that the products stay exact
    cdef Py_ssize_t k,m=0,n=p.shape[0]
    cdef double x
    with nogil:
        for k in range(n):
            x=a-rint(a*pinv[k])*p[k]
            x=x*m1p[k]-b*m0p[k]
            if x-rint(x*pinv[k])*p[k]==0.0:
                out[m]=k
                m+=1
    return m

def _raw_row(st,a,b,G,Fv,lp_alg=0,lp_rat=0):
    # One (a,b) pair as a matrix row: (column numbers, character bits), or None if a value does not factor as expected.
    # lp_alg / lp_rat: a large prime known to divide the norm / the rational value exactly once.
    leading=st['leading']
    R=st['R']
    nq=len(st['cq'])
    ch=1<<(nq+1)                                        # even number of relations in every dependency
    cols=[]
    # the primes that can divide the norm: a = b*r mod p for a root r, or p | c_d and p | b
    fast=b!=0 and abs(a)<(1<<40) and b<(1<<20)
    if fast:
        m=_div_alg(a,b,st['apd'],st['api'],st['ard'],st['dbuf'])
        cand=set(st['ap'][st['dbuf'][:m]].tolist())
    else:
        cand=set(st['ap'][(a-b*st['ar'])%st['ap']==0].tolist())
    for p in st['div_lead']:
        if b%p==0:
            cand.add(p)
    ex,rest=_factor_over(Fv,sorted(cand))
    if rest!=(lp_alg if lp_alg else 1):
        return None
    for p,e in ex.items():
        if not e&1:
            continue
        if p==-1:
            ch|=1<<nq                                   # sign of the norm
            continue
        if leading%p==0 and b%p==0:
            cols.append(NFS_ALG+(p<<28)+p)              # the ideal above p at the projective root (p divides c_d)
        for r in R[p]:
            if (a-b*r)%p==0:
                cols.append(NFS_ALG+(p<<28)+r)          # degree-1 ideal (p, alpha - r)
    if lp_alg:
        if b%lp_alg==0:
            return None
        cols.append(NFS_ALG+(lp_alg<<28)+(a*pow(b,-1,lp_alg))%lp_alg)
    k=0
    for q,s in st['cq']:
        x=(a-b*s)%q
        if x and pow(x,(q-1)>>1,q)!=1:
            ch|=1<<k                                    # quadratic character
        k+=1
    pp=st['pp']
    if fast:
        m=_div_rat(a,b,st['ppd'],st['ppi'],st['m1pd'],st['m0pd'],st['dbuf'])
        rc=pp[st['dbuf'][:m]].tolist()
    else:
        rc=pp[(a%pp*st['m1p']-b%pp*st['m0p'])%pp==0].tolist() if b else pp.tolist()
    rex,rest=_factor_over(G,rc)
    if rest!=(lp_rat if lp_rat else 1):
        return None
    for p,e in rex.items():
        if e&1:
            cols.append(1 if p==-1 else p)
    if lp_rat:
        cols.append(lp_rat)
    return cols,ch

def nfs_params(n,degree):
    # what dep_square needs, as plain values (and what a second sieving process needs to use the same polynomial)
    st=_STATE.get((n,degree))
    if st is None:
        return None
    par={k:st[k] for k in ('n','d','leading','m0','m1','inert','lg_gp','lg_L','f_norm','om')}
    par['g']=[int(x) for x in st['g']]
    par['g_prime_sq']=[int(x) for x in st['g_prime_sq']]
    par['g_prime_eval']=int(st['g_prime_eval'])
    par['poly']=([int(x) for x in st['f']],int(st['m0']),int(st['m1']),list(st['cq']))
    par['nchar']=len(st['cq'])+2
    return par

def nfs_open(n,primeslist,degree,poly=None):
    # set-up without sieving (poly: see _setup); False if it cannot be done
    if (n,degree) not in _STATE:
        _STATE[(n,degree)]=_setup(n,degree,primeslist,poly)
    return _STATE[(n,degree)] is not None

def nfs_state(n,degree):
    # the part of the sieving state a restart needs: lines done, their widths, the free relations already made
    st=_STATE[(n,degree)]
    return {'line':st['line'],'width':dict(st['width']),'free_done':set(st['free_done']),'rat_sieved':st['rat_sieved']}

def nfs_restore(n,degree,saved,pairs):
    st=_STATE[(n,degree)]
    st['line']=saved['line']
    st['width']=dict(saved['width'])
    st['free_done']=set(saved['free_done'])
    st['rat_sieved']=saved['rat_sieved']
    st['seen'].update((int(a),int(b)) for a,b in pairs)

def nfs_sieve(n,primeslist,fb_keep,picks,degree,lines,T,force_T,force_lines,lp_bits=-1,lp2_bits=0,line_mod=1,line_res=0):
    # One round of sieving. primeslist: the algebraic factor base (odd primes not dividing n; fixed after the first call).
    # fb_keep + picks: the rational factor base, any odd primes: they need not be in the algebraic base.
    # picks: primes forced to divide the rational value on a sub-lattice.
    # lines / T: lines b added to the region by this call, and its half-width in a: max(T, skew * lines so far), with
    #   T = 0 standing for the width the polynomial was tuned for. Lines sieved before are widened when the width grows.
    # line_mod / line_res: this process only sieves the lines b = line_res mod line_mod (several processes share a region).
    # lp_bits: a relation may keep a prime up to 2^lp_bits outside the factor base (never above the square of the
    #   largest base prime, so it is certainly prime; 0 = none; -1 = the module default). lp2_bits: a relation may keep one
    #   on each side at once if both are below 2^lp2_bits (0 = never both).
    # Returns the new rows: (pairs as an int64 array [k,2], character bits (list of int), column numbers (int64, all rows
    # one after the other), where each row ends in that array (int64 [k]), number of rows without a large prime), or None.
    st=_STATE.get((n,degree))
    if st is None:
        if (n,degree) in _STATE:
            return None
        st=_setup(n,degree,primeslist)
        _STATE[(n,degree)]=st
        if st is None:
            return None
    t_start=default_timer()
    d,f_x,m0,m1,leading,primes,R=st['d'],st['f'],st['m0'],st['m1'],st['leading'],st['primes'],st['R']
    picks=[q for q in picks if q>2 and m1%q!=0 and n%q!=0]
    rat_primes=sorted(p for p in set(fb_keep)|set(picks) if p>2)
    fb_r=[2]+rat_primes
    if not st['rat_all'].issuperset(fb_r):
        st['rat_all'].update(fb_r)
        allr=sorted(st['rat_all'])
        st['pp']=np.array(allr,dtype=np.int64)
        st['m0p']=np.array([m0%p for p in allr],dtype=np.int64)
        st['m1p']=np.array([m1%p for p in allr],dtype=np.int64)
        st['ppd']=st['pp'].astype(np.float64)                # the same in doubles, for _div_rat
        st['ppi']=1.0/st['ppd']
        st['m0pd']=st['m0p'].astype(np.float64)
        st['m1pd']=st['m1p'].astype(np.float64)
        st['dbuf']=np.zeros(max(len(allr),len(st['ap']))+1,dtype=np.int64)
        st['big_r']=None
    if st.get('big_r') is None or st.get('big_r_n')!=len(fb_r):
        st['big_r']=_big(math.prod(fb_r))
        st['big_r_n']=len(fb_r)
        st['big_a']=_big(math.prod(primes))
    big_r=st['big_r']
    big_a=st['big_a']
    f_dbl=np.array(f_x,dtype=np.float64)
    o_ab=[]
    o_ch=[]
    o_cols=[]
    o_end=[]
    counts=[0,0,0,0,0]                  # free, unforced, forced, with one large prime, with two
    def emit(a,b,row,kind):
        o_ab.append((a,b))
        o_ch.append(row[1])
        o_cols.extend(row[0])
        o_end.append(len(o_cols))
        counts[kind]+=1
    # free relations: p splits completely, so (p) is the product of its d degree-1 ideals; rational side p*m1
    if line_res==0 and _is_smooth(m1,math.prod(fb_r)):
        for p in rat_primes:
            if p in st['free_done'] or p not in R or len(R[p])!=d or leading%p==0:
                continue
            st['free_done'].add(p)
            row=_raw_row(st,p,0,p*m1,leading*pow(p,d) if d&1 else leading*p)
            if row is not None:
                emit(p,0,row,0)
    # sieve set-up
    # the two sides have their own prime lists: the rational primes need not be in the algebraic base
    pickset=set(picks)
    sv_r=[p for p in rat_primes if p>NFS_SMALL and p not in pickset and m1%p!=0]
    sv_a=[p for p in primes if p>NFS_SMALL]
    lg2=lambda p:int(round(2*math.log2(p)))                 # logs in half bits
    key=(len(sv_r),len(picks))
    if st.get('sv_key')!=key:
        st['sv_key']=key
        st['sv']=(np.array(sv_r,dtype=np.int64),
                  np.array([(m0*pow(m1,-1,p))%p for p in sv_r],dtype=np.int64),     # a = b*m0/m1 mod p
                  np.array([lg2(p) for p in sv_r],dtype=np.uint16),
                  np.array([p for p in sv_a for r in R[p]],dtype=np.int64),
                  np.array([r for p in sv_a for r in R[p]],dtype=np.int64),
                  np.array([lg2(p) for p in sv_a for r in R[p]],dtype=np.uint16))
    r_p,r_root,r_lg,a_p,a_root,a_lg=st['sv']
    if lp_bits<0:
        lp_bits=NFS_LP_BITS
    lp_a=min(1<<lp_bits,primes[-1]*primes[-1]) if lp_bits>0 else 0
    # a rational cofactor is only certainly prime below the square of the smallest prime the rational base lacks
    if st.get('miss_n')!=len(fb_r):
        fbs=set(fb_r)
        st['miss_n']=len(fb_r)
        st['miss']=next((p for p in primes if p not in fbs),primes[-1])
    lp_r=min(lp_a,st['miss']*st['miss'])
    slack_g=2*NFS_SLACK
    slack_a=2*NFS_SLACK
    slack_lp=2*(lp_a.bit_length()-1) if lp_a else 0
    lp2=min(1<<lp2_bits,lp_a,lp_r) if lp2_bits>0 else 0
    slack_lp2=2*(lp2.bit_length()-1) if lp2 else -1000000
    ncand=0
    nhit=0
    _T.clear()
    seen=st['seen']
    def take_many(pairs,kind):
        # the candidates of a line: rational side first, the norm only for those that pass
        c=[]
        gs=[]
        for a,b in pairs:
            if b==0 or (a,b) in seen or math.gcd(a,b)!=1:
                continue
            G=a*m1-b*m0
            if G:
                c.append((a,b))
                gs.append(G)
        if not c:
            return
        cgs=_cofactors(gs,big_r)
        c2=[]
        fs=[]
        for k in range(len(c)):
            cg=cgs[k]
            if cg!=1 and cg>lp_r:
                continue
            a,b=c[k]
            Fv=eval_F(a,b,f_x,d)
            if Fv:
                c2.append((a,b,gs[k],cg))
                fs.append(Fv)
        if not c2:
            return
        cfs=_cofactors(fs,big_a)
        for k in range(len(c2)):
            a,b,G,cg=c2[k]
            cf=cfs[k]
            if cf!=1 and (cf>lp_a or (cg!=1 and (cf>lp2 or cg>lp2))):
                continue
            row=_raw_row(st,a,b,G,fs[k],cf if cf!=1 else 0,cg if cg!=1 else 0)
            if row is None:
                continue
            seen.add((a,b))
            emit(a,b,row,kind if cf==1 and cg==1 else (3 if cf==1 or cg==1 else 4))
    # the region: the rectangle |a| <= M, 1 <= b <= L. It grows by `lines` lines per call; M follows skew*L once that
    # exceeds the starting width. Every pair is sieved once: for a line seen before, only the part outside its old width.
    L=st['line']+lines
    st['line']=L
    M=max(T if T>0 else st['M'],st['skew']*L)
    width=st['width']
    if len(fb_r)>1.05*st['rat_sieved']:
        # the rational factor base has grown by more than 5% since the region was sieved: pairs that failed only because
        # of a missing prime are worth another look, so the whole region is sieved again (known pairs are skipped)
        width.clear()
        st['rat_sieved']=len(fb_r)
    FT=(force_T if force_T>0 else primes[-1]) if len(picks)>0 else 0
    sg=np.zeros(NFS_BLOCK,dtype=np.uint16)
    sa=np.zeros(NFS_BLOCK,dtype=np.uint16)
    order=st.get('order')
    if order is None:
        # the lines in the order they are worth: on line b only the a coprime to b count, and a b with small prime
        # factors also makes the norm less likely to be smooth, so such lines wait until the others have grown
        bm=200000
        phi=np.arange(bm+1,dtype=np.float64)
        for p in range(2,bm+1):
            if phi[p]==p:
                phi[p::p]*=(1.0-1.0/p)
        bb=np.arange(1,bm+1,dtype=np.float64)
        order=st['order']=(np.argsort(bb*(bb/phi[1:])**(NFS_LINE_EXP/100.0),kind='stable')+1).tolist()
    for li in range(L):
        if li%line_mod!=line_res%line_mod:
            continue
        b=order[li]
        w=width.get(b,-1)
        if w>=M:
            continue
        width[b]=M
        st_r=(b*r_root-(-M))%r_p                           # where each progression starts, counted from a = -M
        st_a=np.where(b%a_p!=0,(b*a_root+M)%a_p,-1)        # -1: p | b, the norm is not divisible unless p | c_d (left to the slack)
        for lo,hi in ([(-M,M)] if w<0 else [(-M,-w-1),(w+1,M)]):
            sv_len=hi-lo+1
            sh=lo+M
            ncand+=sv_len
            t1=default_timer()
            hits=_sieve_seg(sg,sa,sv_len,NFS_BLOCK,r_p,(st_r-sh)%r_p,r_lg,a_p,np.where(st_a>=0,(st_a-sh)%a_p,sv_len),a_lg,
                            float(lo),1.0,float(b),0.0,float(m1),float(m0),f_dbl,d,0,slack_g,slack_a,slack_lp,slack_lp2)
            _tick('sieve',t1)
            t1=default_timer()
            nhit+=len(hits)
            take_many([(lo+j,b) for j in hits],1)
            _tick('candidates',t1)
    # forced lattice: the pairs with a = b*m0/m1 mod Qf, Qf = product of the picks, so Qf divides a*m1 - b*m0
    Qf=math.prod(picks)
    if Qf>1:
        FL=force_lines
        mm=(m0*pow(m1,-1,Qf))%Qf
        # Gauss reduction of the basis (Qf,0), (mm,1) for the norm a^2/s + b^2*s, s = the skew of the polynomial
        sk=float(st['skew'])
        nq=lambda v:(v[0]*v[0])/sk+(v[1]*v[1])*sk
        bq=lambda v,t:(v[0]*t[0])/sk+(v[1]*t[1])*sk
        b1=(Qf,0)
        b2=(mm,1)
        if nq(b1)<nq(b2):
            b1,b2=b2,b1
        while True:
            mu=round(bq(b1,b2)/nq(b2))
            b1=(b1[0]-mu*b2[0],b1[1]-mu*b2[1])
            if nq(b1)>=nq(b2):
                break
            b1,b2=b2,b1
        e1,e2=b2,b1
        G1=e1[0]*m1-e1[1]*m0
        G2=e2[0]*m1-e2[1]*m0
        # in lattice coordinates: rational value divisible by p when il = rr*jl, norm when il = t*jl (mod p)
        l_rp=[]
        l_rr=[]
        l_ap=[]
        l_ar=[]
        for p in sv_r:
            if G1%p!=0:
                l_rp.append(p)
                l_rr.append((-G2*pow(G1,-1,p))%p)
        for p in sv_a:
            for r in R[p]:
                A1=(e1[0]-r*e1[1])%p
                A2=(e2[0]-r*e2[1])%p
                if A1!=0:
                    l_ap.append(p)
                    l_ar.append((-A2*pow(A1,-1,p))%p)
        l_rp=np.array(l_rp,dtype=np.int64)
        l_rr=np.array(l_rr,dtype=np.int64)
        l_rlg=np.array([lg2(p) for p in l_rp.tolist()],dtype=np.uint16)
        l_ap=np.array(l_ap,dtype=np.int64)
        l_ar=np.array(l_ar,dtype=np.int64)
        l_alg=np.array([lg2(p) for p in l_ap.tolist()],dtype=np.uint16)
        fl_len=2*FT+1
        init_g=int(round(2*math.log2(Qf)))
        for jl in range(1,FL+1):
            ncand+=fl_len
            fp=[]
            for j in _sieve_seg(sg,sa,fl_len,NFS_BLOCK,l_rp,(l_rr*jl+FT)%l_rp,l_rlg,l_ap,(l_ar*jl+FT)%l_ap,l_alg,
                                float(-FT*e1[0]+jl*e2[0]),float(e1[0]),float(-FT*e1[1]+jl*e2[1]),float(e1[1]),
                                float(m1),float(m0),f_dbl,d,init_g,slack_g,slack_a,slack_lp,slack_lp2):
                il=j-FT
                a=il*e1[0]+jl*e2[0]
                b=il*e1[1]+jl*e2[1]
                if b<0:
                    a,b=-a,-b
                fp.append((a,b))
            take_many(fp,2)
    total=default_timer()-t_start
    tot=st['counts']
    for k in range(5):
        tot[k]+=counts[k]
    st['nrows']+=len(o_ab)
    print("[i]NFS(d="+str(d)+")"+(" #"+str(line_res) if line_mod>1 else "")+": region |a|<="+str(M)+", "+str(L)+" lines (b up to "+str(max(order[:L]))+"), "+str(ncand)+" new pairs, "+str(nhit)+" candidates; rational side "+str(len(fb_r))+" primes"+(" (forced "+str(picks)+")" if picks else "")+"; %.1fs: "%total+", ".join(k+" %.1f"%v for k,v in _T.items()))
    print("[i]NFS(d="+str(d)+")"+(" #"+str(line_res) if line_mod>1 else "")+": new rows "+str(counts[0])+" free + "+str(counts[1]+counts[2])+" full"+(" ("+str(counts[2])+" forced)" if picks else "")+" + "+str(counts[3])+" with one large prime + "+str(counts[4])+" with two; so far "+str(tot[0]+tot[1]+tot[2])+" + "+str(tot[3])+" + "+str(tot[4]))
    if not o_ab:
        return None
    return (np.array(o_ab,dtype=np.int64),o_ch,np.array(o_cols,dtype=np.int64),np.array(o_end,dtype=np.int64),counts[0]+counts[1]+counts[2])

def dep_square(par,pairs):
    # A set of (a,b) pairs whose algebraic columns and characters cancel: returns (P, X) with X^2 = P mod n, P the exact
    # product of their rational values a*m1 - b*m0, or None.
    n,d,leading,m0,m1=par['n'],par['d'],par['leading'],par['m0'],par['m1']
    g=[_big(x) for x in par['g']]
    S=len(pairs)
    if S==0:
        return 1,1
    if S&1:
        return None
    xr=par['g_prime_eval']*pow(leading,S>>1,n)%n
    if math.gcd(xr,n)!=1:
        return None
    om=par['om']
    emb=np.zeros(d)
    polys=[]
    gs=[]
    u=1
    for a,b in pairs:
        a=int(a)
        b=int(b)
        polys.append([_big(-b),_big(a*leading)] if b else [_big(a*leading)])
        gs.append(_big(a*m1-b*m0))
        emb=emb+(np.log2(np.abs(a*leading-b*om)) if b else math.log2(abs(a*leading)))
        if abs(a)>u:
            u=abs(a)
        if b>u:
            u=b
    while len(polys)>1:
        nxt=[_mulred(polys[i],polys[i+1],g) for i in range(0,len(polys)-1,2)]
        if len(polys)&1:
            nxt.append(polys[-1])
        polys=nxt
    while len(gs)>1:
        nxt=[gs[i]*gs[i+1] for i in range(0,len(gs)-1,2)]
        if len(gs)&1:
            nxt.append(gs[-1])
        gs=nxt
    PG=int(gs[0])
    gamma=_mulred([_big(x) for x in par['g_prime_sq']],polys[0],g)
    pin=par['inert']
    fq=par.get('_fq')
    if fq is None:
        fq=par['_fq']=_Fq(par['g'],pin)
    root0=fq.inv_sqrt([int(x%pin) for x in gamma])
    if root0 is None:
        return None
    # size of the square root from its embeddings (24 bits of margin); if that fails, the coarse coefficient bound
    lg=par['lg_gp']+0.5*emb
    bits=float(np.max(lg[None,:]+par['lg_L']))+math.log2(d)
    for attempt in range(2):
        if attempt==0:
            bound=1<<max(1,int(bits)+25)
        else:
            fd=int(pow(d,1.5))+1
            bound=max(fd*pow(par['f_norm'],d-1-i)*pow(2*abs(leading)*(2*u)*par['f_norm'],S>>1) for i in range(d))
        root=_lift_sqrt([_big(x) for x in root0],gamma,g,_big(pin),bound)
        root=[0]*(d-len(root))+root
        y=eval_F(leading*m0,m1,root,d-1)*pow(m1,S>>1,n)%n
        X=int(y)*pow(xr,-1,n)%n
        if (X*X-PG)%n==0:
            return PG,X
    return None
