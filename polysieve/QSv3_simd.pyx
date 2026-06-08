#!python
#cython: language_level=3
# cython: profile=False
# cython: overflowcheck=False
###Author: Essbee Vanhoutte
###WORK IN PROGRESS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
###Improved QS Variant 

##References: I have borrowed many of the optimizations from here: https://stackoverflow.com/questions/79330304/optimizing-sieving-code-in-the-self-initializing-quadratic-sieve-for-pypy


###To build: python3 setup.py build_ext --inplace

from datetime import datetime
import random
import sympy
from itertools import chain
import itertools
import sys
import argparse
import multiprocessing
import time
import copy
from timeit import default_timer
import math
import gc
from cpython cimport array
import array
cimport cython
import numpy as np
from numpy.polynomial.legendre import leggauss

key=0                 #Define a custom modulus to factor
build_workers=8
keysize=150           #Generate a random modulus of specified bit length
workers=1 #max amount of parallel processes to use
quad_co_per_worker=1 #Amount of quadratic coefficients to check. Keep as small as possible.
base=1_000
small_base=1000
qbase=10
quad_sieve_size=10
g_debug=0 #0 = No debug, 1 = Debug, 2 = A lot of debug
g_lift_lim=0.5
thresvar=50  ##Log value base 2 for when to check smooths with trial factorization. Eventually when we fix all the bugs we should be able to furhter lower this.
lp_multiplier=2
min_prime=1
g_enable_custom_factors=0
g_p=107
g_q=41
mod_mul=0.25
g_max_exp=10
lin_sieve_size=1000
max_diff=1_000_000
degree=2
##Key gen function ##
def power(x, y, p):
    res = 1;
    x = x % p;
    while (y > 0):
        if (y & 1):
            res = (res * x) % p;
        y = y>>1; # y = y/2
        x = (x * x) % p;
    return res;

def power2(poly, f, p, exp):
    if exp == 1: return poly

    tmp = power2(poly, f, p, exp>>1)
    tmp = poly_prod(tmp, tmp)

    if exp&1:
        tmp = poly_prod(tmp, poly)
        return div_poly_mod(tmp, f, p)
    
    else: return div_poly_mod(tmp, f, p)

def miillerTest(d, n):
    a = 2 + random.randint(1, n - 4);
    x = power(a, d, n);
    if (x == 1 or x == n - 1):
        return True;
    while (d != n - 1):
        x = (x * x) % n;
        d *= 2;
        if (x == 1):
            return False;
        if (x == n - 1):
            return True;
    # Return composite
    return False;

def isPrime( n, k):
    if (n <= 1 or n == 4):
        return False;
    if (n <= 3):
        return True;
    d = n - 1;
    while (d % 2 == 0):
        d //= 2;
    for i in range(k):
        if (miillerTest(d, n) == False):
            return False;
    return True;

def generateLargePrime(keysize = 1024):
    while True:
        num = random.randrange(2**(keysize-1), 2**(keysize))
        if isPrime(num,4):
            return num

def findModInverse(a, m):
    if gcd(a, m) != 1:
        return None
    u1, u2, u3 = 1, 0, a
    v1, v2, v3 = 0, 1, m
    while v3 != 0:
        q = u3 // v3
        v1, v2, v3, u1, u2, u3 = (u1 - q * v1), (u2 - q * v2), (u3 - q * v3), v1, v2, v3
    return u1 % m
   

def generateKey(keySize):
    while True:
        p = generateLargePrime(keySize)
        print("[i]Prime p: "+str(p))
        q=p
        while q==p:
            q = generateLargePrime(keySize)
        print("[i]Prime q: "+str(q))
        n = p * q
        print("[i]Modulus (p*q): "+str(n))
        count=65537
        e =count
        if gcd(e, (p - 1) * (q - 1)) == 1:
            break

    phi=(p - 1) * (q - 1)
    d = findModInverse(e, (p - 1) * (q - 1))
    publicKey = (n, e)
    privateKey = (n, d)
    print('[i]Public key - modulus: '+str(publicKey[0])+' public exponent: '+str(publicKey[1]))
    print('[i]Private key - modulus: '+str(privateKey[0])+' private exponent: '+str(privateKey[1]))
    return (publicKey, privateKey,phi,p,q)
##END KEY GEN##




cdef modinv(n,p):
    p2=p
    n = n % p
    x =0
    u = 1
    while n:
        x, u = u, x - (p // n) * u
        p, n = n, p % n
    return x%p2

def bitlen(int_type):
    int_type=abs(int_type)
    length=0
    while(int_type):
        int_type>>=1
        length+=1
    return length   

def gcd(a,b): # Euclid's algorithm ##To do: Use a version without recursion?
    if b == 0:
        return a
    elif a >= b:
        return gcd(b,a % b)
    else:
        return gcd(b,a)
        
def QS(n,factor_list,sm,xlist,flist):#,jsymbols,testl,primeslist2,disc1_squared_list):#,disc_sr_list,pval_list,pflist):
    g_max_smooths=base+2#+qbase
    if len(sm) > g_max_smooths*10000000: 
        del sm[g_max_smooths:]
        del xlist[g_max_smooths:]
        del flist[g_max_smooths:]  
    M2 = build_matrix(factor_list, sm, flist)#,pflist)
    null_space=solve_bits(M2,factor_list,len(sm))
    f1,f2=extract_factors(n, sm, xlist, null_space)#,disc_sr_list,pval_list,pflist)
    if f1 != 0:
        print("[SUCCESS]Factors are: "+str(f1)+" and "+str(f2))
        return f1,f2   
    print("[FAILURE]No factors found")
    return 0,0

def extract_factors(N, relations, roots, null_space):#,disc_sr_list,pval_list,pflist):
    n = len(relations)
    for vector in null_space:
        prod_left = 1
        prod_right = 1
        pval=1
        disc_sr=1
        for idx in range(len(relations)):
            bit = vector & 1
            vector = vector >> 1
            if bit == 1:
                prod_left *= roots[idx]
                prod_right *= relations[idx]
               # print("adding to smooth:  "+str(relations[idx]))
            idx += 1
        sqrt_right = math.isqrt(prod_right)
        sqrt_left = prod_left
        ###Debug shit, remove for final version
        sqr1=prod_left**2%N 
        sqr2=prod_right%N
        if sqrt_right**2 != prod_right:
            print("not a square in the integers")
            time.sleep(10000)

        if sqr1 != sqr2:
            print("ERROR ERROR")
            time.sleep(10000)
        ###End debug shit#########
        sqrt_left = sqrt_left % N
        sqrt_right = sqrt_right % N
        factor_candidate = gcd(N, abs(sqrt_right+sqrt_left))


        if factor_candidate not in (1, N):
         #   print(str(factor_candidate)+" sm: "+str(sqrt_right)+" root: "+str(sqrt_left))
            other_factor = N // factor_candidate
            return factor_candidate, other_factor

    return 0, 0

def solve_bits(matrix,factor_base,length):
    n=length#len(factor_base)*1#base+2
    lsmap = {lsb: 1 << lsb for lsb in range(n+10000)}
    m = len(matrix)
    marks = []
    cur = -1
    mark_mask = 0
    for row in matrix:
        if cur % 100 == 0:
            print("", end=f"{cur, m}\r")
        cur += 1
        lsb = (row & -row).bit_length() - 1
        if lsb == -1:
            continue
        marks.append(n - lsb - 1)
        mark_mask |= 1 << lsb
        for i in range(m):
            if matrix[i] & lsmap[lsb] and i != cur:
                matrix[i] ^= row
    marks.sort()
    # NULL SPACE EXTRACTION
    nulls = []
    free_cols = [col for col in range(n) if col not in marks]
    k = 0
    for col in free_cols:
        shift2 = n - col - 1
        val = 1 << shift2
        fin = val
        for v in matrix:
            if v & val:
                fin |= v & mark_mask
        nulls.append(fin)
        k += 1
        if k == 10000000000: 
            break
    return nulls

def build_matrix(factor_base, smooth_nums, factors):#,pflist):
    fb_map = {val: i for i, val in enumerate(factor_base)}

    ind=1

    M2=[0]*(base+2)#+2+qbase)
    for i in range(len(smooth_nums)):
        for fac in factors[i]:
            idx = fb_map[fac]
            M2[idx] |= ind
        ind = ind + ind       
    return M2

@cython.profile(False)
def launch(n,primeslist1,primeslist2,small_primeslist):
    manager=multiprocessing.Manager()
    return_dict=manager.dict()
    jobs=[]
    start= default_timer()
    print("[i]Creating iN datastructure... this can take a while...")
    primeslist1c=copy.deepcopy(primeslist1)
    plists=[]


    duration = default_timer() - start
    print("[i]Creating iN datastructure in total took: "+str(duration))

    print("[*]Launching attack with "+str(workers)+" workers\n")
    find_comb(n,primeslist1,primeslist2,small_primeslist)

    return 

def get_gray_code(n):
    gray = [0] * (1 << (n - 1))
    gray[0] = (0, 0)
    for i in range(1, 1 << (n - 1)):
        v = 1
        j = i
        while (j & 1) == 0:
            v += 1
            j >>= 1
        tmp = i + ((1 << v) - 1)
        tmp >>= v
        if (tmp & 1) == 1:
            gray[i] = (v - 1, -1)
        else:
            gray[i] = (v - 1, 1)
    return gray

cdef tonelli(long long n, long long p):  # tonelli-shanks to solve modular square root: x^2 = n (mod p)
    if jacobi(n,p)!=1:
        print("error jacobi")
        return 1
    cdef long long q = p - 1
    cdef long long s = 0
    cdef long long z,c,r,t,m,t2,b,i

    while q % 2 == 0:
        q //= 2
        s += 1
    if s == 1:
        r = pow(n, (p + 1) // 4, p)
        return r
    for z in range(2, p):
        if -1 == jacobi(z, p):
            break
    c = pow(z, q, p)
    r = pow(n, (q + 1) // 2, p)
    t = pow(n, q, p)
    m = s
    t2 = 0
    while (t - 1) % p != 0:
        t2 = (t * t) % p
        for i in range(1, m):
            if (t2 - 1) % p == 0:
                break
            t2 = (t2 * t2) % p
        b = pow(c, 1 << (m - i - 1), p)
        r = (r * b) % p
        c = (b * b) % p
        t = (t * c) % p
        m = i

    return r

@cython.profile(False)
def jacobi(a, n):
    t=1
    while a !=0:
        while a%2==0:
            a //=2
            r=n%8
            if r == 3 or r == 5:
                t = -t
        a, n = n, a
        if a % 4 == n % 4 == 3:
            t = -t
        a %= n
    if n == 1:
        return t
    else:
        return 0    

@cython.boundscheck(False)
@cython.wraparound(False)
cdef factorise_fast(value,long long [::1] factor_base):
    if value == 0:
        return None,None,None,None
    seen_primes=[]
    seen_primes_indexes=[]
    factors = set()
    if value < 0:
        seen_primes.append(-1)
        seen_primes_indexes.append(-1)
        factors ^= {-1}
        value = -value
    while value % 2 == 0:
        seen_primes.append(2)
        seen_primes_indexes.append(-1)
        factors ^= {2}
        value //= 2
    #cdef int factor 
    cdef Py_ssize_t length=factor_base[0]#len(factor_base)#factor_base[0]
    cdef Py_ssize_t i=1
    while i < length:
        factor=factor_base[i]
       # print(factor)
        while value % factor == 0:
            seen_primes.append(factor)
            seen_primes_indexes.append(i-1)
            factors ^= {factor}
            value //= factor
        i+=1
    return factors, value,seen_primes,seen_primes_indexes



#@cython.boundscheck(False)
#@cython.wraparound(False)
def generate_modulus(n,primeslist,seen,tnum,close_range,too_close,LOWER_BOUND_SIQS,UPPER_BOUND_SIQS,tnum_bit,quad):
    const_1=1_000
    const_2=1_000_000

    small_B = base#len(primeslist)
    lower_polypool_index = 2
    upper_polypool_index = small_B - 1
    poly_low_found = False
    
    for i in range(small_B):  ##To do: Can be moved outside mainloop
        if primeslist[i] > LOWER_BOUND_SIQS and not poly_low_found:
            lower_polypool_index = i
            poly_low_found = True
            break
        if primeslist[i] > UPPER_BOUND_SIQS:
            upper_polypool_index = i - 1
            break
    small_B=upper_polypool_index
    counter4=0
    while counter4 < const_1:
        counter4+=1
        cmod = 1
        cfact = []#[0]*base
        indexes=[]
        counter2=0
        while counter2 < const_2:
            found_a_factor = False
            counter=0
            while(found_a_factor == False) and counter < const_2:
                randindex = random.randint(lower_polypool_index, upper_polypool_index)
                if  jacobi((-quad*n)%primeslist[randindex],primeslist[randindex])!=1:
                    counter+=1
                    continue
                potential_a_factor = primeslist[randindex]
                found_a_factor = True
                it=0
                length=len(cfact)
                while it < length:
                    if potential_a_factor ==cfact[it]:
                        found_a_factor = False
                        break
                    it+=1
                counter+=1
            if counter == const_2:
                cmod = 1
                s = 0
                cfact = []#[0]*base
                indexes=[]
                continue                
            cmod = cmod * potential_a_factor
            cfact.append(potential_a_factor)
            if  jacobi((-quad*n)%primeslist[randindex],primeslist[randindex])!=1:
                print("THE FUC")
                time.sleep(1000000)
            indexes.append(randindex)
            j = tnum_bit - cmod.bit_length()
            counter2+=1
            if j < too_close:
                cmod = 1
                s = 0
                cfact = []#[0]*base
                indexes=[]
                continue
            elif j < (too_close + close_range):
                break
        a1 = tnum // cmod
        mindiff = 100000000000000000
        randindex = 0
        for i in range(small_B):
            if abs(a1 - primeslist[i]) < mindiff:
                randindex = i
                mindiff = abs(a1 - primeslist[i])
                
        

        found_a_factor = False
        counter3=0
        while not found_a_factor and counter3< const_2 and randindex <base:
            if  jacobi((-quad*n)%primeslist[randindex],primeslist[randindex])!=1:
                randindex += 1
                continue
            potential_a_factor = primeslist[randindex]

            found_a_factor = True
            it=0
            length=len(cfact)
            while it < length:
                if potential_a_factor ==cfact[it]:
                    found_a_factor = False
                    break
                it+=1
            if not found_a_factor:
                randindex += 1
            counter3+=1
        if randindex > small_B:
            continue
        if counter3==const_2:
            continue

        cmod = cmod * potential_a_factor
        if  jacobi((-quad*n)%primeslist[randindex],primeslist[randindex])!=1:
            print("THE FUC: ",randindex)
            time.sleep(1000000)
        cfact.append(potential_a_factor)
        indexes.append(randindex)

        diff_bits = (tnum - cmod).bit_length()
        if diff_bits < tnum_bit:
            if cmod in seen:
                continue
            else:
                seen.append(cmod)
                return cmod,cfact,indexes
    return 0,0,0

def solve_lin_con(a,b,m):
    ##ax=b mod m
    #g=gcd(a,m)
    #a,b,m = a//g,b//g,m//g
    return pow(a,-1,m)*b%m  





def poly_search(n, primes, nb_roots, prime_bound, c, M, d, NB_POLY_COARSE_EVAL, NB_POLY_PRECISE_EVAL):
    print("Starting polynomial search")


    return Kleinjung_poly_search(n, primes, nb_roots, prime_bound, c, M, d, NB_POLY_COARSE_EVAL,
                                     NB_POLY_PRECISE_EVAL)


# Kleinjung first polynomial search algorithm
def Kleinjung_poly_search(n, primes, NB_ROOTS, PRIME_BOUND, MULTIPLIER, M, d, NB_POLY_COARSE_EVAL, NB_POLY_PRECISE_EVAL):
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
      #  print("roots: ",roots)
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
                    for j in range(d): 
                        root_used[i][j] = root_used[i][j]*tmp2%prod

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
                       # print("root_used: ",root_used)
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
                            print("")
                          #  print("Polynomial search done in "+format_duration(t2-t1)+".\n")
                            print(str(cpt)+" polynomials created, "+str(len(polys))+" kept for ranking")
                            print("Average L2 score = "+str(avg/cpt))
                            print("Ranking polynomials")

                            return select_best_poly_candidate(polys, primes)

                        z += 1
        a_d += MULTIPLIER

def poly_prod(a, b):
    res = [0]*(max(len(a), len(b))+min(len(a), len(b))-1)

    for i in range(len(a)):
        for j in range(len(b)):
            res[i+j] += a[i]*b[j]

    return res

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

def shift(poly,k):
    res = [i for i in poly]
    for i in range(len(poly)-1):
        for j in range(i+1, len(poly)):
            res[j] += binom(j-i, len(poly)-1-i)*pow(k, j-i)*poly[i]
    return res

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

def binom(k, n):
    res = 1

    for i in range(n-k+1, n+1): res *= i
    for i in range(2, k+1): res //= i

    return res

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

def compute_f(n, a_d, m0, d, prod, root_used, NB_ROOTS, e):
    f0 = (n-a_d*pow(m0, d))/(prod*prod*pow(m0, d-1))
    f = []

    for i in range(NB_ROOTS):
        line = [0]*d
        for j in range(d):
            line[j] = -(a_d*d*root_used[i][j]/pow(prod, 2)+e[i][j]/prod)
        f.append(line)

    return f, f0

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

def invmod(a, m):
    (r, u, R, U) = (a, 1, m, 0)
    while R:
        q = r//R
        (r, u, R, U) = (R, U, r - q *R, u - q*U)
    return u%m

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

def gcd_mod(f, poly, p):
    while poly != [0]*len(poly):
        (f,poly) = (poly, div_poly_mod(f, poly, p))

    return f

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
        h = power2([1, a], g, p, (p-1)>>1)
        for k in range(len(h)):
            if h[k]:
                h = h[k:]
                break
        h[-1] -= 1
        h = gcd_mod(h, g, p)
    r = roots(h, p)
    h = quotient_poly_mod(g, h, p)
    return r+roots(h, p)
    

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

# Compute x such that x² = n (mod p)
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

def eval_F(x, y, f, d):
    tmp = 0
    tmp2 = 1

    for k in range(d):
        tmp += f[k]*tmp2
        tmp *= x
        tmp2 *= y

    tmp += f[d]*tmp2

    return tmp

def dickman(x, table):
    k,res = math.ceil(x),0
    delta = k-x
    if k-1 > len(table): table = get_dickman_table(k)
    tmp = 1
    for i in range(len(table[k-2])):
        res += table[k-2][i]*tmp
        tmp *= delta

    return res, table

def non_central(k, l ,x):
    if x <= 0: return 0
    res = 0
    for i in range(100): res += central(k+2*i, x)*pow(l/2, i)/fac(i)
    
    return res*math.exp(-l/2)

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

def select_best_poly_candidate(polys, primes):
  #  print("polys: ",polys)
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

def fac(n):
    res = 1
    for i in range(2, n+1): res *= i
    return res

def central(k,x): return pow(x, k/2-1)/(math.exp(x/2)*pow(2, k/2)*math.gamma(k/2))

def eval_mod(f, x, n):
    res = 0

    for i in range(len(f)-1):
        res += f[i]
        res *= x
        res %= n

    res += f[-1]
    res %= n

    return res

def get_derivative(f):
    res = [0]*(len(f)-1)
    for i in range(len(f)-1):
        res[i] = (len(f)-1-i)*f[i]
    return res

def evaluate(f, x):
    x = int(x)
    res = 0
    for coeff in f:
        res = res * x + int(coeff)
    return res

def get_partials(mod,list1):
    i=0
    new_list=[]
    while i < len(list1):
        prime=list1[i]
        new_list.append(prime)
        new_list.append([])
        k=0
        while k < len(list1[i+1]):
            r1=list1[i+1][k]
            aq = mod // prime
            invaq = modinv(aq%prime, prime)
            gamma = r1 * invaq % prime
            new_list[-1].append(aq*gamma)
           # lin+=aq*gamma
           # all_lin_parts.append(aq*gamma)
            k+=1
        i+=2
    

    return new_list

def f_minimum_x(f_x):
    return -f_x[1] / (2 * f_x[0])

def make_polys(n, y, offset, b):
    n, y, offset, b = int(n), int(y), int(offset), int(b)
    c = -(n - y * offset + offset ** 2) * b ** 2
    f_x = [1, int(y * b), int(c)]
    g_x = [1, int((y - offset) * b)]
    h_x = [1, int(offset * b)]
    return f_x, g_x, h_x

def sieve_interval(f_x, g_x, h_x, M):
    num = -f_x[1]
    den = 2 * f_x[0]
    return num, den, M 

def int_log(n):
    if n <= 0:
        raise ValueError
    return abs(n).bit_length() * 0.6931471805599453  

def e_score(n, y, offset, b, B, M, n_gauss=32, table=None):
    if table is None:
        table = get_dickman_table(20)

    f_x, g_x, h_x = make_polys(n, y, offset, b)
    center_num, center_den, half_width = sieve_interval(f_x, g_x, h_x, M)
    log_B = math.log(B)

    nodes, weights = leggauss(n_gauss)
    scale = half_width 

    score = 0.0
    for xi, wi in zip(nodes, weights):
        x = int(round(center_num / center_den + half_width * xi))
        fval = evaluate(f_x, x)
        ghval = evaluate(g_x, x) * evaluate(h_x, x)

        abs_fval = abs(fval)
        abs_ghval = abs(ghval)

        if abs_fval == 0 or abs_ghval == 0:
            continue

        try:
            log_f  = int_log(abs_fval)
            log_gh = int_log(abs_ghval)
        except ValueError:
            continue

        if log_f <= 0:
            rho_f = 1.0
        else:
            u_f = log_f / log_B
            rho_f, table = dickman(u_f, table)

        if log_gh <= 0:
            rho_gh = 1.0
        else:
            u_gh = log_gh / log_B
            rho_gh, table = dickman(u_gh, table)
       # if first:
          #  print(f"DEBUG: u_f={u_f:.2f} u_gh={u_gh:.2f} "
          #      f"fval_bits={abs_fval.bit_length()} "
          #      f"ghval_bits={abs_ghval.bit_length()} "
         #       f"B={B} log_B={log_B:.2f}")
          #  first = False
        if rho_f > 0 and rho_gh > 0:
            contribution = scale * wi * math.exp(math.log(rho_f) + math.log(rho_gh))
          #  print(f"contribution={contribution:.6e} scale={scale:.6e} wi={wi:.6e}")
            score += contribution

    return score, table

def search_params(n,B,M,y_range,offset_range,b_range,n_gauss=32,verbose=True):
    table=get_dickman_table(20)
    best=None
    best_params=None
    total=len(y_range)*len(offset_range)*len(b_range)
    done=0
    for y in y_range:
        for offset in offset_range:
            for b in b_range:
                if y == offset:
                    done+=1
                    continue
                try:
                    score, table = e_score(n, y, offset, b, B, M, n_gauss=n_gauss, table=table)
                except Exception as e:
                    print(f"Error at y={y} offset={offset} b={b}: {e}")
                    done += 1
                    continue
                done += 1
                if best is None or score > best:
                    best = score
                    best_params = (y, offset, b)
                    if verbose:
                        print(f"  [{done}/{total}] New best: y={y} offset={offset} "
                              f"b={b}  E={score:.6e}")
    return best_params, best

def optimal_B(n):
    ln_n = int(n).bit_length() * 0.6931
    ln_ln_n = math.log(ln_n)
    return int(math.exp(math.sqrt(0.5 * ln_n * ln_ln_n)))

def get_sieve_region2(n, B):
    B=optimal_B(n) ##to do: ??? maybe remove it eventually
    y0 = math.isqrt(n)
    y_range = list(range(y0, y0 + 10))
    b_range = list(range(1, 8))
    offset_range = [int(x) for x in np.logspace(0, math.log10(float(y0)), 50)]
    best_params, best_score = search_params(n, B, lin_sieve_size, y_range, offset_range, b_range, n_gauss=16)
    baseline_score, _ = e_score(n, y0, 1, 1, B, lin_sieve_size, n_gauss=16)
    print(f"Improvement factor: {best_score/baseline_score:.3f}x")
    return best_params, best_score

def binomial_coeffs_fast(y, n):
    coeffs=[1]                
    c=1
    yk=1
    for k in range(1, n):
        c=c*(n-k+1)//k 
        yk*=y                
        coeffs.append(c*yk)
    return coeffs

def new_coeffs(f, x):
    b = 1
    tmp = [i for i in f]
    for i in range(len(f)-1):
        tmp[i] *= b
        b *= x
    tmp[-1] *= b
    return tmp

cdef construct_interval(list ret_array,partials,n,primeslist,large_prime_bound,primeslist2,small_primeslist):


    grays = get_gray_code(20)
    primelist_f=copy.copy(primeslist)


    primelist_f.insert(0,len(primelist_f)+1)
    primelist_f=array.array('q',primelist_f)
    cdef Py_ssize_t size
    primelist=copy.copy(primeslist)
    primelist.insert(0,2) ##To do: remove when we fix lifting for powers of 2
    primelist.insert(0,-1)

    primeslist_a=copy.copy(primeslist)
    primeslist_a=array.array('q',primeslist_a)
 

    smooths=[]
    coefficients=[]
    factors=[]
    jsymbols=[]
    testl=[]
    disc1_squared_list=[]

    x1=1
    y0=1
    o=1

    close_range=10
    too_close=5
    LOWER_BOUND_SIQS=1
    UPPER_BOUND_SIQS=40000
    tnum=int(((n)**0.25) /1)#(lin_sieve_size))
    
    threshold = keysize-thresvar#to do: fix this# int(math.log2((lin_sieve_size)*math.sqrt(abs(n))) - thresvar)
    
    
    LARGE_PRIME_CONST=10000
    BLOCK_SIZE=8
    NB_POLY_COARSE_EVAL=100
    NB_POLY_PRECISE_EVAL=10
    PRIME_BOUND=300
    NB_ROOTS=2
    MULTIPLIER=1
    d=degree
   # f_x,m0,m1,tmp,_ = poly_search(n, primeslist, NB_ROOTS, PRIME_BOUND, MULTIPLIER,int(pow(n, 1/(d+1))), d, NB_POLY_COARSE_EVAL,NB_POLY_PRECISE_EVAL)
    best_params,best_score=get_sieve_region2(n,primeslist[-1])
    print("best_params: "+str(best_params)+" best_score: "+str(best_score))
    
    
    y_start=best_params[0]#round(n**0.70)#132
    y_start=y_start//2
    offset_start=best_params[1]
    y_ind=0
    while y_ind < 1000:
        y=y_start+y_ind
        co=binomial_coeffs_fast(y, d)
   # print("f_x: "+str(f_x)+" m0: "+str(m0)+" m1: "+str(m1)+" g_x: "+str(g_x))
        seen=[]

        while 1:

            
            new_mod,cfact,indexes=generate_modulus(n,primeslist,seen,tnum,close_range,too_close,LOWER_BOUND_SIQS,UPPER_BOUND_SIQS,bitlen(tnum),1)
         #   new_mod=3
         #   cfact=[3]
         #   indexes=[0]
            print("new_mod: "+str(new_mod))
            if new_mod ==0:
                break

            offset_ind=0
            while offset_ind < 500:
                offset=offset_start+offset_ind
                b=1
                while b < 5:
                    partials=[]
                    fail=0
                    
                    f_x_temp=co+[-(n-(y**d-(y-offset)**d))]
                    f_x=new_coeffs(f_x_temp,b)
                    m1=f_x[0]
                    m0=-((y*2)-offset)*b
                    g_x=[m1,-m0]
                    i=0
                    while i < len(cfact):
            #    print("cfact: "+str(cfact[i]))
                        p=cfact[i]
                        f_x_c=copy.deepcopy(f_x)##To do: Fix this later...
                     #   print("f_x_c: "+str(f_x_c))
                        root = find_roots_poly(f_x_c, p)
                        if len(root)>0:
                            partials.append(p)
                            partials.append(root)
                        else:
                            fail=1
                            break
                        i+=1
         #   print("offset: "+str(offset)+" fail: "+str(fail))
                    if fail==1:
                        b+=1
                        continue
          #  print("partials: "+str(partials))
                    partials=get_partials(new_mod,partials)
          #  print("partials: "+str(partials))
                    enum=[]
                    i=0
                    while i < len(partials):
                        enum.append(partials[i+1])

                        i+=2
                    for combo2 in itertools.product(*enum):
                        tot=0
                        for r in combo2:
                            tot+=r
                        tot%=new_mod
                
                        interval=np.zeros(lin_sieve_size,dtype=np.int16)
                        t=0
                        while t < len(primeslist):  ###To do: Can move up the root calculations..
                            p=primeslist[t]
                            if new_mod%p ==0:
                                t+=1
                                continue
                            root = find_roots_poly(f_x_c, p)
                            for r in root:
                                s=solve_lin_con(new_mod,r-tot,p)
                                fval=evaluate(f_x,tot+s*new_mod)
                                if fval%(new_mod*p) != 0:
                                    print("what the fuck")

                                interval[s::p]+=round(math.log2(p))    
                            t+=1

                        np.putmask(interval, interval < 1, 0)
                        indexlist=np.nonzero(interval)

                        indexlist_x=indexlist[0]
                        ind=0
                        length=len(indexlist_x)
            
                        while ind < length:# length: 
                            x_ind=int(indexlist_x[ind])

                            x=tot+x_ind*new_mod
                            if math.gcd(x, b) != 1:
                                ind+=1
                                continue
                            fval=evaluate(f_x,x)
                            if fval%new_mod !=0:
                                print("ergh fatal")
                                sys.exit()
                            gval=evaluate(g_x,x)

                            f_x2=[1,0,n*fval*b**d*f_x[0]]
                            fval2=evaluate(f_x2,fval*f_x[0])
                 #   print("fval2: "+str(fval2))
                            local_factors, value,seen_primes,seen_primes_indexes = factorise_fast(x+offset*b,primelist_f)
                            local_factors2, value2,seen_primes2,seen_primes_indexes2 = factorise_fast(fval,primelist_f)
                            local_factors3, value3,seen_primes3,seen_primes_indexes3 = factorise_fast(gval,primelist_f)
                           # print("n: "+str(n)+" f_x: "+str(f_x)+" x+offset: "+str(bitlen(x+offset))+" fval//new_mod: "+str(bitlen(fval//new_mod))+" gval: "+str(bitlen(gval))+" modulus: "+str(new_mod))
                           
                            ##NOTE: For degree 2 when we divide fval by (x+offset*b)*fval*gval we get 1, but for degree 4, another part is added.. got to figure out the math for it
                            if fval2%(x+offset*b)*fval*gval != 0:
                                print("missing something..."+str(offset)+" f_x: "+str(f_x)+" g_x: "+str(g_x)+" fval: "+str(fval)+" gval: "+str(gval)+" b: "+str(b)+" fval2: "+str(fval2/fval)+" x: "+str(x))
                                sys.exit()
                            if value == 1 and value2 == 1 and value3 == 1:
                                if fval*f_x[0] not in coefficients and fval2 not in smooths:
                                    local_factors4, value4,seen_primes4,seen_primes_indexes4 = factorise_fast(fval2,primelist_f)
                                    if value4 ==1:
                                       # print("catastrophic fcking failure, quit math")
                                        smooths.append(fval2)
                                        factors.append(local_factors4)
                                        coefficients.append(fval*f_x[0])
                                        print("Smooths #: "+str(len(smooths))+"/"+str(base+2)+" b: "+str(b)+" offset: "+str(offset)+" f_x: "+str(f_x)+" g_x: "+str(g_x))
                                        if len(smooths)>(base+2):
                                            f1,f2=QS(n,primelist,smooths,coefficients,factors)
                                            if f1 !=0:
                                                sys.exit()            
                            ind+=1
                    b+=1
                offset_ind+=1
            y_ind+=1
    sys.exit()

   
    return      



def find_comb(n,primeslist1,primeslist2,small_primeslist):

    ret_array=[[],[],[],[]]

    

    partials={}


    large_prime_bound = primeslist1[-1] ** lp_multiplier
    construct_interval(ret_array,partials,n,primeslist1,large_prime_bound,primeslist2,small_primeslist)

    return 0

def get_primes(start,stop):
    return list(sympy.sieve.primerange(start,stop))

@cython.profile(False)
def main(l_keysize,l_workers,l_debug,l_base,l_key,l_lin_sieve_size,l_quad_sieve_size,sbase,l_degree):
    global key,keysize,workers,g_debug,base,key,quad_sieve_size,small_base,lin_sieve_size,degree
    key,keysize,workers,g_debug,base,quad_sieve_size,small_base,lin_sieve_size,degree=l_key,l_keysize,l_workers,l_debug,l_base,l_quad_sieve_size,sbase,l_lin_sieve_size,l_degree
    start = default_timer() 
    if degree%2 !=0:
        print("Fatal, degree must be even")
        sys.exit()
    if g_p !=0 and g_q !=0 and g_enable_custom_factors == 1:
        p=g_p
        q=g_q
        key=p*q
    if key == 0:
        print("\n[*]Generating rsa key with a modulus of +/- size "+str(keysize)+" bits")
        publicKey, privateKey,phi,p,q = generateKey(keysize//2)
        n=p*q
        key=n
    else:
        print("[*]Attempting to break modulus: "+str(key))
        n=key

    sys.set_int_max_str_digits(1000000)
    sys.setrecursionlimit(1000000)
    bits=bitlen(n)
    primeslist=[]
    primeslist1=[]
    primeslist2=[]
    small_primeslist=[]
    mod_primeslist=[]
    print("[i]Modulus length: ",bitlen(n))
    count = 0
    num=n
    while num !=0:
        num//=10
        count+=1
    print("[i]Number of digits: ",count)
    print("[i]Gathering prime numbers..")
    primeslist.extend(get_primes(3,20000000))
    i=0
    while len(primeslist1) < base:
        if n%primeslist[i]==0:
            i+=1
            continue
        if len(small_primeslist)<small_base:
            small_primeslist.append(primeslist[i])
    
        primeslist1.append(primeslist[i])

        i+=1
    
    while len(primeslist2) < qbase:
        if n%primeslist[i]==0:
            i+=1
            continue
        primeslist2.append(primeslist[i])
        i+=1
    launch(n,primeslist1,primeslist2,small_primeslist)     
    duration = default_timer() - start
    print("\nFactorization in total took: "+str(duration))