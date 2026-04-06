#!python
#cython: language_level=3
# cython: profile=False
# cython: overflowcheck=False
###Author: Essbee Vanhoutte
###WORK IN PROGRESS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
###Improved QS Variant 

##References: 
##-I have borrowed many of the SIQS related optimizations from here: https://stackoverflow.com/questions/79330304/optimizing-sieving-code-in-the-self-initializing-quadratic-sieve-for-pypy
####NFS related code is borrowed from: https://github.com/basilegithub/General-number-field-sieve-Python (note: Very impressively written, helped me big time, thanks)


###To build: python3 setup.py build_ext --inplace

from sympy import symbols, Poly, Matrix,discriminant
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

key=0                 #Define a custom modulus to factor
build_workers=8
keysize=150           #Generate a random modulus of specified bit length
workers=1 #max amount of parallel processes to use
quad_co_per_worker=1 #Amount of quadratic coefficients to check. Keep as small as possible.
base=1_000
small_base=1000
qbase=1
quad_sieve_size=10
g_debug=0 #0 = No debug, 1 = Debug, 2 = A lot of debug
g_lift_lim=0.5
thresvar=10  ##Log value base 2 for when to check smooths with trial factorization. Eventually when we fix all the bugs we should be able to furhter lower this.
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

def power2(poly, f, p, exp):
    if exp == 1: return poly

    tmp = power2(poly, f, p, exp>>1)
    tmp = poly_prod(tmp, tmp)

    if exp&1:
        tmp = poly_prod(tmp, poly)
        return div_poly_mod(tmp, f, p)
    
    else: return div_poly_mod(tmp, f, p)

def power(x, y, p):
    res = 1;
    x = x % p;
    while (y > 0):
        if (y & 1):
            res = (res * x) % p;
        y = y>>1; # y = y/2
        x = (x * x) % p;
    return res;

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
        
def QS(n,factor_list,sm,flist,x_list,x_f_list):#,jsymbols,testl,primeslist2,disc1_squared_list):#,disc_sr_list,pval_list,pflist):
    g_max_smooths=base+2#+qbase
    if len(sm) > g_max_smooths*10000000: 
        del sm[g_max_smooths:]
       # del xlist[g_max_smooths:]
        del flist[g_max_smooths:]  
    M2 = build_matrix(factor_list, sm, flist,x_f_list)#,pflist)
    null_space=solve_bits(M2,factor_list,len(sm))
    f1,f2=extract_factors(n, sm, null_space,x_list)#,disc_sr_list,pval_list,pflist)
    if f1 != 0:
        print("[SUCCESS]Factors are: "+str(f1)+" and "+str(f2))
        return f1,f2   
   # print("[FAILURE]No factors found")
    return 0,0

def extract_factors(N, relations, null_space,x_list):#,disc_sr_list,pval_list,pflist):
    n = len(relations)
    for vector in null_space:
        prod_left = 1
        prod_right = 1
        pval=1
        disc_sr=1
        xy=1
        x=1
        for idx in range(len(relations)):
            bit = vector & 1
            vector = vector >> 1
            if bit == 1:
                prod_left *= relations[idx]
                x*=x_list[idx]
             #   print("polyval:  "+str(relations[idx])+" disc constant "+str(x_list[idx]))
            idx += 1
        prod_right=x
        sqrt_right = math.isqrt(x)
        sqrt_left = math.isqrt(prod_left)#prod_left
        if sqrt_left**2 != prod_left:
            print("horrible error")
            sys.exit()
       # print(" polyval sqrt: "+str(sqrt_left)+" disc constant sqrt: "+str(sqrt_right))#+" zx*zxy: "+str(x))
        ###Debug shit, remove for final version
        sqr1=prod_left%N 
        sqr2=prod_right%N
        if sqrt_right**2 != prod_right:
            print("not a square in the integers")
            sys.exit()
         #   time.sleep(10000)

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

def build_matrix(factor_base, smooth_nums, factors,x_f_list):#,pflist):
    fb_map = {val: i for i, val in enumerate(factor_base)}

    ind=1

    M2=[0]*((base+2)*2)#+qbase)#+2+qbase)
    for i in range(len(smooth_nums)):
        for fac in factors[i]:
            idx = fb_map[fac]
            M2[idx] |= ind
        ind = ind + ind       

  #  offset=(base+2)-1
  #  ind=1
    
 #   for i in range(len(xy_f_list)):
 #       j=0
 #       while j < len(xy_f_list[i]):
 #           fac=xy_f_list[i][j]
 #           if fac == -1 or fac ==0:
 #               M2[j+offset] |= ind
 #           j+=1
 #       ind = ind + ind

    offset=(base+2)
    ind=1
    for i in range(len(x_f_list)):
        for fac in x_f_list[i]:
            idx = fb_map[fac]
            M2[idx+offset] |= ind
        ind = ind + ind
    return M2

@cython.profile(False)
def launch(n,primeslist1,primeslist2,small_primeslist):
    manager=multiprocessing.Manager()
    return_dict=manager.dict()
    jobs=[]
    start= default_timer()
    print("[i]Creating coefficient datastructure... this can take a while...")
    primeslist1c=copy.deepcopy(primeslist1)
    plists=[]

    hmap,hmap2=create_hashmap(n,primeslist1)
    duration = default_timer() - start
    print("[i]Creating coefficient datastructure in total took: "+str(duration))

    print("[*]Launching attack with "+str(workers)+" workers\n")
    find_comb(n,hmap,hmap2,primeslist1,primeslist2,small_primeslist)

    return 


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

def solve_congruence(a, b, p):
    if a == 0:
        if b == 0:
            return [] 
        else:
            return None 
    else:
        return b * modinv(a, p) % p



def solve_quadratic_congruence(a, b, c, p):
    if a == 0:
        return [solve_congruence(b, -c, p)]
    else:
        a_inv = modinv(a, p)
        ba = (b * a_inv) % p
        ca = (c * a_inv) % p
        b_div_2 = (ba * modinv(2, p)) % p
        alpha = (pow(b_div_2, 2, p) - ca) % p
        if jacobi(alpha,p)==-1:
            return -1
        if jacobi(alpha,p)==0:
            ##To do: Sort out the number theory for this
            i=0
            while i < p:
                if (a*i**2+b*i+c)%p ==0:
                    ret=i 
                    break
                i+=1

            return [ret]
        y = tonelli(alpha, p)
        if y is None:
            return -1
        x1 = (y - b_div_2) % p
        x2 = (p - y - b_div_2) % p
        return [x1,x2]

def lift_root(z,y,n, root, p, exp):
    r =(z*root**2+root*y-n)%p**(exp+1)
    deriv = 2*z*root+y
    inv = modinv(deriv, p)
    t = ((-r//p**exp)*inv)%p
    new_root = (root+t*p**exp)%p**(exp+1)
    return new_root

def roots(g, p):
    if len(g) == 1: return []
    if len(g) == 2: return [-g[1]*modinv(g[0], p)%p]
    if len(g) == 3:
        tmp = (g[1]*g[1]-4*g[0]*g[2])%p
        if tmp == 0: return [-g[1]*modinv(2*g[0], p)%p]
        if compute_legendre_character(tmp, p) == -1: return []
        tmp = compute_sqrt_mod_p(tmp, p)*modinv(2*g[0], p)%p
        return [(-g[1]*modinv(g[0]<<1, p)+tmp)%p, (-g[1]*modinv(g[0]<<1, p)-tmp)%p]
    
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

def poly_prod(a, b):
    res = [0]*(max(len(a), len(b))+min(len(a), len(b))-1)

    for i in range(len(a)):
        for j in range(len(b)):
            res[i+j] += a[i]*b[j]

    return res

def div_poly_mod(a, tmp_b, p):
    remainder = [i%p for i in a]
    b = [i%p for i in tmp_b]
    
    #print(remainder, b)
    while not b[0]: del b[0]
    
    difference = len(a)-len(b)+1
    coeff = modinv(-b[0], p)
    for j in range(difference):
        if remainder[j]:
            quotient = remainder[j]*coeff%p
            remainder[j] = 0
            for k in range(1,len(b)): remainder[j+k] = (remainder[j+k]+quotient*b[k]%p)%p
            
    for k in range(len(remainder)):
        if remainder[k]: return remainder[k:]
        
    return [0]

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
        

def gcd_mod(f, poly, p):
    while poly != [0]*len(poly):
        (f, poly) = (poly, div_poly_mod(f, poly, p))
    return f

def quotient_poly_mod(a, b, p):
    remainder = [i%p for i in a]
    b = [i%p for i in b]
    
    while not b[0]: del b[0]
    
    difference = len(a)-len(b)+1
    coeff = modinv(-b[0], p)
    res = [0]*difference

    for j in range(difference):
        quotient = remainder[j]*coeff%p
        res[j] = -quotient
        for k in range(len(b)):
            remainder[j+k] = (remainder[j+k]+quotient*b[k]%p)%p
            
    for k in range(len(res)):
        if res[k]: return res[k:]
        
    return [0]

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

def solve_roots(prime,n):
    hmap_p=[]
    hmap_p2=[]
    k=0
    while k < prime and k < quad_sieve_size+1:
        hmap_p2.append([])
        coeff=[(-n*k)%prime]


        d=degree
        d_ind=0
        while d_ind < d:
            coeff.insert(0,0)

            d_ind+=1
        coeff[0]+=1
        ranges = [range(start, prime) for start in coeff[:-1]]
        for combo in itertools.product(*ranges):
            current = list(combo) +  [coeff[-1]]
            if degree > 2:
                root = find_roots_poly(current, prime)
                if len(root)> 0:
                    hmap_p2[-1].append(list(combo))
            else:
                test=current[1]**2+4*n*k*current[0]
                if jacobi(test%prime,prime) != -1:
                    hmap_p2[-1].append(list(combo))
            #print("root: "+str(root)+" current: "+str(current)+" prime: "+str(prime))
        k+=1

    return hmap_p2

def create_hashmap(n,primeslist):
    i=0
    hmap=[]
    hmap2=[]
    while i < len(primeslist):
        hmap_p=solve_roots(primeslist[i],n)
        hmap.append(hmap_p)
        hmap2.append([])
        i+=1

    return hmap,hmap2


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

@cython.boundscheck(False)
@cython.wraparound(False)
cdef factorise_fast_quads(value,long long [::1] factor_base):
    factors = set()
    if value % 2 == 0:
        factors ^= {2}
        value //= 2
        if value % 2 == 0:
            return -1, -1
    length=factor_base[0]#len(factor_base)#factor_base[0]
    cdef Py_ssize_t i=1
    while i < length:
        factor=factor_base[i]
        if value % factor == 0:
            factors ^= {factor}
            value //= factor
            if value % factor == 0:
                return -1, -1
        i+=1
    return factors, value

def filter(qbase,n,start,size):
    ###Note: We look for quadratic coefficients that factor over the factor base but have no even exponents. This garantuees unique results
    interval= np.zeros(size+1,dtype=np.int32)
    interval2= np.zeros(size+1,dtype=np.int32)

    length=qbase[0]
    i=1
    while i < length:
        prime=qbase[i]
        log=round(math.log2(prime))
        dist=start%prime
        dist=-dist%prime
        interval[dist::prime]+=log
        if (start+dist)%prime !=0:
            print("fatal error")
        i+=1
    
    
    valid_quads=[]
    valid_quads_factors=[]

    roots=[]
    i=0
    while i < size:
        if start+i > 999:
            threshold=round(math.log2((start+i)*0.50))
        if start+i < 1000 or interval[i]>threshold:
            quad_local_factors, quad_value = factorise_fast_quads(start+i,qbase) 
         #   print("checing")
            if quad_value != 1:
                i+=1
                continue

            interval2[i]=1
         #   print("found")
            valid_quads.append(start+i)
            valid_quads_factors.append(quad_local_factors)
        
        i+=1
    return valid_quads,valid_quads_factors,interval2

def solve_quad_integers(a,b,c):
    ##To do: Could I use this for the sieving process?
    disc=b**2+4*a*c

    if disc < 0:
        return -1
    if disc == 0:
        return -1
    test=math.isqrt(disc)
    if test**2!=disc:
      #  print("fail")
        return -1

    return [(-b+test)//(2*a),(-b-test)//(2*a)] 
#@cython.boundscheck(False)
#@cython.wraparound(False)
def generate_modulus(n,primeslist,seen,tnum,close_range,too_close,LOWER_BOUND_SIQS,UPPER_BOUND_SIQS,tnum_bit,hmap,quad):
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
            if  jacobi((-quad*n)%primeslist[randindex],primeslist[randindex])!=1:#hmap[randindex][1]!=quad%primeslist[randindex]:
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

def find_xy(primeslist,n,x,collected,z,lin):
    xy=x+lin

    i=0
    while i < len(collected):
        
        prime=collected[i][0]
        i+=1
    return xy

def retrieve(hmap,primeslist,x,lin,n):

    collected=[1]*(quad_sieve_size+1)
    mult_list=[]
    interval=array.array("i",[0]*lin_sieve_size)
    i=0
    while i < quad_sieve_size:
        mult_list.append(array.array("i",[0]*lin_sieve_size))
        i+=1
    i=0
    while i < len(hmap):
        
        prime=primeslist[i]#**2
        log=round(math.log2(prime))
        c=1
        while c <  len(hmap[i])+1:
            z=c
            xy=(z*x)+lin
            if (z*x)%prime !=0 and xy%prime !=0:
                mul=solve_lin_con(xy,hmap[i][c-1][(x)%prime],prime)
                t=mul   
                y=(xy*mul)-x*z
                if (z*x**2+y*x-n)%prime !=0:
                    print("fatal error: "+str(z)+" p: "+str(prime)+" x: "+str(x))
                    sys.exit()
                while t < lin_sieve_size:
                    mult_list[c-1][t]+=log
                    t+=prime
       
              #  mult_list[c].append([prime,mul])

             #   collected.extend([prime,c])
            c+=1
        i+=1

    return mult_list

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

def coefficients(P):
    Pc = []
    for i in range(P.degree(), -1, -1):
        Pc.append(P.nth(i))
    return tuple(Pc)
def sylvester(P, Q):
    rows = []
    m = P.degree()
    n = Q.degree()
    size = m + n
    CP = coefficients(P)
    CQ = coefficients(Q)
    for i in range(size):
        tail = []
        row = []
        if i in range(0, n):
            row = list(CP)
            row.extend((n-1-i)*[0])
            row[:0] = [0]*i
            rows.append(row)
        if i in range(n, size):
            row = list(CQ)
            row.extend((size-1-i)*[0])
            row[:0] = [0]*(size-len(row))
            rows.append(row)
    return Matrix(rows)
def resultant(P, Q):
    return sylvester(P,Q).det()

def evaluate(f, x):
    res = 0

    for i in range(len(f)-1):
        res += f[i]
        res *= x

    res += f[-1]

    return res
def new_coeffs(f, x):
    b = 1
    tmp = [i for i in f]
    for i in range(len(f)-1):
        tmp[i] *= b
        b *= x
    tmp[-1] *= b
    return tmp
    
def sieve(length, f_x, rational_base, algebraic_base, m0, m1, b, leading):
    ##Code borrowed from: https://github.com/basilegithub/General-number-field-sieve-Python
 
    pairs = []
    tmp_poly = new_coeffs(f_x, b)
  #  print("b: ",b)
    sieve_array = [0]*(length<<1)
    for q, p in enumerate(rational_base):
        tmp_len = length%p
        log = round(math.log2(p))#logs[q]

        if m1%p:
            root = (tmp_len+b*m0*modinv(m1, p))%p
            for i in range(root, len(sieve_array), p): sieve_array[i] += log

        for r in algebraic_base[q]:
            root = (tmp_len+b*r)%p
            for i in range(root, len(sieve_array), p): sieve_array[i] += log

    if b&1:
     #   print("hit")
        eval2 = -length*m1-b*m0
        a = -length


        tmp = [0]*len(tmp_poly)
        for j in range(len(tmp_poly)): 
            tmp[j] = evaluate(tmp_poly, -length+j)
        tmp_debug=copy.copy(tmp)
        for q in range(1, len(tmp_poly)):
            for k in range(len(tmp_poly)-1, q-1, -1): 
                tmp[k] -= tmp[k-1]

        eval1 = tmp[0]
        for k in range(len(sieve_array)):
            if a and math.gcd(a, b) == 1 and eval2:
                eval = abs(eval1*eval2)
               # print("eval: ",eval)
                if eval != 0 and sieve_array[k]:# > eval.bit_length()-offset:
                    pairs.append([[-b, a*leading], eval1, [[1 ,1]], eval2, [1], [-b,a], 1])
            a += 1
            eval2 += m1
            eval1 += tmp[1]
            for q in range(1, len(tmp_poly)-1): 
                tmp[q] += tmp[q+1]
    else:
     #   print("hit2")
        init = 0
        if not length&1:
            length -= 1
            init = 1
        eval2 = -length*m1-b*m0
        a = -length
        tmp = [0]*len(tmp_poly)
        for j in range(len(tmp_poly)): 
            tmp[j] = evaluate(tmp_poly, -length+2*j)
        for q in range(1, len(tmp_poly)):
            for k in range(len(tmp_poly)-1, q-1, -1): 
                tmp[k] -= tmp[k-1]

        eval1 = tmp[0]
        for k in range(init, len(sieve_array), 2):
            if math.gcd(a, b) == 1 and eval2:
                eval = abs(eval1*eval2)
                if eval != 0 and sieve_array[k]:# > eval.bit_length()-offset:
                  #  if b == 366:
                     #   print("ADDING eval2: "+str(eval2)+" a: "+str(a)+" m1: "+str(m1)+" b: "+str(b)+" m0: "+str(m0))

                    pairs.append([[-b, a*leading], eval1, [[1, 1]], eval2, [1], [-b,a], 1])
            a += 2
            eval2 += m1<<1
            eval1 += tmp[1]
            for q in range(1, len(tmp_poly)-1): 
                tmp[q] += tmp[q+1]
    return pairs

def trial(pair, primes, div):
    large1 = 1
    large2 = 1
    
    result = abs(pair[3])
    for p in primes:
        while not result%p: 
            result //= p
        if result == 1: 
            break

    if result > 1:
            return False, 1, 1, 1, 1
    else: 
       # print("wohooo1")
        large1, large2 = 1, 1
   # print("hit")
    result = abs(pair[1])//div
   # print("result: ",result)
    for p in primes:
        while not result%p: 
            result //= p
        if result == 1: 
          #  print("WOHOOO!!")
            return True, 1, 1, large1, large2

    if result > 1:
        return False, 1, 1, 1, 1
        
    return True, 1, 1, 1, 1
    

def evaluate(f, x):
    res = 0

    for i in range(len(f)-1):
        res += f[i]
        res *= x

    res += f[-1]

    return res

def evaluate_x2(f, x):
    res = 0
    i=0
    while i < len(f)-2:
        res += f[i]
        res *= x
        i+=1
    res+=f[i]

    return res



cdef construct_interval(list ret_array,partials,n,primeslist,hmap,hmap2,large_prime_bound,primeslist2,small_primeslist):
    


    print("[i]Entering attack loop")
    k=1
    while k < quad_sieve_size+1:
        print("[i]Building interval")
        interval=np.zeros([lin_sieve_size,lin_sieve_size],dtype=np.int16)
        t=0
        while t < len(hmap):
            prime=primeslist[t]
            coeff=[(-n*k)%prime]


            d=degree
            d_ind=0
            while d_ind < d:
                coeff.insert(0,0)

                d_ind+=1
            coeff[0]+=1
            ranges = [range(start, prime) for start in coeff[:-1]]
            for combo in itertools.product(*ranges):
                cur=list(combo)
               # if cur not in hmap[t][k]:
                test=cur[1]**2+4*n*k*cur[0]
                if jacobi(test,prime)==-1:
                    if cur[0] < lin_sieve_size and cur[1] < lin_sieve_size:
                        interval[cur[0]::prime,cur[1]::prime]=1
            t+=1
       # print("interval: ",interval)
        print("[i]Checking interval")
        coeff=[-n*k]
        d=degree
        d_ind=0
        while d_ind < d:
            coeff.insert(0,0)
            d_ind+=1
        coeff[0]+=1
        ranges = [range(start, lin_sieve_size) for start in coeff[:-1]]
        for combo in itertools.product(*ranges):
            current = list(combo) +  [coeff[-1]]
    
            if interval[current[0],current[1]]==1:
                continue
            ##This wont work for higher degrees. For higher degrees I suppose we'll need to calculate roots mod p and use CRT.
            test=current[1]**2+4*n*k*current[0]

            print("found one: "+str(current)+" my_for: "+str(test))
            test_sqr=math.isqrt(test)
            if test_sqr**2 == test:
                gcdtest=gcd(test_sqr+current[1],n)
                if gcdtest != 1 and gcdtest != n:
                    print("factors of "+str(n)+" are: "+str(gcdtest)+" and: "+str(n//gcdtest))
                    sys.exit()
                else:
                    print("[i]Might need to increase factor base size")
            else:
                print("[i]Might need to increase factor base size")
        k+=1
    return 




def find_comb(n,hmap,hmap2,primeslist1,primeslist2,small_primeslist):

    ret_array=[[],[],[],[]]

    

    partials={}


    large_prime_bound = primeslist1[-1] ** lp_multiplier
    construct_interval(ret_array,partials,n,primeslist1,hmap,hmap2,large_prime_bound,primeslist2,small_primeslist)

    return 0

def get_primes(start,stop):
    return list(sympy.sieve.primerange(start,stop))

@cython.profile(False)
def main(l_keysize,l_workers,l_debug,l_base,l_key,l_lin_sieve_size,l_quad_sieve_size,l_d):
    global key,keysize,workers,g_debug,base,key,quad_sieve_size,small_base,lin_sieve_size,qbase,degree
    key,keysize,workers,g_debug,base,quad_sieve_size,lin_sieve_size,degree=l_key,l_keysize,l_workers,l_debug,l_base,l_quad_sieve_size,l_lin_sieve_size,l_d
    start = default_timer() 

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
        keysize=bitlen(n)
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
   # qbase=base
    while num !=0:
        num//=10
        count+=1
    print("[i]Number of digits: ",count)
    print("[i]Gathering prime numbers..")
    print("[i]!!!!!!!! Only -d 2 (degree) polynomials supported for now.. work in progresss")
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