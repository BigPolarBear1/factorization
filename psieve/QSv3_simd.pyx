#!python
#cython: language_level=3
#cython: profile=False
#cython: overflowcheck=False
###Author: Essbee Vanhoutte
###WORK IN PROGRESS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
###Cuda_QS_Variant
###This code is pro-LGBTQ.

##References: I have borrowed many of the optimizations from here: https://stackoverflow.com/questions/79330304/optimizing-sieving-code-in-the-self-initializing-quadratic-sieve-for-pypy


import pstats, cProfile
import sympy
import random
import itertools
import sys
import argparse
import time
import copy
from timeit import default_timer
import math
import gc
import array
import numpy as np
import cython
cimport cython
import os
from sympy.ntheory.residue_ntheory import jacobi_symbol
from sympy import kronecker_symbol


min_lin_sieve_size=10_000
max_bound=10_000_000
key=0                 #Define a custom modulus to factor
build_workers=8
keysize=150           #Generate a random modulus of specified bit length
workers=1 #max amount of parallel processes to use
quad_co_per_worker=1 #Amount of quadratic coefficients to check. Keep as small as possible.
base=1_000
qbase=10
lin_sieve_size=1
lin_sieve_size2=10_000_000
quad_sieve_size=10
g_debug=0 #0 = No debug, 1 = Debug, 2 = A lot of debug
g_lift_lim=0.5
thresvar=30  ##Log value base 2 for when to check smooths with trial factorization. Eventually when we fix all the bugs we should be able to furhter lower this.
thresvar2=30
dupe_max_prime=1_000_000
lp_multiplier=2
min_prime=1
g_max_diff_similar=5
g_enable_custom_factors=0
g_p=107
g_q=41
mod_mul=0.5
g_max_exp=20
quad_per_interval=1
lift_lim=1
k_max = 12#_000
b_max=10_000
quad_sign="neg" 
##Key gen function##
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
def modinv(n,p):
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

def gcd(a,b): # Euclid's algorithm ##
    if b == 0:
        return a
    elif a >= b:
        return gcd(b,a % b)
    else:
        return gcd(b,a)

def formal_deriv(y,x,z):
    result=(z*2*x)+(y)
    return result


        
def QS(n,factor_list,sm,flist,x_list,factor_list2):#,jsymbols,testl,primeslist2,disc1_squared_list):#,disc_sr_list,pval_list,pflist):
    g_max_smooths=base+2#+qbase
    if len(sm) > g_max_smooths*10000000: 
        del sm[g_max_smooths:]
       # del xlist[g_max_smooths:]
        del flist[g_max_smooths:]  
    M2 = build_matrix(factor_list, sm, flist,factor_list2)#,pflist)
    null_space=solve_bits(M2,factor_list,len(sm))
    f1,f2=extract_factors(n, sm, null_space,x_list,flist)#,disc_sr_list,pval_list,pflist)
    if f1 != 0:
        print("[SUCCESS]Factors are: "+str(f1)+" and "+str(f2))
        sys.exit()
        return f1,f2   
   # print("[FAILURE]No factors found")
    return 0,0

def extract_factors(N, relations, null_space,x_list,factor_list):#,disc_sr_list,pval_list,pflist):
    n = len(relations)
    for vector in null_space:
        prod_left = 1
        prod_right = 1
        pval=1
        disc_sr=1
        xy=1
        x=1
        count=0
        for idx in range(len(relations)):
            bit = vector & 1
            vector = vector >> 1
            if bit == 1:
                count+=1
                prod_left *= relations[idx]
                prod_right *=x_list[idx]
                x*=x_list[idx]
               # print("polyval:  "+str(relations[idx])+" disc constant "+str(x_list[idx])+" factors: "+str(factor_list[idx]))
            idx += 1

        sqrt_right = math.isqrt(prod_right)
        sqrt_left = math.isqrt(prod_left)#prod_left
        if sqrt_left**2 != prod_left:
            print("horrible error")
            sys.exit()
       # print(" polyval sqrt: "+str(sqrt_left%N)+" disc constant sqrt: "+str(sqrt_right%N))#+" zx*zxy: "+str(x))
        ###Debug shit, remove for final version
        sqr1=prod_left%N 
        sqr2=prod_right%N
        if sqrt_right**2 != prod_right:
            print("not a square in the integers")
            sys.exit()
         #   time.sleep(10000)

        if sqr1 != sqr2:
            print("ERROR ERROR")
            #time.sleep(10000)
        ###End debug shit#########
        sqrt_left = sqrt_left % N
        sqrt_right = sqrt_right % N
        factor_candidate = gcd(N, abs(sqrt_right+sqrt_left))


        if factor_candidate not in (1, N):
          #  print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!: "+str(factor_candidate))#+" sm: "+str(sqrt_right)+" root: "+str(sqrt_left))
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

def build_matrix(factor_base, smooth_nums, factors,factor_list2):#,pflist):
    fb_map = {val: i for i, val in enumerate(factor_base)}

    ind=1

    M2=[0]*((len(factor_base)+2)*2)#+qbase)#+2+qbase)
    for i in range(len(smooth_nums)):
        for fac in factors[i]:
            idx = fb_map[fac]
            M2[idx] |= ind
        ind = ind + ind       

    offset=(len(factor_base)+2)-1
    ind=1
    for i in range(len(factor_list2)):
        for fac in factor_list2[i]:
            idx = fb_map[fac]
            M2[idx+offset] |= ind
        ind = ind + ind
    return M2

def launch(n,primeslist,primeslist2):
    print("[*]Launching attack")# with "+str(workers)+" workers\n")
    ret_array=[[],[],[],[]]
    sbase=copy.deepcopy(primeslist[0:50])
   # print("fbase: "+str(fbase))
    print("[i]Building psieve Residue Map (to do: some duplication here from merging two algos, fix later)")

    found=0

    primelist_f=copy.copy(primeslist)
    primelist_f.insert(0,len(primelist_f)+1)
    primelist_f=array.array('q',primelist_f)

    primelist=copy.copy(primeslist)
    primelist.insert(0,2) ##To do: remove when we fix lifting for powers of 2
    primelist.insert(0,-1)

    
    qlist=copy.copy(primeslist2)
    qlist.insert(0,len(qlist)+1)
    qlist=array.array('q',qlist)

    primeslist_a=copy.copy(primeslist)
    primeslist_a=array.array('q',primeslist_a)

    div=1
    while div < 100:

        psievefound=psieve(n,ret_array,primelist_f,primeslist,div,sbase)
        if psievefound !=0:
            print("[*](Psieve)Trying linear algebra after succesful psieve run")
            test,test2=QS(n,primelist,ret_array[0],ret_array[2],ret_array[1],ret_array[3])
                    

            if test !=0:
                print("\n\n\n\nFound at: ",len(ret_array[0]))
                sys.exit()

        div+=1 
    
    return 

def squareRootExists(n,p,b,a):
    b=b%p
    c=n%p
    bdiv = (b*modinv(2,p))%p
    alpha = (pow_mod(bdiv,2,p)-c)%p
    if alpha == 0:
        return 1
    
    if jacobi(alpha,p)==1:
        return 1
    return 0

def pow_mod(base, exponent, modulus):
    return pow(base,exponent,modulus)  

def lift(exp,co,r,n,z,z2,prime):
    z_inv=modinv(z,prime**exp)
    c=(-z*n)%prime**exp
    temp_r=r*z
    zz=2*temp_r
    zz=pow(zz,-1,prime)
    x=((c-temp_r**2)//prime)%prime
    y=(x*zz)%prime
    new_r=(temp_r+y*prime)%prime**exp
    root2=(new_r*z_inv)%prime**exp

    co2=(formal_deriv(0,new_r,z))%(prime**exp) 
    ret=[]
    ret.extend([co2,new_r])
    new_eq=(z*new_r**2+n)%prime**exp
    if new_eq !=0:
        print("Something went wrong while lifting z: "+str(z)+" new_r: "+str(new_r))
        time.sleep(10000)
    return ret

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

def equation2(y,x,n,mod,z,z2):
    rem=z*(x**2)+y*x-n*z2
    rem2=rem%mod
    return rem2,rem

@cython.boundscheck(False)
@cython.wraparound(False)
cdef factorise_fast(value,long long [::1] factor_base):
    if value == 0:
        print("blah")
        return [],-1
    factors = set()
    if value < 0:
        factors ^= {-1}
        value = -value
    while value % 2 == 0:
        factors ^= {2}
        value //= 2

    length=factor_base[0]#len(factor_base)#factor_base[0]
    cdef Py_ssize_t i=1
    while i < length:
        factor=factor_base[i]
        while value % factor == 0:
            factors ^= {factor}
            value //= factor
        i+=1
    return factors, value

@cython.boundscheck(False)
@cython.wraparound(False)
cdef factorise_fast2(value,long long [::1] factor_base):
    seen_primes=[]
    if value == 0:
        print("blah")
        return [],-1,[]
    factors = set()
    if value < 0:
        factors ^= {-1}
        value = -value
    while value % 2 == 0:
        factors ^= {2}
        value //= 2

    length=factor_base[0]#len(factor_base)#factor_base[0]
    cdef Py_ssize_t i=1
    while i < length:
        factor=factor_base[i]
        while value % factor == 0:
            seen_primes.append(factor)
            factors ^= {factor}
            value //= factor
        i+=1
    return factors, value,seen_primes

def get_root(p,b,a):
    a_inv=modinv((a%p),p)
    if a_inv == None:
        return -1
    ba=(b*a_inv)%p 
    bdiv = (ba*modinv(2,p))%p
    return bdiv%p

def solve_lin_con(a,b,m):
    ##ax=b mod m
    #g=gcd(a,m)
    #a,b,m = a//g,b//g,m//g
    return pow(a,-1,m)*b%m  
    
#@cython.boundscheck(False)
#@cython.wraparound(False)
cdef factorise_fast_debug(value,long long [::1] factor_base,max_exp):
    factors = set()
    seen_primes=[]
    if value < 0:
        factors ^= {-1}
        value = -value
    while value % 2 == 0:
        factors ^= {2}
        value //= 2
    #    seen_primes.append(2)
    length=factor_base[0]#len(factor_base)#factor_base[0]
    cdef Py_ssize_t i=1
    while i < length:
        factor=factor_base[i]
        exp=1
        while value % factor == 0:
            factors ^= {factor}
            value //= factor
            if exp < (max_exp+1):
                seen_primes.append(factor)
            exp+=1
        i+=1
    return factors, value,seen_primes



def evaluate(f, x):
    res = 0

    for i in range(len(f)-1):
        res += f[i]
        res *= x

    res += f[-1]

    return res


def power2(poly, f, p, exp):
    if exp == 1: return poly

    tmp = power2(poly, f, p, exp>>1)
    tmp = poly_prod(tmp, tmp)

    if exp&1:
        tmp = poly_prod(tmp, poly)
        return div_poly_mod(tmp, f, p)
    
    else: return div_poly_mod(tmp, f, p)

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

def find_r(mod,total):
    mo,i=mod,0
    while (total%mod)==0:
        mod=mod*mo
        i+=1
    return i

def brute_force_padic_solutions2(prime,k,n,a,sqr):
    ##to do: delete later.. just for verifications
    solutions=[]
    blist=[]
    exp=1
    b=0
    while b < prime**exp:
       # print("b: "+str(b)+" prime: "+str(prime))
        poly=[1,b*sqr,-n*k]
        #poly=[1,-b,n*k]
        roots=[]
        x=0
        while x < prime**exp:
            if evaluate(poly,x)%prime**exp ==0:
                roots.append(x)

            x+=1
        if len(roots)>0:
            solutions.extend([b,roots,0])
        b+=1
   # print("solutions: "+str(solutions))
    max_lift=6
    exp+=1
    while exp < max_lift:
        am_to_lift=0
        sqr_lifted=lift_root2([1,0,-a],sqr,prime,exp)
        new_solutions=[]
        i=0
        while i < len(solutions):
            b=solutions[i]
            while b < prime**exp:
                #print("b: "+str(b)+" prime: "+str(prime)+" exp: "+str(exp)+" prime**exp: "+str(prime**exp))
                poly=[1,(b*sqr_lifted)%prime**exp,(-n*k)%prime**exp]
                der=get_derivative(poly)
                #poly=[1,-b,n*k]
                new_roots=[]
                hit=0
                if solutions[i+2]==0:
                    roots=solutions[i+1]
                    j=0
                    while j < len(roots):
                        r=roots[j]
                 #   print("r: "+str(r)+" b: "+str(b))
                        while r < prime**exp:
                            if evaluate(poly,r)%prime**exp ==0:
                                new_roots.append(r)
                                deriv=evaluate(der,r)
                                r0=find_r(prime,deriv)
                                ceiling=(prime*r0)
                                ceiling=prime**ceiling
                                if prime**exp > ceiling:
                                    hit=1
                              #  print("******************prime**exp: "+str(prime**exp)+" b: "+str(b)+" r: "+str(r)+" "+str(deriv)+" ceiling: "+str(ceiling))
                           # else:
                             #   print("prime**exp: "+str(prime**exp)+" b: "+str(b)+" r: "+str(r)+" "+str(deriv)+" ceiling: "+str(ceiling))
                            r+=prime**(exp-1)
                        j+=1
                    if len(new_roots)>0:
                    #print("b: "+str(b)+" "+str(len(new_roots)))
                        if hit==1:
                            new_solutions.extend([b,new_roots,0])
                        else:
                            new_solutions.extend([b,new_roots,0])
                            am_to_lift+=1
                
                #if solutions[i+2]==1 and len(roots)==0:
                #    print("fatal error")
                #    sys.exit()
                #if solutions[i+2]==1 and hit==0:
                #    print("fata; error2")
                #    sys.exit()



                b+=prime**(exp-1)
            i+=3
        solutions=new_solutions
       # print("solutions: "+str(prime**exp)+" "+str(len(solutions)//3))
        if am_to_lift==0:
          #  print("!!!!!!!!!!!!!!!!FULLY LIFTED AT: "+str(exp))
            break
        if exp+1 == max_lift:
            break
        exp+=1

    
    blist.append(prime**(exp))
    blist.append([])
    i=0
    while i < len(solutions):
        if solutions[i+2]==1:
            texp=solutions[i+1]
            tempb=solutions[i]
            while tempb < prime**(exp):
                blist[-1].append(tempb)
                ##To do: can add a verification loop here later if I'm bored to make sure everything added is correct

                tempb+=prime**texp
        else:    
            blist[-1].append(solutions[i])

        i+=3
    blist[-1].sort()
    #print("solutions len: "+str(len(blist[-1]))+" k: "+str(k))

    return blist

def brute_force_padic_solutions(prime,k,n,a,sqr):
    ##Add proper hensel later and use this to verify
    solutions=[]
    blist=[]
    exp=1
    b=0
    while b < prime**exp:
       # print("b: "+str(b)+" prime: "+str(prime))
        poly=[1,b*sqr,-n*k]
        #poly=[1,-b,n*k]
        roots=[]
        x=0
        while x < prime**exp:
            if evaluate(poly,x)%prime**exp ==0:
                roots.append(x)

            x+=1
        if len(roots)>0:
            solutions.extend([b,roots,0])
        b+=1
  #  print("solutions: "+str(solutions))
    if prime == 2:
        max_lift=15
    elif prime == 3:
        max_lift=8
    elif prime == 5:
        max_lift=5
    elif prime == 7:
        max_lift=4
    else:
        max_lift=3
    exp+=1
    while exp < max_lift:
        sqr_lifted=lift_root2([1,0,-a],sqr,prime,exp)
        am_to_lift=0
        new_solutions=[]
        i=0
        while i < len(solutions):
            b=solutions[i]
            while b < prime**exp:
                #print("b: "+str(b)+" prime: "+str(prime)+" exp: "+str(exp)+" prime**exp: "+str(prime**exp))
                poly=[1,(b*sqr_lifted)%prime**exp,(-n*k)%prime**exp]
                der=get_derivative(poly)
                #poly=[1,-b,n*k]
                new_roots=[]
                hit=0
                if solutions[i+2]==0:
                    roots=solutions[i+1]
                    j=0
                    while j < len(roots):
                        r=roots[j]
                 #   print("r: "+str(r)+" b: "+str(b))
                        while r < prime**exp:
                            if evaluate(poly,r)%prime**exp ==0:
                                new_roots.append(r)
                                deriv=evaluate(der,r)
                                r0=find_r(prime,deriv)
                                ceiling=(prime*r0)
                                ceiling=prime**ceiling
                                if prime**exp > ceiling:
                                    hit=1
                              #  print("******************prime**exp: "+str(prime**exp)+" b: "+str(b)+" r: "+str(r)+" "+str(deriv)+" ceiling: "+str(ceiling))
                           # else:
                             #   print("prime**exp: "+str(prime**exp)+" b: "+str(b)+" r: "+str(r)+" "+str(deriv)+" ceiling: "+str(ceiling))
                            r+=prime**(exp-1)
                        j+=1
                    if len(new_roots)>0:
                    #print("b: "+str(b)+" "+str(len(new_roots)))
                        if hit==1:
                            new_solutions.extend([b,exp,1])
                        else:
                            new_solutions.extend([b,new_roots,0])
                            am_to_lift+=1

                else:
                    new_solutions.extend([solutions[i],solutions[i+1],solutions[i+2]])
                    break
                #if solutions[i+2]==1 and len(roots)==0:
                #    print("fatal error")
                #    sys.exit()
                #if solutions[i+2]==1 and hit==0:
                #    print("fata; error2")
                #    sys.exit()



                b+=prime**(exp-1)
            i+=3
        solutions=new_solutions
       # print("solutions: "+str(prime**exp)+" "+str(len(solutions)//3))
      #  if am_to_lift==0:
          #  print("!!!!!!!!!!!!!!!!FULLY LIFTED AT: "+str(exp))
          #  break
        if exp+1 == max_lift:
            break
        exp+=1

    
    blist.append(prime**(exp))
    blist.append([])
    i=0
    while i < len(solutions):
        if solutions[i+2]==1:
            texp=solutions[i+1]
            tempb=solutions[i]
            while tempb < prime**(exp):
                blist[-1].append(tempb)
                ##To do: can add a verification loop here later if I'm bored to make sure everything added is correct

                tempb+=prime**texp
        else:    
            blist[-1].append(solutions[i])

        i+=3
    blist[-1].sort()
    #print("solutions len: "+str(len(blist[-1]))+" k: "+str(k))

    return blist

def enumerated_product(*args):
    yield from itertools.product(*(range(len(x)) for x in args))

def psieve_build_interval(sbase,n,k,mod,root,a):
    interval=array.array("i",[1]*1_000)
  #  interval=np.ones(10_000,dtype=np.uint8)
    i=0
    while i < len(sbase):
        prime=sbase[i]
        if mod%prime==0:
            i+=1
            continue
        j=0
        while j <prime:
            disc=a*j**2+4*n*k
            if kronecker_symbol(disc,prime)==-1:
                
                dist=solve_lin_con(mod,j-root,prime)
             #   diff=(dist-start_ind)%prime
             #   dist2=start_ind+diff
              #  if dist2%prime != dist:
             #       print("fatal error")

                while dist < len(interval):
                    interval[dist]=0
                    dist+=prime
            j+=1
        i+=1
    return interval
   
def psieve(n,ret_array,primelist_f,fbase,a,sbase):#(n,fbase,div,hmap2,ret_array):
    sbase_trunc=20
    found=0
    primes_to_check=[]
    for prime in fbase:
        if a%prime==0:
            primes_to_check.append(prime)

    k=1
    while k < 1000: #to do: can also just precalculate residues here... but probably want to consider mostly small-ish k values...
        if math.gcd(a,k)!=1 or isPrime(k,5)!=1: #to do: does not need to be prime.. just need to avoid squares in the factorization... fix later
            k+=1
            continue
        
        skip=0
        for prime in primes_to_check:
            if kronecker_symbol(4*n*k,prime)==-1:
                skip=1

        if skip == 0:
            #print("[i]Checking k: "+str(k))
            sol_mod=1
            solutions=[]
            pcan=3 ##to do: prime=2 is most powerful but the find_roots_poly(poly,pcan) wont work on it.
            total_combo=1
            primes_added=[]
            while bitlen(sol_mod)<keysize//3 and pcan < 40:
                
                if isPrime(pcan,5)==1 and a%pcan !=0:
                    poly=[1,0,-a]
                    sqr=find_roots_poly(poly,pcan)
                    if len(sqr)>0:
  
           # max_filter=100
            ##THIS I WILL REFER TO AS THE FILTER AND WE WILL SET THE STEP SIZE FOR INTERVAL TO THIS!
                   # print("building pcan: "+str(pcan))
                        temp_solutions=brute_force_padic_solutions(pcan,k,n,a,sqr[0]) ##Using this as a filter... got to expand on this concept and add actual hensel too..
                       # print("pcan: "+str(pcan)+" len(temp_solutions[1]): "+str(len(temp_solutions[1])))
                        density=(temp_solutions[0]/len(temp_solutions[1]))
                    
                        if temp_solutions != -1 and len(temp_solutions[-1]) > 0 and density>3:
                      #  print("density: "+str(density)+" prime: "+str(pcan)+" prime^e: "+str(temp_solutions[0]))
                            solutions.extend(temp_solutions)
                            sol_mod*=temp_solutions[0]
                            total_combo*=len(temp_solutions[1])
                            primes_added.append(pcan)
                  #  solutions2=brute_force_padic_solutions2(pcan,k,n,a) ##Using this as a filter... got to expand on this concept and add actual hensel too..
                  #  if temp_solutions != solutions2:
                  #      print("solutions: "+str(solutions))
                  #      print("solutios2: "+str(solutions2))
                  #      print("wtf")
                  #      sys.exit()

                pcan+=1
            if bitlen(sol_mod)<keysize//3:# or total_combo > 10_000:
                k+=1
                continue

            solutions=get_partials(sol_mod,solutions)
            if 1:

                mod=1
                enum=[]
                total_combo=1
                i=0
                while i < len(solutions):
                    mod*=solutions[i]
                    enum.append(solutions[i+1])
                    total_combo*=len(solutions[i+1])
                    i+=2
                if mod != sol_mod:
                    print("catasrophic error1: "+str(a))
                    sys.exit()
                #print("enum: "+str(enum))
                print("total_combo: "+str(total_combo)+" mod: "+str(mod))
                for idx in enumerated_product(*enum):

                    root=0
                    i=0
                    while i < len(idx):
                        ind=idx[i]
                        root+=enum[i][ind]
                        i+=1
                    root%=mod
                    diff=abs((4*n*k)//a)
                    diff=math.isqrt(diff)
                    rdiff=(root-diff)%mod
                    diff+=rdiff
                    if diff%mod != root:
                        print("fatal error")
                        sys.exit()
                    bstart=diff-(mod*500)

                    interval=psieve_build_interval(sbase[:sbase_trunc],n,k,mod,root,a)
                    indexlist=np.nonzero(interval)[0]

           # print("Checking lin: "+str(lin)+" quad: "+str(quad_can)+" cmod: "+str(cmod)+" u2: "+str(u2)+" u: "+str(u)+" temp: "+str(temp))
                  #  print("k: "+str(k)+" root: "+str(root)+" mod: "+str(mod))
                  #  ind=0
                  #  length=len(indexlist)
                  #  print(indexlist)
                 #   while ind < length:# length:  
                  #      i=int(indexlist[ind])
                       # print("Found one at index: "+str(i))
                  #  print(interval)
                    nsqr=0
                    i=0
                    while i < len(interval):
                        if interval[i]==0:
                            i+=1
                            continue
                       # print("hallo??")
                        disc_otherside=a*(root+mod*i)**2+4*n*k 
                        krons=[]
                        for sprime in sbase:
                          #  if math.gcd(sprime, a)!=1 or mod%sprime ==0:
                          #      continue
                            #a_inv=modinv(a,sprime)

                          #  disc_otherside=(a*(root+mod*i)**2+4*n*k)%sprime    
                          #  disc_otherside*=a_inv
                           # disc_otherside%=sprime
                            sym=kronecker_symbol(disc_otherside,sprime)
                            if sym==-1:
                                nsqr+=1
                            krons.append(sym)#==-1:
                               # print("super catastrophic error core logic went wrong: "+str(sprime))
                               # sys.exit()
                         #   krons.append(kronecker_symbol(disc_otherside,sprime))
                       # if nsqr < 10:
                      #      print("symbols: "+str(krons))
                        disc_otherside=a*(root+mod*i)**2+4*n*k 
                        for prime in primes_added:
                            if kronecker_symbol(disc_otherside,prime)==-1:
                                print("catastrophic error: "+str(prime)+" a: "+str(a)+" primes added: "+str(primes_added))
                                sys.exit()
                      #  if i==500:
                       #     print(str(bitlen(disc_otherside))+" i: "+str(i))#+" k: "+str(k)+" a: "+str(a)+" mod: "+str(mod)+" krons: "+str(krons))


                      #  if disc_otherside%mod_otherside !=0:
                      #      print('fatal error')
                      #      sys.exit()
                        if 1:#disc_otherside%a == 0 and disc_otherside > 0:
                           # disc_otherside//=a
                            test=math.isqrt(abs(disc_otherside))





                            if test**2 == disc_otherside:# and (root+mod_otherside*i) != o_b:
                                ##To do: use blist for marking an interval.. we can use the small prime that we lifted as step size
                               # t=0
                               # while t < len(blist_otherside): ##I'm still thinking on how to incorporate this otherside.. some meet in the middle type algo? I odn't know..
                               #     if test%blist_otherside[t] not in blist_otherside[t+1]:
                               #         print("SUPER FATAL ERROR!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ")#+str(blist[t])+" "+str(test%blist[t])+" test: "+str(test))
                               #         sys.exit()
                               #     t+=2
                                #print("test: "+str(test)+" root: "+str((root+mod_otherside*i))+" mod_otherside: "+str(mod_otherside))
                                #print("found i (blist): "+str(i)+" new_root (otherside): "+str(new_root)+" mod_otherside: "+str(mod_otherside))
                                print("****************************************************************************Found one with psieve: "+str(test)+" k: "+str(k)+" a: "+str(a)+" interval index: "+str(i)+" interval[i]: "+str(interval[i]))#,krons)
                                new_root=a*(root+mod*i)
                                poly_val=(new_root)**2+4*n*k*a 
                                local_factors, value = factorise_fast(poly_val,primelist_f)

                                ret_array[1].append(new_root**2)
                                ret_array[0].append(poly_val)
                                ret_array[2].append(local_factors)
                                ret_array[3].append([])
                                found+=1


                        i+=1
                      #  ind+=1
        if found > 1:
            return found #should be enouhg..

        k+=1

    return found

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

def get_derivative(f):
    res = [0]*(len(f)-1)
    for i in range(len(f)-1):
        res[i] = (len(f)-1-i)*f[i]
    return res

def lift_root2(coeffs, root, p, k):

   # coeffs = list(reversed(coeffs))
    deriv_coeffs = get_derivative(coeffs)

    fprime_at_root = evaluate(deriv_coeffs, root)%p
    r=root%p
    modulus=p
    for _ in range(1, k):
        modulus*=p
        f_val = evaluate(coeffs, r)%modulus
        fp_val = evaluate(deriv_coeffs, r)%modulus
        fp_inv = modinv(fp_val,modulus)
        r = (r - f_val * fp_inv) % modulus
    return r

def get_primes(start,stop):
    return list(sympy.sieve.primerange(start,stop))


def main(l_keysize,l_workers,l_debug,l_base,l_key,l_lin_sieve_size,l_quad_sieve_size):
    global key,keysize,workers,g_debug,base,key,lin_sieve_size,quad_sieve_size,max_bound,lin_sieve_size2
    key,keysize,workers,g_debug,base,lin_sieve_size,quad_sieve_size=l_key,l_keysize,l_workers,l_debug,l_base,l_lin_sieve_size,l_quad_sieve_size
    lin_sieve_size2=lin_sieve_size
    start = default_timer() 
    if max_bound > lin_sieve_size:
        max_bound = lin_sieve_size
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
        if n%primeslist[i] !=0:
            primeslist1.append(primeslist[i])
        i+=1
    primeslist2.append(2)
    i=0
    while len(primeslist2) < 100:
        if n%primeslist[i] !=0:
            primeslist2.append(primeslist[i])
        i+=1
    launch(n,primeslist1,primeslist2)     
    duration = default_timer() - start
    print("\nFactorization in total took: "+str(duration))