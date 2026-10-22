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
import multiprocessing
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
    a_mul=1
    a_mul_list=[]
    while a_mul< 10_000:
        amul=a_mul
        
        a_mul_factors=[]
        for prime in primeslist:
            while amul%prime==0:
                amul//=prime
                if prime not in a_mul_factors:
                    a_mul_factors.append(prime)
        if amul==1:
            a_mul_list.append([a_mul,a_mul_factors])
        a_mul+=2
    hmap=create_hashmap(n,primeslist)
    ret_array=[[],[],[],[]]
    partials={}
    large_prime_bound = primeslist[-1] ** lp_multiplier
    sbase=copy.deepcopy(primeslist[0:30])
    resmaps,resmaps2=build_residues(sbase,n)
   # print("fbase: "+str(fbase))
    print("[i]Building psieve Residue Map (to do: some duplication here from merging two algos, fix later)")

    #sys.exit()

    grays = get_gray_code(20)
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
    #opt values:
    #close_range=10
    #too_close=5
    #LOWER_BOUND_SIQS=400
    #UPPER_BOUND_SIQS=4000

    close_range=20
    too_close=1
    LOWER_BOUND_SIQS=1
    UPPER_BOUND_SIQS=4000
   # tnum=int(((n)**0.5) /(lin_sieve_size))
    tnum=int(((n)**0.5) / 1)
    seen=[]

    print("[i]Building 2d root map")
    roots2d=build_2drootmap(primeslist_a,hmap,n)
    print("[i]Entering attack loop")
    valid_quads,valid_quads_factors=filter_quads(qlist,n)
    fb_map = {val: i for i, val in enumerate(primeslist)}
   # print("fb_map: ",fb_map)
   # factor_ranking=[]#np.zeros(len(primeslist),dtype=np.uint16)
    seen_factors=array.array('i',len(primeslist)*[0])
    retry=0
    while 1:
        factor_ranking=[]
        quad=1
        new_mod,cfact,indexes=generate_modulus(n,primeslist,seen,tnum,close_range,too_close,LOWER_BOUND_SIQS,UPPER_BOUND_SIQS,bitlen(tnum),quad)
       # print("mod: "+str(new_mod)+" cfact: "+str(cfact)+" indexes: "+str(indexes))


      #  new_mod=961
      #  cfact=[961]
      #  indexes=[9]
       # new_mod=37**2
       # cfact=[37**2]
       # indexes=[10]
        if new_mod ==0:
            retry+=1
            if retry > 5:
                print("failed to generate modulus..")
                return 0,0
            continue

        retry=0
        lin,lin_parts=get_lin(cfact,new_mod,quad,n,1)
        z=quad#quadlist[j]
        lin_co_array=[]
        q=0
 
        lin2=lin#lin3%new_mod
        poly_ind=0
        end = 1 << (len(cfact) - 1)
        lin=0
        while poly_ind < end:
            if poly_ind != 0:
                v,e=grays[poly_ind]
                lin=(lin + 2 * e * lin_parts[v])%new_mod
            else:
                lin=lin2
            if quad_sign == "neg":
                if (z*lin**2-n)%new_mod !=0:
                    print("super big error")
                    sys.exit()
            else:
                if (z*lin**2+n)%new_mod !=0:
                    print("super big error")
                    sys.exit()            
            interval=build_database2interval(primeslist_a,quad,n,lin,new_mod,roots2d,0,factor_ranking)
            found+=process_interval2d(n,ret_array,quad,primelist_f,large_prime_bound,partials,lin,new_mod,factor_ranking,fb_map,0,seen_factors,interval,primeslist,resmaps,resmaps2,valid_quads,valid_quads_factors,qlist,primelist,sbase,a_mul_list)#,lin,new_mod,sum_list)
           # if found > 100 or len(ret_array[0]) > base+10:
           #     if g_debug ==1:
           #         print("seen_factors: ",seen_factors)
           #     print("[i]Performing linear algebra")
           #     test,test2=QS(n,primelist,ret_array[0],ret_array[2],ret_array[1],ret_array[3])
      #      test,test2=QS(n,primelist,ret_array[0],ret_array[1],ret_array[2]) 
           #     found=0 
           #     if test !=0:
           #         print("\n\n\n\nFound at: ",len(ret_array[0]))
           #         return 

            poly_ind+=1
      
     
        ######To do: Return a list of all odd exponent factors for any smooths found. Now iterate that and call gen_modulus_and_interal for those quadratic coefficients.
    

    return 


def equation(y,x,n,mod,z,z2):
    rem=z*(x**2)-y*x+n*z2
    rem2=rem%mod
    return rem2,rem 

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


def tonelli(n,p):  # tonelli-shanks to solve modular square root: x^2 = n (mod p)
    q = p - 1
    s = 0

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

def lift_b(prime,n,co,z,max_exp):
    z2=1
    k=0
    ret=[]
    cos=[]
    step_size=[]
    new=[]
    r=get_root(prime,co%prime,z) 
    if r==-1:
        return 0
    exp=2
    ret=[co,r]
    cos.append([co,(prime-co)%prime])
    while exp < max_exp+1:

        ret=lift(exp,ret[0],ret[1],n,z,z2,prime)
        co2=(prime**exp)-ret[0]
        cos=([ret[0],co2])

        exp+=1
    return cos[0]

def solve_roots(prime,n): 
    if quad_sign == "neg":
        s=1  
        while jacobi((s*n)%prime,prime)!=1:
            s+=1
    else:
        s=1  
        while jacobi((-s*n)%prime,prime)!=1:
            s+=1
    z_div=modinv(s,prime)  
    dist=(n*z_div)%prime
    if quad_sign != "neg":
        dist=(-dist)%prime
    main_root=tonelli(dist,prime)
    if main_root**2%prime != dist:
        print("what the fuck")

    if quad_sign == "neg":
        if (s*main_root**2-n)%prime !=0:
            print("fatal error123: "+str(prime)+" s: "+str(s)+" root: "+str(main_root))
            sys.exit()
    else:
        if (s*main_root**2+n)%prime !=0:
            print("fatal error123: "+str(prime)+" s: "+str(s)+" root: "+str(main_root))
            sys.exit()
    try:
     #   size=prime*2+1
     #   if size > quad_sieve_size*2:
     #       size= quad_sieve_size*2+1
        size=3
        temp_hmap = array.array('I',[0]*size) ##Got to make sure the allocation size doesn't overflow.... 
        temp_hmap[0]=1

        s_inv=modinv(s*z_div,prime)
        if s_inv == None or jacobi(s_inv,prime)!=1:
            print("should this ever happen?")
            return temp_hmap
        root_mult=tonelli(s_inv,prime)
        new_root=((main_root*root_mult))%prime
        if quad_sign == "neg":
            if (s*new_root**2-n)%prime !=0:
                print("error2")
        else:
            if (s*new_root**2+n)%prime !=0:
                print("error2")
       # new_co=(2*s*new_root)%prime
       # if (new_co**2+n*4*s)%prime !=0:
           # print("error")
      #  if (s*new_root**2-new_co*new_root+n)%prime !=0: ###To do: For debug delete later
           # print("error")
        if new_root > prime // 2:
            new_root=(prime-new_root)%prime  

        end=temp_hmap[0]
        temp_hmap[end]=s
        if bitlen(new_root)>32:
            print("fatal error, increase element size of array in solve_roots")
            sys.exit()
        temp_hmap[end+1]=new_root
        temp_hmap[0]+=2
        s+=1   
    except Exception as e:
        print(e)
    return temp_hmap


def create_hashmap(n,primeslist):
    i=0
    hmap=[]
    while i < len(primeslist):
        hmap_p=solve_roots(primeslist[i],n)
        hmap.append(hmap_p)
        i+=1
    return hmap



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

cdef process_interval2d(n,ret_array,quad_can,primelist_f,large_prime_bound,partials,lin,cmod,factor_ranking,fb_map,bSeenOnly,seen_factors,interval,primeslist,resmaps,resmaps2,valid_quads,valid_quads_factors,qlist,primelist,sbase,a_mul_list):#,lin,cmod,sum_list):
    linsize=lin_sieve_size
    if bSeenOnly==1:
        linsize=lin_sieve_size2
        threshold = int(math.log2((linsize)*math.sqrt(abs(n))) - thresvar2)
    else:
        threshold = int(math.log2((linsize)*math.sqrt(abs(n))) - thresvar)
    #threshold = int(math.log2(abs(n)) - thresvar)

    found=0
   
    u=0
    cdef Py_ssize_t k
        
    np.putmask(interval, interval<threshold, 0)
    indexlist=np.nonzero(interval)[0]
    temp=interval
    root=lin
           # print("Checking lin: "+str(lin)+" quad: "+str(quad_can)+" cmod: "+str(cmod)+" u2: "+str(u2)+" u: "+str(u)+" temp: "+str(temp))

    k=0
    length=len(indexlist)
    while k < length:# length:  
        i=indexlist[k]
        if temp[i]>threshold:       
            x=root+cmod*int(i)
            if quad_sign=="neg":
                poly_val=quad_can*x**2-n
            else: 
                poly_val=quad_can*x**2+n
            if poly_val%cmod !=0:
                print("ERGH")
                time.sleep(1000)                    
            new_root=quad_can*x
            if quad_sign=="neg":
                poly_val=new_root**2-n*quad_can
            else: 
                poly_val=new_root**2+n*quad_can
            if poly_val%(cmod*quad_can)!=0:
                print("fatal")
            local_factors, value = factorise_fast(poly_val,primelist_f)
                #if g_debug ==1: note: need to fix is we skip odd exponents for large primes when lifting
                    #########START DEBUG###############
                    #poly_val2=quad_can*x**2+n
                    #local_factors2, value2,seen_primes2 = factorise_fast_debug(poly_val2//cmod,primelist_f,g_max_exp)
                    
                   # seen_log=0
                  #  for prime in seen_primes2:
                       # if cmod%prime !=0 and quad_can%prime !=0:
                          #  seen_log+=round(math.log2(prime))
 
        
                #    if seen_log != temp[i]:
                     #   print("error:"+str(seen_primes2)+" seen_log: "+str(seen_log)+" assumed log: "+str(temp[i])+" quad: "+str(quad_can)+" cmod: "+str(cmod))
                    #########END DEBUG###############
            #if value != 1:
            #    if value < large_prime_bound:
            #        if value in partials:
            #            rel, lf, pv = partials[value]
            #            if rel == new_root:
            #                k+=1
            #                continue
            #            new_root *= rel
            #            local_factors ^= lf
            #            poly_val *= pv
            #        else:
            #            partials[value] = (new_root, local_factors, poly_val)
            #            k+=1
            #            continue
            #    else:
            #        k+=1 
            #        continue     
            if value != 1 and isPrime(value,5)==0:#math.isqrt(abs(value))**2 != value:
                k+=1
                continue    
            if new_root not in ret_array[1]:
                factor_ranking.append([])
                local_factors=list(local_factors)
                local_factors.sort()
    
                for fac in local_factors:
                    if fac != -1 and fac !=2:
                        idx=fb_map[fac]
                        seen_factors[idx]+=1
                        if bSeenOnly==0:
                               
                            factor_ranking[-1].append(idx)
                if g_debug == 1:
                    if bSeenOnly==0:
                        print("***seen_primes: "+str(local_factors)+" cmod: "+str(cmod))
                    if bSeenOnly==1:
                        print("seen_primes: "+str(local_factors)+" cmod: "+str(cmod))
                found+=1

                #To do: uncomment later
                #ret_array[1].append(new_root**2)
                #ret_array[0].append(poly_val)
                #ret_array[2].append(local_factors)
                #ret_array[3].append([])
                div_fac=[]
                faclist=list(local_factors)
                faclist.sort()
                faclist.reverse()
             
              #  print("faclist: "+str(faclist))
                div=1
                for odd_exp_factor in faclist:
                   # if odd_exp_factor != -1:
                    div*=odd_exp_factor
                    div_fac.append(odd_exp_factor)
                
                #div*=-1
                #    break
             #   div*=3
               # div_fac.append(3)
           
                #To do: Sieve around the coefficient the b-smooth was found at
               # print("[*]Smooths: "+str(len(ret_array[0]))+" / "+str(base)+" b: "+str(new_root)+" k: "+str(quad_can))#+" square: "+str((abs(poly_val//div))**0.5))
               # if bitlen(div)<keysize*0.50: ##Dont know if this matters.. another parameter to test with..
                  #  print("PSIEVE1")
                local_factors2, value2 = factorise_fast(new_root,primelist_f)
                #To do: Fix this for when poly_val is smaller then 0... for some reason my calculations dont always hold true in that case
                if value==1 and div!=1 and div%2!=0 and bitlen(div)<(keysize*0.4) and poly_val >0:# and isPrime(value,5)==1:# and isPrime(div,5)==1:# and value2==1:# and len(div_fac)==1:
                    
                    print("[i]Trying psieve b: "+str(2*new_root)+" a: "+str(div)+" bitlen a: "+str(bitlen(div)))
                    psievefound=psieve(n,ret_array,primelist_f,primeslist,div,sbase,2*new_root,resmaps,resmaps2,a_mul_list)
                    if psievefound !=0:
                        ret_array[1].append(new_root**2)
                        ret_array[0].append(poly_val)
                        ret_array[2].append(local_factors)
                        ret_array[3].append([])
                        print("[*](Psieve)Trying linear algebra after succesful psieve run")
                        test,test2=QS(n,primelist,ret_array[0],ret_array[2],ret_array[1],ret_array[3])
                    

                        if test !=0:
                            print("\n\n\n\nFound at: ",len(ret_array[0]))
                            sys.exit()
                #    sys.exit()
                    #break



                  #  sys.exit()
            #    print("***seen_primes: "+str(local_factors)+" cmod: "+str(cmod)+" root: "+str(new_root)+" polyval: "+str(poly_val))
                if len(ret_array[0])>(base+10):
                    return found
        k+=1

  

    if g_debug ==1:
        if bSeenOnly == 1:
            print("found: "+str(found)+" total: "+str(len(ret_array[0]))+" cmod bits: "+str(bitlen(cmod))+" quad: "+str(quad_can))
        else:
            print("***found: "+str(found)+" total: "+str(len(ret_array[0]))+" cmod bits: "+str(bitlen(cmod))+" quad: "+str(quad_can))
    return found

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
        exp=0
        factor=factor_base[i]
        while value % factor == 0:
            exp+=1
            factors ^= {factor}
            value //= factor
        i+=1
        if exp !=0 and exp%2==0:
            return -1,-1
    return factors, value


def filter_quads(qbase,n):
    ###Note: We look for quadratic coefficients that factor over the factor base but have no even exponents. This garantuees unique results
    valid_quads=[]
    valid_quads_factors=[]

    roots=[]
    i=1#2**(keysize//2)
    while len(valid_quads) < quad_sieve_size+1:

        quad_local_factors, quad_value = factorise_fast(i,qbase) 
        if quad_value != 1:
            i+=1
            continue
       # print("found quad: "+str(len(valid_quads)))
        valid_quads.append(i)
        valid_quads_factors.append(quad_local_factors)
        i+=1
    print("valid_quads: "+str(valid_quads))
    return valid_quads,valid_quads_factors


#@cython.boundscheck(False)
#@cython.wraparound(False)
def generate_modulus(n,primeslist,seen,tnum,close_range,too_close,LOWER_BOUND_SIQS,UPPER_BOUND_SIQS,tnum_bit,quad):
    const_1=1_000
    const_2=10_000

    small_B = base#len(primeslist)
    lower_polypool_index = 2
    upper_polypool_index = small_B - 1
    poly_low_found = False
    
    for i in range(small_B):  ##To do: Can be moved outside mainloop
        if primeslist[i]**2 > LOWER_BOUND_SIQS and not poly_low_found:
            lower_polypool_index = i
            poly_low_found = True
            break
        if primeslist[i]**2 > UPPER_BOUND_SIQS:
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
        while counter2 < const_1:
            counter2+=1
            found_a_factor = False
            counter=0
            while(found_a_factor == False) and counter < const_2:
                randindex = random.randint(lower_polypool_index, upper_polypool_index)
                if quad_sign == "neg":
                    if  jacobi((quad*n)%primeslist[randindex],primeslist[randindex])!=1:
                        counter+=1
                        continue
                else:
                    if  jacobi((-quad*n)%primeslist[randindex],primeslist[randindex])!=1:
                        counter+=1
                        continue
                potential_a_factor = primeslist[randindex]**2
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
            if quad_sign == "neg":
                if  jacobi((quad*n)%primeslist[randindex],primeslist[randindex])!=1:#hmap[randindex][1]!=quad%primeslist[randindex]:
                    print("THE FUC")
                    time.sleep(1000000)
            else:
                if  jacobi((-quad*n)%primeslist[randindex],primeslist[randindex])!=1:#hmap[randindex][1]!=quad%primeslist[randindex]:
                    print("THE FUC")
                    time.sleep(1000000)                    
            indexes.append(randindex)
            j = tnum_bit - cmod.bit_length()
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
            if abs(a1 - primeslist[i]**2) < mindiff:
                randindex = i
                mindiff = abs(a1 - primeslist[i]**2)
                
        

        found_a_factor = False
        counter3=0
        while not found_a_factor and counter3< const_1 and randindex <base:
            if quad_sign=="neg":
                if  jacobi((quad*n)%primeslist[randindex],primeslist[randindex])!=1:
                    randindex += 1
                    counter3+=1
                    continue
            else:
                if  jacobi((-quad*n)%primeslist[randindex],primeslist[randindex])!=1:
                    randindex += 1
                    counter3+=1
                    continue                    
            potential_a_factor = primeslist[randindex]**2

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
        if quad_sign=="neg":
            if  jacobi((quad*n)%primeslist[randindex],primeslist[randindex])!=1:
                print("THE FUC: ",randindex)
                time.sleep(1000000)
        else:
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

#@cython.boundscheck(False)
#@cython.wraparound(False)
def generate_modulus2(n,primeslist,seen,tnum,close_range,too_close,LOWER_BOUND_SIQS,UPPER_BOUND_SIQS,tnum_bit,quad):
    const_1=1_000
    const_2=1_000_000

    small_B = base#len(primeslist)
    lower_polypool_index = 2
    upper_polypool_index = small_B - 1
    poly_low_found = False
    
    for i in range(small_B):  ##To do: Can be moved outside mainloop
        if primeslist[i]**2 > LOWER_BOUND_SIQS and not poly_low_found:
            lower_polypool_index = i
            poly_low_found = True
            break
        if primeslist[i]**2 > UPPER_BOUND_SIQS:
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
                potential_a_factor = primeslist[randindex]**2
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
            if abs(a1 - primeslist[i]**2) < mindiff:
                randindex = i
                mindiff = abs(a1 - primeslist[i]**2)
                
        

        found_a_factor = False
        counter3=0
        while not found_a_factor and counter3< const_2 and randindex <base:
            if  jacobi((-quad*n)%primeslist[randindex],primeslist[randindex])!=1:
                randindex += 1
                continue
            potential_a_factor = primeslist[randindex]**2

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

def get_lin(cfact,local_mod,k,n,a):
    all_lin_parts=[]
    j=0
    lin=0

    while j < len(cfact):

        prime=math.isqrt(cfact[j])
        if quad_sign == "neg":
            try:
                r1=find_roots_poly([a,0,-n*k],prime)
                r1=r1[0]
                sq=[1,0,-a]
                sqr=find_roots_poly(sq, prime)
                r1=lift_root2([a,0,-n*k],r1,prime,2)#(r1,prime,-n,quad_co,2) #to do: use the cleaner lift_root2()
                if evaluate([a,0,-n*k],r1)%prime**2 !=0:
                    print("something fatal happened")
                    sys.exit()
            except Exception as e:
                print(str(" a: "+str(a)+" k: "+str(k))+" prime: "+str(prime))#+" r1: "+str(find_roots_poly([a,0,-n*4*k],prime))+" sqr: "+str(sqr))
                print(e)
                sys.exit()
        else:
            r1=find_roots_poly([a,0,n*k],prime)
            r1=r1[0]
            r1=lift_root2([a,0,n*k],r1,prime,2)
            if evaluate([a,0,n*k],r1)%prime**2 !=0:
                print("something fatal happened")
                sys.exit()
        prime=prime**2
        aq = local_mod // prime
        invaq = modinv(aq%prime, prime)
        gamma = r1 * invaq % prime
        lin+=aq*gamma
        all_lin_parts.append(aq*gamma)
        j+=1
    lin%=local_mod
   # print("all_lin_parts: "+str(all_lin_parts))
    return lin,all_lin_parts

def get_lin3(cfact,local_mod,k,n,a):
    ##Merge all these variants later on..
    all_lin_parts=[]
    j=0
    lin=0

    while j < len(cfact):

        prime=math.isqrt(cfact[j])
        if quad_sign == "neg":
            try:
                sq=[1,0,-a]
                sqr=find_roots_poly(sq, prime)
                r1=find_roots_poly([a,0,-n*4*k],prime)
                r1=r1[0]

                r1=lift_root2([a,0,-n*4*k],r1,prime,2)#(r1,prime,-n,quad_co,2) #to do: use the cleaner lift_root2()
                if evaluate([a,0,-n*4*k],r1)%prime**2 !=0:
                    print("something fatal happened")
                    sys.exit()
            except Exception as e:
                print(str(" a: "+str(a)+" k: "+str(k))+" prime: "+str(prime)+" r1: "+str(find_roots_poly([a,0,-n*4*k],prime))+" sqr: "+str(sqr))
                print(e)
                sys.exit()
        else:
            try:
                sq=[1,0,-a]
                sqr=find_roots_poly(sq, prime)
                r1=find_roots_poly([a,0,n*4*k],prime)
                r1=r1[0]

                r1=lift_root2([a,0,n*4*k],r1,prime,2)#(r1,prime,-n,quad_co,2) #to do: use the cleaner lift_root2()
                if evaluate([a,0,n*4*k],r1)%prime**2 !=0:
                    print("something fatal happened")
                    sys.exit()
            except Exception as e:
                print(str(" a: "+str(a)+" k: "+str(k))+" prime: "+str(prime)+" r1: "+str(find_roots_poly([a,0,n*4*k],prime))+" sqr: "+str(sqr))
                print(e)
                sys.exit()
        prime=prime**2
        aq = local_mod // prime
        invaq = modinv(aq%prime, prime)
        gamma = r1 * invaq % prime
        lin+=aq*gamma
        all_lin_parts.append(aq*gamma)
        j+=1
    lin%=local_mod
   # print("all_lin_parts: "+str(all_lin_parts))
    return lin,all_lin_parts

def get_lin2(cfact,local_mod,k,n,a):
    all_lin_parts=[]
    j=0
    lin=0

    while j < len(cfact):

        prime=cfact[j]
        if quad_sign == "neg":
            try:
                r1=find_roots_poly([a,0,-n*4*k],prime)
                sq=[1,0,-a]
                sqr=find_roots_poly(sq, prime)

                r1=r1[0]
            except Exception as e:
                print(str(" a: "+str(a)+" k: "+str(k))+" prime: "+str(prime)+" r1: "+str(find_roots_poly([a,0,-n*4*k],prime))+" sqr: "+str(sqr))
                print(e)
                sys.exit()

        else:
            r1=find_roots_poly([a,0,n*4*k],prime)
            r1=r1[0]


        aq = local_mod // prime
        invaq = modinv(aq%prime, prime)
        gamma = r1 * invaq % prime
        lin+=aq*gamma
        all_lin_parts.append(aq*gamma)
        j+=1
    lin%=local_mod
   # print("all_lin_parts: "+str(all_lin_parts)+" cfact: "+str(cfact)+" local_mod: "+str(local_mod))
    return lin,all_lin_parts

@cython.boundscheck(False)
@cython.wraparound(False)
cdef build_database2interval(long long [:] primeslist,quad,n,lin,cmod,roots2d,bSeenOnly,factor_ranking):
    

    cdef Py_ssize_t i
    if bSeenOnly==1:
        linsize=lin_sieve_size2
    else:
        linsize=lin_sieve_size
    interval_single=np.zeros(linsize,dtype=np.uint16)    

    i=0
    while i < len(primeslist):
        prime=primeslist[i]
        log=round(math.log2(prime))
        xc = int(roots2d[quad,i])
        if xc ==0:
            i+=1
            continue

        root=lin
        if cmod%prime == 0 or quad%prime ==0:
            i+=1
            continue
        lcmod=cmod%prime
        modi=modinv(lcmod,prime)
        x=xc
        exp=1
        while exp < g_max_exp+1:
                
            p=prime**exp
            if exp > 1:
                if quad_sign=="neg":
                    x=lift_root(x,prime**(exp-1),-n,quad,exp)#replace with lift_root2
                else:
                    x=lift_root(x,prime**(exp-1),n,quad,exp)
            root_dist1=solve_lin_con(cmod,x-root,p)
            x_b=(p-x)%p 
            root_dist2=solve_lin_con(cmod,x_b-root,p)
            if bSeenOnly==0 or (bSeenOnly ==1 and prime < dupe_max_prime):
                hit=0
                if root_dist1 < linsize+1:
                    interval_single[root_dist1::p]+=log
                    hit=1
                if root_dist2 < linsize+1:
                    if root_dist1 != root_dist2:
                        interval_single[root_dist2::p]+=log  
                        hit=1
                if quad_sign=="neg":
                    CAN=quad*(root+root_dist1*cmod)**2-n
                    if CAN % (cmod*p)  != 0:

                        print("EROROREROR",prime)
                        sys.exit()
                else:
                    CAN=quad*(root+root_dist1*cmod)**2+n
                    if CAN % (cmod*p)  != 0:

                        print("EROROREROR",prime)
                        sys.exit()                    
                if hit ==0:
                    break
            elif bSeenOnly ==1 and prime > dupe_max_prime-1 and exp%2 ==0:
                hit=0
                if root_dist1 < linsize+1:
                    interval_single[root_dist1::p]+=(log*2)
                    hit=1
                if root_dist2 < linsize+1:
                    if root_dist1 != root_dist2:
                        interval_single[root_dist2::p]+=(log*2)
                        hit=1
                if quad_sign=="neg":
                    CAN=quad*(root+root_dist1*cmod)**2-n
                    if CAN % (cmod*p)  != 0:

                        print("EROROREROR",prime)
                        sys.exit()
                else:
                    CAN=quad*(root+root_dist1*cmod)**2+n
                    if CAN % (cmod*p)  != 0:

                        print("EROROREROR",prime)
                        sys.exit()                    
                if hit ==0:
                    break

            exp+=1   
        i+=1
    return interval_single

def build_2drootmap(primeslist,hmap,n):       
    roots2d=np.zeros([quad_sieve_size+1,base],dtype=np.int32)
    i=0
    while i < len(primeslist):
        prime=primeslist[i]
        j=0
        while j < prime and j < 1+1: ###To do: Only go up to prime
            quad=j#+offset
            #print("Building quad: "+str(quad))
            z=hmap[i][1]
            z_div=modinv(z,prime)
            z_inv=modinv(quad*z_div,prime)
            if z_inv == None or jacobi(z_inv,prime)!=1:
                j+=1
                continue  
            
            

            x=hmap[i][2]
            root_mult=tonelli(z_inv,prime)
            x=(x*root_mult)%prime
            if bitlen(x)>32:
                print("big root, increase dtype in build_2drootmap()")
            roots2d[quad::prime,i]=x

            if quad_sign == "neg":
                if (quad*x**2-n)%prime !=0:
                    print("fatal error, sad face :(")
                    sys.exit()
            else:
                if (quad*x**2+n)%prime !=0:
                    print("fatal error, sad face :(")
                    sys.exit()                
            j+=1
        i+=1
    i=0
    return roots2d

def find_r(mod,total):
    mo,i=mod,0
    while (total%mod)==0:
        mod=mod*mo
        i+=1
    return i


def brute_force_padic_solutions2(prime,k,n,a,sqr,blift):
    ##to do: delete later.. just for verifications
    solutions=[]
    blist=[]
    exp=1
    b=0
    while b < prime**exp:
       # print("b: "+str(b)+" prime: "+str(prime))
        if n > 0:
            poly=[1,-b*sqr,n*k]
        else:
            poly=[1,b*sqr,n*k]
        #poly=[1,-b,n*k]
        roots=[]
        x=0
        while x < prime**exp:
            if evaluate(poly,x)%prime**exp ==0:
                roots.append(x)

            x+=1
        if len(roots)>0:
            solutions.extend([b,roots])
       # if len(roots)==1:
       #     return -1
        b+=1
   # print("solutions: "+str(solutions))
    if prime == 2:
        max_lift=12
    elif prime == 3:
        max_lift=8
    elif prime == 5:
        max_lift=5
    elif prime == 7:
        max_lift=3
    elif prime < 40:
        max_lift=2
    else:
        max_lift=1
    if blift==0:
        max_lift=1
    if max_lift > 1:
        exp+=1
        while exp < max_lift:
            sqr_lifted=lift_root2([1,0,-a],sqr,prime,exp)
            if evaluate([1,0,-a],sqr_lifted)%prime**exp !=0:
                print("fatal error")
                sys.exit()
            new_solutions=[]
            i=0
            while i < len(solutions):
                b=solutions[i]
                while b < prime**exp:
                #print("b: "+str(b)+" prime: "+str(prime)+" exp: "+str(exp)+" prime**exp: "+str(prime**exp))
                    if n > 0:
                        poly=[1,(-b*sqr_lifted)%prime**exp,(n*k)%prime**exp]
                    else:
                        poly=[1,(b*sqr_lifted)%prime**exp,(n*k)%prime**exp]
                    der=get_derivative(poly)
                #poly=[1,-b,n*k]
                    new_roots=[]
      

                    roots=solutions[i+1]
                    j=0
                    while j < len(roots):
                        r=roots[j]
                 #   print("r: "+str(r)+" b: "+str(b))
                        while r < prime**exp:
                            if evaluate(poly,r)%prime**exp ==0:
                                new_roots.append(r)
                            r+=prime**(exp-1)
                        j+=1
                    if len(new_roots)>0:

                        new_solutions.extend([b,new_roots])
                    b+=prime**(exp-1)
                i+=2
            solutions=new_solutions
            if exp+1 == max_lift:
                break
            exp+=1

    
    blist.append([prime,exp])
    blist.append([])
    i=0
    while i < len(solutions):

        blist[-1].append(solutions[i])

        i+=2
    blist[-1].sort()
   # print("blist: "+str(blist))
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
        poly=[1,-b*sqr,n*k]
        #poly=[1,-b,n*k]
        roots=[]
        x=0
        while x < prime**exp:
            if evaluate(poly,x)%prime**exp ==0:
                roots.append(x)
                #disc=a*(b**2)-4*n*k*a
                #disc//=a

            x+=1
        if len(roots)>0:
            solutions.extend([b,roots,0])
        b+=1
  #  print("solutions: "+str(solutions))
    if prime == 2:
        max_lift=12
    elif prime == 3:
        max_lift=8
    elif prime == 5:
        max_lift=5
    elif prime == 7:
        max_lift=3
    elif prime < 20:
        max_lift=2
    else:
        max_lift=1
    max_lift=1
    if max_lift > 1:
        exp+=1
        while exp < max_lift:
            sqr_lifted=lift_root2([1,0,-a],sqr,prime,exp)
            if evaluate([1,0,-a],sqr_lifted)%prime**exp !=0:
                print("fatal error")
                sys.exit()
            am_to_lift=0
            new_solutions=[]
            i=0
            while i < len(solutions):
                b=solutions[i]
                while b < prime**exp:
                #print("b: "+str(b)+" prime: "+str(prime)+" exp: "+str(exp)+" prime**exp: "+str(prime**exp))
                    poly=[1,-b*sqr_lifted,n*k]
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




def psieve_build_interval(resmaps,n,k,root,a,sbase,primes_to_mark,sqr_list):
    interval=array.array("i",[1]*50_000)

  #  interval=np.ones(10_000,dtype=np.uint8)
    i=0
    while i < len(primes_to_mark):
        ind=primes_to_mark[i]
        prime=sbase[ind]#blist_otherside2[i][0]**blist_otherside2[i][1]
        if a%prime == 0 or k%prime==0:
            i+=1
            continue

     #   sqr2=sqr_list[i] 
        sqr=find_roots_poly([1,0,-a], prime) 
    #    sqr2.sort()
    #    sqr.sort()
       # if sqr2 != sqr:
          #  print("sqr2: "+str(sqr2)+" sqr: "+str(sqr))
        blist=resmaps[ind][k%prime]
           # print("blist: "+str(blist))

        rootlist=blist[1]
        p=0
        while p < prime:
            p2=(p*sqr[0])%prime
            if p2 not in rootlist:
                    
                dist=solve_lin_con(a,p2-root,prime)
                while dist < len(interval):
                    interval[dist]=0
                    dist+=prime

            p+=1
        i+=1
    return interval

def debug_find_residues3(prime,n,a,exp,k):
 #   hmap={}
   # klist=[prime**exp,[]]
    blist=[prime**exp,[]]
  
    b=0
    sq=[1,0,-a]
    sqr=find_roots_poly(sq, prime) 

    t=0
    while t < len(sqr):
        sqr_r=sqr[t]
   # for sqr_r in sqr:
   # sqr_r=sqr[0]
        sqr_r=lift_root2([1,0,-a], sqr_r, prime, 2)
   # print("sqr_r: "+str(sqr_r))
        while b < prime**exp:
            poly=[1,-b*sqr_r,n*k]
            polyc=copy.deepcopy(poly)
            roots=find_roots_poly(polyc,prime)
            if len(roots)<2:
                b+=1
                continue
      #  if b == 2 and len(roots)==2:
      #  print("roots: "+str(roots)+" b: "+str(b)+" prime**exp: "+str(prime**exp))
            if len(roots)>0: ##to do: singular??
                i=0
                while i < len(roots):
                #deriv=get_derivative(poly)
               # der=evaluate(deriv,roots[i])
               # poly2=[1,-der*sqr_r,n]
                    roots[i]=lift_root2(poly, roots[i], prime, exp)

                    if evaluate(poly,roots[i])%prime**exp!=0:
             #   if (roots[i]**2-4*n*k)%prime**exp !=0:
                        print("something screwed up: "+str(roots)+" k: "+str(k))
                        sys.exit(0)
                    disc1=a*b**2-4*n*k
                    if kronecker_symbol(disc1,prime)==-1:
                        print("kroneckerfail:"+str(prime))
                        sys.exit()
                    disc2=(a*b)**2-4*n*a*k 
                    if disc2%a!=0:
                        print("fatal")
                        sys.exit()
                    i+=1
           # klist[-1].append(k%prime**exp)
                if t==0:
                    blist[-1].append(b)
                if t == 1 and b not in blist[-1]:
                    print("WHAT THE FUCK")
                    sys.exit()
          #  hmap[k]=roots
           # hmap.append(roots)
   
            b+=1
        t+=1
   # print("prime: "+str(prime)+" hmap: "+str(blist))
    return blist

def psieve_padic_solutions(prime,exp,n,k):
    blist=[]
    blist2=[]
    blist.append([prime,exp])
    blist2.append([prime,exp])
    rootlist=[] ##To do: precalculate..................
    rootlist2=[]
    p=0
    while p < prime**exp:
        disc=p**2-4*n*k
        nroots=find_roots_poly([1,0,-disc],prime)
        if len(nroots)<2:
            if len(nroots)==1:
                if p not in rootlist:
                    rootlist.append(p)
                if nroots[0] not in rootlist2:
                    rootlist2.append(nroots[0])
            p+=1
            continue
        for b in nroots:
            b=lift_root2([1,0,-disc],b,prime,exp)
            if b**2%prime**exp != disc%prime**exp:
                print("fatal: "+str(roots))
                sys.exit()
            roots=find_roots_poly([1,b,-n*k],prime)
            if len(roots)>1:
                for r in roots:
                 
                    r=lift_root2([1,b,-n*k],r,prime,exp)
                    if evaluate([1,b,-n*k],r)%prime**exp == 0:
                        if p not in rootlist:
                            rootlist.append(p)
                        if b not in rootlist2:
                            rootlist2.append(b)
            elif len(roots)==1:
                der=get_derivative([1,b,-n*k])
                nb=evaluate(der,roots[0])
                r=find_roots_poly([1,nb,n*k],prime)
                for r2 in r:
                    r2=lift_root2([1,nb,n*k],r2,prime,exp)
                    if evaluate([1,nb,n*k],r2)%prime**exp == 0:
                        if p not in rootlist:
                            rootlist.append(p)
                        if b not in rootlist2:
                            rootlist2.append(b)
        p+=1    

    
    blist.append(rootlist)
    blist2.append(rootlist2)

    return blist,blist2

def build_residues(sbase,n):
    resmaps=[]
    resmaps2=[]
    for prime in sbase:
        resmaps.append([])
        resmaps2.append([])
        k=0
        while k < prime and k < 1000:
            ##Hensel here isnt doing shit. I need to have a better look here first..........
            exp=1
            temp_plist,temp_plist2=psieve_padic_solutions(prime,exp,n,k)
           # if prime == 5:
           #     print("resmaps: "+str(temp_plist)+" k: "+str(k))
            if len(temp_plist)>0:
                resmaps[-1].append(temp_plist)
                resmaps2[-1].append(temp_plist2)
            k+=1
        #if prime == 5:
           # print("resmaps: "+str(resmaps))
   # sys.exit()
    return resmaps,resmaps2

def debug_find_residues4(prime,n,exp,k):

 #   hmap={}
    blist=[prime**exp,[]]


    poly=[1,0,-4*n*k]
    polyc=copy.deepcopy(poly)
    roots=find_roots_poly(polyc,prime)
   # print("roots: "+str(roots))
    if len(roots)>0:
      #  if len(roots)==1:
     #       print("one")
        i=0
        while i < len(roots):
            roots[i]=lift_root2(poly, roots[i], prime, exp)
            if (roots[i]**2-4*n*k)%prime**exp !=0:
                print("something screwed up: "+str(roots)+" k: "+str(k))
                sys.exit(0)
            i+=1

        blist[-1].extend(roots)
    return blist

def build_disc_residues(fbase,a,n,k):
    blist_otherside=[]
    mod_otherside=1
    for prime in fbase:
        if a%prime==0:
            acpy=a
            exp=0
            while acpy%prime==0:
                acpy//=prime
                exp+=1

            #print("hit: "+str(prime))
            temp_blist=debug_find_residues4(prime,n,exp,k) ##TO DO: when prime == 2    
            if len(temp_blist)>0:
                if len(temp_blist[-1])==0:
                    print("fail: "+str(a)+" k: "+str(k)+" prime: "+str(prime)+" exp: "+str(exp))
                    sys.exit()
                mod_otherside*=prime**exp
                blist_otherside.extend(temp_blist)   

    return blist_otherside,mod_otherside

def lift_disc_residues(a,k,n,mod_otherside,blist_otherside):
    #note: this doesn't work... delete later... maybe there isnt some cycle/repeating patten here but it requires more investigation
    new_blist_otherside=[]
    new_mod_otherside=1
    i=0
    while i < len(blist_otherside): 
               
        prime=blist_otherside[i]
      #  if prime > 20:
      #      i+=2
      #      continue
        new_blist_otherside.append(prime**2)
        new_blist_otherside.append([])
        new_mod_otherside*=prime**2
        
        j=0
        while j < len(blist_otherside[i+1]):
            b=blist_otherside[i+1][j]
            while b < prime**2:
                disc=b**2-4*n*k

                if disc%prime !=0:
                    print("fatal error: "+str(prime)+" disc: "+str(disc)+" k: "+str(k))
                    sys.exit()

                disc//=prime
                roots=find_roots_poly([1,0,-disc],prime)
                if len(roots)>0:
                    new_blist_otherside[-1].append(b)
                b+=prime                
            j+=1
        i+=2
    #print(new_blist_otherside)
    return new_blist_otherside,new_mod_otherside

gain = lambda p: 2*p/(p-1) if p % 4 == 3 else 2*p/(p+1)

def best_a_mul(a, n, k, good):          # good = prefiltered odd primes: gcd(k,p)==1 and kronecker(n*k,p)==1
    #Disclaimer: This is one of the few/only functions generated with claude, as it gave a better implementation of my own which just used random sampling.
    tb = (n*k).bit_length()
    tb = tb*0.45
    best = (0.0, -1)
    for r in range(1, len(good) + 1):
        for s in itertools.combinations(good, r):
            m = math.prod(s)
            if abs((a*m*m).bit_length() - tb) < 2:
                g = math.prod(gain(p) for p in s)
                if g > best[0]:
                    best = (g, m)
    return best[1]
   
def psieve(n,ret_array,primelist_f,fbase,a_o,sbase,original_b,resmaps,resmaps2,a_mul_list):#(n,fbase,div,hmap2,ret_array):
    found=0
    a=a_o
    primes_to_mark=[]
    primes_to_mark_debug=[]
    sqr_list=[]
    i=0 
    while i < len(sbase):
        prime=sbase[i]
        sqr=find_roots_poly([1,0,-a], prime) 
        if len(sqr)>1: 
            primes_to_mark.append(i)
            sqr_list.append(sqr)
            primes_to_mark_debug.append(prime)
            if len(primes_to_mark)==9:
                break
        i+=1

    
    primes_to_check=[]
    for prime in fbase:
        if a%prime==0:
            primes_to_check.append(prime)



    k=2
    while k < 100: #To do: I know how to calculate possible "k" values for a modulus.. but there seems to be something else also going on.. kronecker(a,-n) must be 1.. but thats still not enough. Investigate later. Can add squares to a instead to optimize the interval.
        if isPrime(k,5)!=1 and k != 1:
            k+=1
            continue
       # krons=[]
        skip=0
        for prime in primes_to_check:
        #    krons.append(jacobi(n*k,prime))
            if kronecker_symbol(4*n*k,prime)==-1:
                skip=1
        if skip == 1:
            k+=1
            continue
        if jacobi(a_o,n*k) == -1: ##To do: Explore this a bit deeper eventually... does seem to hold true
            k+=1
            continue
        a_mul_ind=0
      #  if 1:
        
        
        while a_mul_ind < 1:
            sbase_temp=[]
            for prime in sbase:
                if kronecker_symbol(n*k,prime)==1 and math.gcd(prime,k)==1 and a_o%prime!=0 and kronecker_symbol(a_o, prime) == 1:
                    sbase_temp.append(prime)
            a_mul=best_a_mul(a_o,n,k,sbase_temp)
             #   a_mul=optimize_a_mul(a,n,k,sbase)
            if a_mul==-1:
                return found
            skip=0
          #  print("a_mul: "+str(a_mul))
            a=a_o*(a_mul**2)

       # print("trying k: "+str(k))
        
            if math.gcd(a,k)!=1:# or kronecker_symbol(a,n) != 1: #to do: does not need to be prime.. just need to avoid squares in the factorization... fix later
                a_mul_ind+=1
                continue


       # print("*trying k: "+str(k))
            if skip == 0:

               
                blist_otherside,mod_otherside=build_disc_residues(fbase,a,n,k)
                if mod_otherside!=a:
                    #print("skipping: "+str(a_mul))
                    a_mul_ind+=1
                    continue

                blist_disc=get_partials(mod_otherside,blist_otherside)            
                enum=[]            
                i=0
                while i < len(blist_disc):
                    enum.append(blist_disc[i+1])
                    i+=2
               # print("checking k: "+str(k)+" a: "+str(a_mul)+" "+str(jacobi(a,n*k))+" "+str(jacobi(a,k))+" "+str(jacobi(n*k,a))+" mod_otherside: "+str(mod_otherside)+" krons: "+str(krons))#+" enum: "+str(enum))
                for idx in enumerated_product(*enum):
                    b=0
                    i=0
                    while i < len(idx):
                        ind=idx[i]
                        b+=enum[i][ind]
                        i+=1
                    b_temp=b%mod_otherside

                    if (b_temp**2-4*n*k)%a != 0:
                        print("something weird went wrong: "+str(b_temp))
                        sys.exit()

                    interval=psieve_build_interval(resmaps,n,k,b_temp,a,sbase,primes_to_mark,sqr_list)
                
                ##note: Just for debugging....
                    icounter=0
                    q=0
                    while q < len(interval):
                        if interval[q]==1:
                            icounter+=1
                        q+=1
                    q=0
                    ####
                   # print("sols: "+str(icounter))

                    while q < len(interval):
                        if interval[q]==0:
                            q+=1
                            continue
                  
                        b=b_temp+mod_otherside*q
                        if b == original_b:
                            q+=1
                            continue
                        disc=b**2-4*n*k
                        if disc%a!=0:
                            print("fatal error")
                            sys.exit()
                        disc//=a
 
                        new_root=math.isqrt(abs(disc))
                        if new_root**2!=disc:
                            q+=1
                            continue


                       # new_root=math.isqrt(abs(disc))
                        if b!=original_b and b**2 not in ret_array[1]:
                            poly_val=(b)**2-4*n*k 
                            local_factors, value = factorise_fast(poly_val,primelist_f)

                            ret_array[1].append((b)**2)
                            ret_array[0].append(poly_val)
                            ret_array[2].append(local_factors)
                            ret_array[3].append([])
                           # debug_blist_otherside,debug_mod=lift_disc_residues(a,k,n,mod_otherside,blist_otherside)
                            i=0
                            while i < len(primes_to_mark):
                                ind=primes_to_mark[i]
                                prime=sbase[ind]#blist_otherside2[i][0]**blist_otherside2[i][1]
                                colist=resmaps2[ind][k%prime][1]
                                sqr=find_roots_poly([1,0,-a], prime) 
                                disc=b**2-4*n*k
                                a_inv=modinv(a,prime)
                                disc=(disc*a_inv)%prime

                                ###Important: This line below will fail for invalid solutions.. this gives a clue on how to solve what I'm trying to do here....
                                nroots=find_roots_poly([1,0,-disc], prime) 
                                
                                sqr=find_roots_poly([1,0,-a], prime) 
                                if (nroots[0]*sqr[0])%prime not in colist:
                                    print("fatal error should neer happen. Bear fail: "+str(resmaps2[ind][k%prime])+" prime: "+str(prime)+" new_root: "+str(new_root)+" k: "+str(k)+" sqr: "+str(sqr))
                                    sys.exit()
                                i+=1  
                          #  print("a*b**2+4*n*k: "+str(a*new_root**2+4*n*k)+" (a*b)**2+4*n*k*a: "+str((a*new_root)**2+4*n*k*a)+" b: "+str(b)+" a: "+str(a)+" k: "+str(k)+" mod_otherside: "+str(mod_otherside)+" new_root: "+str(new_root)+" b**2-4*n*k: "+str(b**2-4*n*k)+" b_temp: "+str(b_temp)+" primes used to mark: "+str(primes_to_mark_debug))
                            print("[i]Found one with psieve()!!!!!!!!!!!!! b: "+str(b)+" k: "+str(k)+" #smooths: "+str(len(ret_array[0]))+" index: "+str(q)+" interval[q]: "+str(interval[q])+" sols in interval: "+str(icounter)+" a_mul: "+str(a_mul)+" kronecker_symbol(a,n*k): "+str(kronecker_symbol(a,n*k))+" "+str(kronecker_symbol(a,k)))#+" interval2: "+str(interval2[q])+" k: "+str(k))
                            if kronecker_symbol(a,n*k) != 1:
                                print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!THOUSIDOADJOASDHODSA")
                                sys.exit()
                            found+=1
                            


                        q+=1
            a_mul_ind+=1
 
        if found > 0:
            return found #should be enouhg..

        k+=1


    return found

def construct_interval(ret_array,partials,n,primeslist,hmap,large_prime_bound,primeslist2):


    return

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



def lift_root(r,prime,n,quad_co,exp):
    z_inv=modinv(quad_co,prime**exp)
    c=(-quad_co*n)%prime**exp
    temp_r=r*quad_co
    zz=2*temp_r
    zz=pow(zz,-1,prime)
    x=((c-temp_r**2)//prime)%prime
    y=(x*zz)%prime
    new_r=(temp_r+y*prime)%prime**exp
    root2=(new_r*z_inv)%prime**exp
    return root2

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