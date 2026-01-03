#!python
#cython: language_level=3
#cython: profile=False
#cython: overflowcheck=False
###Author: Essbee Vanhoutte
###WORK IN PROGRESS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
###Cuda_QS_Variant
 
##References: I have borrowed many of the optimizations from here: https://stackoverflow.com/questions/79330304/optimizing-sieving-code-in-the-self-initializing-quadratic-sieve-for-pypy


import pstats, cProfile
import sympy
import random
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
import array
import numpy as np
import cupy as cp
from cupyx.scipy.ndimage import shift 
import cython
cimport cython
import os
import h5py


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
lin_sieve_size2=1_000_000
quad_sieve_size=10
g_debug=0 #0 = No debug, 1 = Debug, 2 = A lot of debug
g_lift_lim=0.5
thresvar=30  ##Log value base 2 for when to check smooths with trial factorization. Eventually when we fix all the bugs we should be able to furhter lower this.
thresvar2=10
dupe_max_prime=10_000
lp_multiplier=2
min_prime=1
g_max_diff_similar=5
g_enable_custom_factors=0
g_p=107
g_q=41
mod_mul=0.5
g_max_exp=20
quad_per_interval=1

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


        
def QS(n,factor_list,sm,xlist,flist):
    g_max_smooths=base*1+2
    if len(sm) > g_max_smooths*100: 
        del sm[g_max_smooths:]
        del xlist[g_max_smooths:]
        del flist[g_max_smooths:]  
    M2 = build_matrix(factor_list, sm, flist)
    null_space=solve_bits(M2)
    f1,f2=extract_factors(n, sm, xlist, null_space)
    if f1 != 0:
        print("[SUCCESS]Factors are: "+str(f1)+" and "+str(f2))
        return f1,f2   
    print("[FAILURE]No factors found")
    return 0,0

def extract_factors(N, relations, roots, null_space):
    n = len(relations)
    for vector in null_space:
        prod_left = 1
        prod_right = 1
        for idx in range(len(relations)):
            bit = vector & 1
            vector = vector >> 1
            if bit == 1:
                prod_left *= roots[idx]
                prod_right *= relations[idx]
            idx += 1
        sqrt_right = math.isqrt(prod_right)
        sqrt_left = prod_left# math.isqrt(prod_left)
        ###Debug shit, remove for final version
        sqr1=prod_left**2%N 
        sqr2=prod_right%N
        if sqrt_right**2 != prod_right:
            print("something fucked up1")
            time.sleep(10000)

        if sqr1 != sqr2:
            print("ERROR ERROR")
            time.sleep(10000)
        ###End debug shit#########
        sqrt_left = sqrt_left % N
        sqrt_right = sqrt_right % N
        factor_candidate = gcd(N, abs(sqrt_right+sqrt_left))
       # print(factor_candidate)
        if factor_candidate not in (1, N):
            other_factor = N // factor_candidate
            return factor_candidate, other_factor

    return 0, 0

def solve_bits(matrix):
    n=base+2
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

def build_matrix(factor_base, smooth_nums, factors):
    fb_map = {val: i for i, val in enumerate(factor_base)}

    ind=1

    M2=[0]*(base+2)
    for i in range(len(smooth_nums)):
        for fac in factors[i]:
            idx = fb_map[fac]
            M2[idx] |= ind
        ind = ind + ind


  #  offset=(small_base+2)-1
  #  ind=1
 #  for i in range(len(quad_flist)):
  #      for fac in quad_flist[i]:
   #         idx = fb_map2[fac]
   #         M2[idx+offset] |= ind
   #     ind = ind + ind
    return M2

def launch(n,primeslist1,primeslist2):
    print("[i]Total lin_sieve_size: ",lin_sieve_size)
    manager=multiprocessing.Manager()
    return_dict=manager.dict()
    jobs=[]
    start= default_timer()
    print("[i]Creating iN datastructure... this can take a while...")
    primeslist1c=copy.deepcopy(primeslist1)
    plists=[]
    i=0
    while i < build_workers:
        plists.append([])
        i+=1
    i=0
    while i < len(primeslist1):
        plists[i%build_workers].append(primeslist1c[0])
        primeslist1c.pop(0)
        i+=1
    z=0
    while z < build_workers:
        p=multiprocessing.Process(target=create_hashmap, args=(n,z,return_dict,plists[z]))
        jobs.append(p)
        p.start()
        z+=1 
    for proc in jobs:
        proc.join(timeout=0)  
    complete_hmap=[]
    while 1:
        time.sleep(1)
        check123=return_dict.values()
        if len(check123)==build_workers:
            check123.sort()
            copy1=[]
            i=0
            while i < len(check123):
                copy1.append(check123[i][1])
                i+=1
            while 1:
                found=0
                m=0
                while m < len(copy1):
                    if len(copy1[m])>0:
                        complete_hmap.append(copy1[m][0])
                        copy1[m].pop(0)
                        found=1
                    m+=1
                if found ==0:
                    break
            break
    del check123
    del return_dict
    #complete_hmap=create_hashmap(n,primeslist1)
    duration = default_timer() - start
    print("[i]Creating iN datastructure in total took: "+str(duration))

    print("[*]Launching attack with "+str(workers)+" workers\n")
    find_comb(n,complete_hmap,primeslist1,primeslist2)

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
    s=1  
    while jacobi((-s*n)%prime,prime)!=1:
        s+=1
    z_div=modinv(s,prime)  
    dist=(n*z_div)%prime
    dist=(-dist)%prime
    main_root=tonelli(dist,prime)
    if main_root**2%prime != dist:
        print("what the fuck")
    if (s*main_root**2+n)%prime !=0:
        print("fatal error123: "+str(prime)+" s: "+str(s)+" root: "+str(main_root))
        time.sleep(10000000)
    try:
        size=prime*2+1
        if size > quad_sieve_size*2:
            size= quad_sieve_size*2+1
        temp_hmap = array.array('Q',[0]*size) ##Got to make sure the allocation size doesn't overflow.... 
        temp_hmap[0]=1

        s_inv=modinv(s*z_div,prime)
        if s_inv == None or jacobi(s_inv,prime)!=1:
            print("should this ever happen?")
            return temp_hmap
        root_mult=tonelli(s_inv,prime)
        new_root=((main_root*root_mult))%prime
        if (s*new_root**2+n)%prime !=0:
            print("error2")
        new_co=(2*s*new_root)%prime
        if (new_co**2+n*4*s)%prime !=0:
            print("error")
      #  if (s*new_root**2-new_co*new_root+n)%prime !=0: ###To do: For debug delete later
           # print("error")
        if new_co > prime // 2:
            new_co=(prime-new_co)%prime  

        end=temp_hmap[0]
        temp_hmap[end]=s
        temp_hmap[end+1]=new_co
        temp_hmap[0]+=2
        s+=1   
    except Exception as e:
        print(e)
    return temp_hmap


def create_hashmap(n,procnum,return_dict,primeslist):
    i=0
    hmap=[]
    while i < len(primeslist):
        hmap_p=solve_roots(primeslist[i],n)
        hmap.append(hmap_p)
        i+=1
    return_dict[procnum]=[procnum,hmap]
    return 



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

def sumrowsby_index(a, index):
    lens=np.array([len(i) for i in index])
    id_ar = cp.zeros((len(lens), a.shape[0]),dtype=np.uint16)
    c=np.concatenate(index)
    r=np.repeat(np.arange(len(index)),lens)
    id_ar[r,c]=1
    return id_ar.dot(a)[0]

cdef process_interval2d(n,ret_array,quad_can,primelist_f,large_prime_bound,partials,root_list_complete,cmod,factor_ranking,fb_map,bSeenOnly):#,lin,cmod,sum_list):
    linsize=lin_sieve_size
    if bSeenOnly==1:
        linsize=lin_sieve_size2
        threshold = int(math.log2((linsize)*math.sqrt(abs(n))) - thresvar2)
    else:
        threshold = int(math.log2((linsize)*math.sqrt(abs(n))) - thresvar)
    #threshold = int(math.log2(abs(n)) - thresvar)
    h5f = h5py.File('res.hdf5','r')
    found=0
    linlen=len(root_list_complete)
    u=0
    cdef Py_ssize_t k


    u2=0
    interval2d=cp.asarray(h5f["intervals/"+str(cmod)+"/"+str(quad_can)][:])
    while u2 < linlen:
        interval=interval2d[u2]
            
        cp.putmask(interval, interval<threshold, 0)
        indexlist=cp.nonzero(interval)
        indexlist=cp.asnumpy(indexlist[0])
        temp=cp.asnumpy(interval)
        root=root_list_complete[u2]
           # print("Checking lin: "+str(lin)+" quad: "+str(quad_can)+" cmod: "+str(cmod)+" u2: "+str(u2)+" u: "+str(u)+" temp: "+str(temp))

        k=0
        length=len(indexlist)
        while k < length:# length:  
            i=indexlist[k]
            if temp[i]>threshold:       
                x=root+cmod*int(i)
                poly_val=quad_can*x**2+n

                if poly_val%cmod !=0:
                    print("ERGH")
                    time.sleep(1000)                    
                new_root=quad_can*x
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
                if value != 1:
                    if value < large_prime_bound:
                        if value in partials:
                            rel, lf, pv = partials[value]
                            if rel == new_root:
                                k+=1
                                continue
                            new_root *= rel
                            local_factors ^= lf
                            poly_val *= pv
                        else:
                            partials[value] = (new_root, local_factors, poly_val)
                            k+=1
                            continue
                    else:
                        k+=1 
                        continue         
                if new_root not in ret_array[1]:
                    factor_ranking.append([])
                    local_factors=list(local_factors)
                    local_factors.sort()
    
                    for fac in local_factors:
                        if fac != -1 and fac !=2 and bSeenOnly==0:
                            idx=fb_map[fac]
                            factor_ranking[-1].append(idx)
                    if g_debug == 1:
                        if bSeenOnly==0:
                            print("***seen_primes: "+str(local_factors)+" cmod: "+str(cmod))
                        if bSeenOnly==1:
                            print("seen_primes: "+str(local_factors)+" cmod: "+str(cmod))
                    found+=1
                    ret_array[1].append(new_root)
                    ret_array[0].append(poly_val)
                    ret_array[2].append(local_factors)
            k+=1
        u2+=1
  
    h5f.close()
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
        factor=factor_base[i]
        if value % factor == 0:
            factors ^= {factor}
            value //= factor
            if value % factor == 0:
                return -1, -1
        i+=1
    return factors, value


def filter_quads(qbase,n):
    ###Note: We look for quadratic coefficients that factor over the factor base but have no even exponents. This garantuees unique results
    valid_quads=[]
    valid_quads_factors=[]

    roots=[]
    i=1
    while i < quad_sieve_size+1:
        quad_local_factors, quad_value = factorise_fast_quads(i,qbase) 
        if quad_value != 1:
            i+=1
            continue
        valid_quads.append(i)
        valid_quads_factors.append(quad_local_factors)
        i+=1
    return valid_quads,valid_quads_factors


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

def get_lin(hmap,cfact,local_mod,indexes,quad_co,n,roots2d):
    all_lin_parts=[]
    j=0
    lin=0
    #print("indexes: ",indexes)
    while j < len(indexes):
        ind=indexes[j]
        prime=math.isqrt(cfact[j])
 

        r1=int(roots2d[quad_co,ind])#hmap[ind][2]
    #    print("prime: "+str(prime)+" quad_co: "+str(quad_co))
        r1=lift_root(r1,prime,n,quad_co,2)
       # r1=quad_co*r1*2
       # r1=lift_b(prime,n,r1,quad_co,2)
      #  if (r1**2+n*4*quad_co)%prime**2 !=0:
           # print("fatal error")
          #  print("prime: "+str(prime)+" index: "+str(ind)+" hmap[ind][2]: "+str(hmap[ind][1]))
          #  time.sleep(1000)
        prime=prime**2
        aq = local_mod // prime
        invaq = modinv(aq%prime, prime)
        gamma = r1 * invaq % prime
        lin+=aq*gamma
        all_lin_parts.append(aq*gamma)
        j+=1
    lin%=local_mod
    return lin,all_lin_parts


def get_lin2(hmap,cfact,local_mod,indexes,quad_co,n,roots2d):
    all_lin_parts=[]
    j=0
    lin=0
    #print("indexes: ",indexes)
    while j < len(indexes):
        ind=indexes[j]
        prime=cfact[j]
 

        r1=int(roots2d[quad_co,ind])#hmap[ind][2]
    #    print("prime: "+str(prime)+" quad_co: "+str(quad_co))
      #  r1=lift_root(r1,prime,n,quad_co,2)
       # r1=quad_co*r1*2
       # r1=lift_b(prime,n,r1,quad_co,2)
      #  if (r1**2+n*4*quad_co)%prime**2 !=0:
           # print("fatal error")
          #  print("prime: "+str(prime)+" index: "+str(ind)+" hmap[ind][2]: "+str(hmap[ind][1]))
          #  time.sleep(1000)

        aq = local_mod // prime
        invaq = modinv(aq%prime, prime)
        gamma = r1 * invaq % prime
        lin+=aq*gamma
        all_lin_parts.append(aq*gamma)
        j+=1
    lin%=local_mod
    return lin,all_lin_parts


@cython.boundscheck(False)
@cython.wraparound(False)
cdef build_database2interval(long long [:] primeslist,quad,hmap,n,root_list_complete,cmod,roots2d,bSeenOnly,factor_ranking):
    
    h5f = h5py.File('res.hdf5','a')
    linlen=len(root_list_complete) 
    cdef Py_ssize_t i
    cdef Py_ssize_t v
    if bSeenOnly==1:
        linsize=lin_sieve_size2
    else:
        linsize=lin_sieve_size
    interval_single=np.zeros([linlen,linsize],dtype=np.uint16)    

    h5f.create_dataset("intervals/"+str(cmod)+"/"+str(quad), data=interval_single,chunks=True)

        
    interval_single = cp.asarray(h5f["intervals/"+str(cmod)+"/"+str(quad)][:])
    i=0
    while i < len(primeslist):
        prime=primeslist[i]
        log=round(math.log2(prime))
        xc = int(roots2d[quad,i])
        if xc ==0:
            i+=1
            continue
        v=0
        while v < len(root_list_complete):
            root=root_list_complete[v]
            if cmod%prime == 0 or quad%prime ==0:
                v+=1
                continue
            lcmod=cmod%prime
            modi=modinv(lcmod,prime)
            x=xc
            exp=1
            while exp < g_max_exp+1:
                
                p=prime**exp
                if exp > 1:
                    x=lift_root(x,prime**(exp-1),n,quad,exp)

                root_dist1=solve_lin_con(cmod,x-root,p)
                x_b=(p-x)%p 
                root_dist2=solve_lin_con(cmod,x_b-root,p)
                if bSeenOnly==0 or (bSeenOnly ==1 and prime < dupe_max_prime):
                    hit=0
                    if root_dist1 < linsize+1:
                        interval_single[v,root_dist1::p]+=log
                        hit=1
                    if root_dist2 < linsize+1:
                        if root_dist1 != root_dist2:
                            interval_single[v,root_dist2::p]+=log  
                            hit=1
                    if g_debug ==1:
                        CAN=quad*(root+root_dist1*cmod)**2+n
                        if CAN % (cmod*p)  != 0:

                            print("EROROREROR",prime)
                            sys.exit()
                    if hit ==0:
                        break
                elif bSeenOnly ==1 and prime > dupe_max_prime-1 and exp%2 ==0:
                    hit=0
                    if root_dist1 < linsize+1:
                        interval_single[v,root_dist1::p]+=(log*2)
                        hit=1
                    if root_dist2 < linsize+1:
                        if root_dist1 != root_dist2:
                            interval_single[v,root_dist2::p]+=(log*2)
                            hit=1
                    if g_debug ==1:
                        CAN=quad*(root+root_dist1*cmod)**2+n
                        if CAN % (cmod*p)  != 0:

                            print("EROROREROR",prime)
                            sys.exit()
                    if hit ==0:
                        break

                exp+=1
            v+=1    
        i+=1
    h5f["intervals/"+str(cmod)+"/"+str(quad)][...]=cp.asnumpy(interval_single)
    del interval_single   


    h5f.close()
    return 

def build_2drootmap(primeslist,hmap,n):
    roots2d=np.zeros([quad_sieve_size+1,base],dtype=np.int32)
    i=0
    while i < len(primeslist):
        prime=primeslist[i]
        j=0
        while j < prime and j < quad_sieve_size+1: ###To do: Only go up to prime
            quad=j#+offset
            #print("Building quad: "+str(quad))
            z=hmap[i][1]
            z_div=modinv(z,prime)
            z_inv=modinv(quad*z_div,prime)
            if z_inv == None or jacobi(z_inv,prime)!=1:
                j+=1
                continue  
            
            y=hmap[i][2]
            root_mult=tonelli(z_inv,prime)
            x=get_root(prime,y,z)
            x=(x*root_mult)%prime
            if bitlen(x)>32:
                print("big root, increase dtype in build_2drootmap()")
            roots2d[quad::prime,i]=x

            if (quad*x**2+n)%prime !=0:
                print("fatal error, sad face :(")
                sys.exit()
            j+=1
        i+=1
   # print(roots2d)
    return roots2d

def gen_modulus_and_interval(quad,n,primeslist,seen,tnum,close_range,too_close,LOWER_BOUND_SIQS,UPPER_BOUND_SIQS,hmap,primeslist_a,grays,roots2d,bSeenOnly,factor_ranking):
    new_mod,cfact,indexes=generate_modulus(n,primeslist,seen,tnum,close_range,too_close,LOWER_BOUND_SIQS,UPPER_BOUND_SIQS,bitlen(tnum),hmap,quad)
    if new_mod ==0:
        return 0,0
    root_list_complete=[]
    lin,lin_parts=get_lin(hmap,cfact,new_mod,indexes,quad,n,roots2d)
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
        if (z*lin**2+n)%new_mod !=0:
            print("super big error")
            time.sleep(1000)
        lin_co_array.append(lin)
        poly_ind+=1

    
    build_database2interval(primeslist_a,quad,hmap,n,lin_co_array,new_mod,roots2d,bSeenOnly,factor_ranking)
    return lin_co_array,new_mod


def construct_interval(ret_array,partials,n,primeslist,hmap,large_prime_bound,primeslist2):
    grays = get_gray_code(20)
    cp.set_printoptions(
    linewidth=500,   # Max characters per line before wrapping
    precision=1,     # Decimal places for floats
    suppress=True,   # Avoid scientific notation for small numbers
    edgeitems=20,#np.inf # Show full array without truncation
    threshold=10
)    
    found=0
    file_path = "res.hdf5"
    if os.path.isfile(file_path):
        os.remove(file_path)

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


    close_range=10
    too_close=5
    LOWER_BOUND_SIQS=400
    UPPER_BOUND_SIQS=4000
    tnum=int(((n)**0.5) /(lin_sieve_size))
    seen=[]
    i=0
    while i < len(primeslist):
        prime=primeslist[i]
        h5f = h5py.File('res.hdf5','a')
        dat=np.zeros(prime,dtype=np.int64)
        h5f.create_dataset("tiles/"+str(prime), data=dat,chunks=True)  
        h5f.close()
        i+=1
    print("[i]Building 2d root map")
    roots2d=build_2drootmap(primeslist_a,hmap,n)
    print("[i]Entering attack loop")
    valid_quads,valid_quads_factors=filter_quads(qlist,n)
    fb_map = {val: i for i, val in enumerate(primeslist)}
   # print("fb_map: ",fb_map)
   # factor_ranking=[]#np.zeros(len(primeslist),dtype=np.uint16)
    while 1:
        factor_ranking=[]
        quad=1
        root_list_complete,new_mod=gen_modulus_and_interval(quad,n,primeslist,seen,tnum,close_range,too_close,LOWER_BOUND_SIQS,UPPER_BOUND_SIQS,hmap,primeslist_a,grays,roots2d,0,factor_ranking)
        if new_mod ==0:
            continue            
        found+=process_interval2d(n,ret_array,quad,primelist_f,large_prime_bound,partials,root_list_complete,new_mod,factor_ranking,fb_map,0)#,lin,new_mod,sum_list)
      #  print("", end=f"[i]Smooths: {len(ret_array[2])}\r")
        i=len(factor_ranking)-1
        while i >-1:
            l=0
            while l < len(valid_quads):
                #print("v: "+str(v))
                new_quad=valid_quads[l]#primeslist[factor_ranking[i][-1]]
                target=int(((n/new_quad)**0.5) /(lin_sieve_size2))
                # new_quad=primeslist[i]
            
                skip=0
                k=len(factor_ranking[i])-1
                new_cfact=[]
                new_indexes=[]
                lmod=1
                while k >-1:
                    if jacobi((-new_quad*n)%primeslist[factor_ranking[i][k]],primeslist[factor_ranking[i][k]])==1:
                   # skip=1
                   # break
                        new_cfact.append(primeslist[factor_ranking[i][k]])
                        new_indexes.append(factor_ranking[i][k])
                        lmod*=primeslist[factor_ranking[i][k]]
                        if bitlen(lmod)>(keysize//2):
                            break
                        if abs(bitlen(lmod)-bitlen(target))<3:
                            break
                #    elif new_quad == 1:
                    #    print("WHY????????????????????????????????????????????????????",jacobi((-new_quad*n)%primeslist[factor_ranking[i][k]],primeslist[factor_ranking[i][k]]))
                     #   print(primeslist[factor_ranking[i][k]])
                    k-=1
            
                if abs(bitlen(lmod)-bitlen(target))>2:
                    l+=1
                    continue
                if g_debug ==1:
                    print("Modulus factors: "+str(new_cfact))
         #   print("trying quad: "+str(new_quad)+" mod: "+str(lmod)+" cfact_new: "+str(new_cfact))
                lin,lin_parts=get_lin2(hmap,new_cfact,lmod,new_indexes,new_quad,n,roots2d)
                z=new_quad#quadlist[j]
                lin_co_array=[]
                q=0
           
                lin2=lin#lin3%new_mod
                poly_ind=0
                end = 1 << (len(new_cfact) - 1)
                lin=0
                while poly_ind < end:
                    if poly_ind != 0:
                        v,e=grays[poly_ind]
                        lin=(lin + 2 * e * lin_parts[v])%lmod
                    else:
                        lin=lin2
                    if (z*lin**2+n)%lmod !=0:
                        print("super big error")
                        time.sleep(1000)
                    lin_co_array.append(lin)
                    poly_ind+=1
                if len(lin_co_array) > 10000:
                    print("fatal error: "+str(new_cfact)+" lin_co_array: "+str(len(lin_co_array))+" cmod bits: "+str(bitlen(lmod)))#+" linlen: "+str(len))
                    sys.exit()
                build_database2interval(primeslist_a,new_quad,hmap,n,lin_co_array,lmod,roots2d,1,factor_ranking)
                found+=process_interval2d(n,ret_array,new_quad,primelist_f,large_prime_bound,partials,lin_co_array,lmod,factor_ranking,fb_map,1)#,lin,new_mod,sum_list)
                print("", end=f"[i]Smooths: {len(ret_array[2])}\r")
                h5f = h5py.File('res.hdf5','a')
                del h5f["intervals/"+str(lmod)]
                h5f.close()
                    #   print("", end=f"[i]Smooths: {len(ret_array[2])}\r")
                if found > 500 or len(ret_array[0]) > base+2:
                    test,test2=QS(n,primelist,ret_array[0],ret_array[1],ret_array[2]) 
                    found=0 
                    if test !=0:
                        return
               
                l+=1
            i-=1 
        ######To do: Return a list of all odd exponent factors for any smooths found. Now iterate that and call gen_modulus_and_interal for those quadratic coefficients.
        if found > 500 or len(ret_array[0]) > base+2:
            test,test2=QS(n,primelist,ret_array[0],ret_array[1],ret_array[2]) 
            found=0 
            if test !=0:
                print("\n\n\n\nFound at: ",len(ret_array[0]))
                return
        h5f = h5py.File('res.hdf5','a')
        del h5f["intervals/"+str(new_mod)]
        h5f.close()
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

def find_comb(n,hmap,primeslist1,primeslist2):
    #To do: Pointless function
    ret_array=[[],[],[],[]]
    partials={}
    large_prime_bound = primeslist1[-1] ** lp_multiplier
    construct_interval(ret_array,partials,n,primeslist1,hmap,large_prime_bound,primeslist2)
    return 0

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
    i=0
    while len(primeslist2) < base:
        if n%primeslist[i] !=0:
            primeslist2.append(primeslist[i])
        i+=1
    launch(n,primeslist1,primeslist2)     
    duration = default_timer() - start
    print("\nFactorization in total took: "+str(duration))

