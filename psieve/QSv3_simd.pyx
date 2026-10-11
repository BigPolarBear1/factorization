#!python
#cython: language_level=3
#cython: profile=False
#cython: overflowcheck=False
###Author: Essbee Vanhoutte
###WORK IN PROGRESS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
###Cuda_QS_Variant
###This code is pro-LGBTQ.

##References: I have borrowed many of the optimizations from here: https://stackoverflow.com/questions/79330304/optimizing-sieving-code-in-the-self-initializing-quadratic-sieve-for-pypy

#To do: If we're going to restrict ourselves to quadratics.. then tonelli is better for root finding
#To do: Lots of overly complex code that can be simplified in the SIQS variant, especially since we restrict k to 1.


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
import almostsq
from quadroots import QuadRoots
from psievecore import PSieveCore, KOrder


k_max=10_000
psieve_trigger=20           # try psieve after this many new SIQS relations
psieve_window=8192          # positions sieved per root class, starting at the crossing X = sqrt(4nK)
psieve_sieve_primes=24      # small primes filtering each window (a mod-64 condition is applied on top)
psieve_share=1.0            # psieve may use this much time relative to the SIQS time since the last attempt ...
psieve_min_budget=1.0       # ... but at least this many seconds per attempt
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
    fb_map = {val: i for i, val in enumerate(primeslist)}
    psieve_cores={}
    siqs_clock=default_timer()
    hmap=create_hashmap(n,primeslist)
    ret_array=[[],[],[],[]]
    partials={}
    large_prime_bound = primeslist[-1] ** lp_multiplier
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
    qr=QuadRoots(primeslist,hmap,100_000,k_max)
    print("[i]Building 2d root map")
    roots2d=build_2drootmap(primeslist_a,hmap,n,qr)
    print("[i]Entering attack loop")
    siqs_clock=default_timer()
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
            found+=process_interval2d(n,ret_array,quad,primelist_f,large_prime_bound,partials,lin,new_mod,factor_ranking,fb_map,0,seen_factors,interval,primeslist,[],[],valid_quads,valid_quads_factors,qlist,primelist,[],[])#,lin,new_mod,sum_list)
            if found > psieve_trigger or len(ret_array[0]) > base+10:
                print("[i]#B-smooths: "+str(len(ret_array[0])))
                almostsq_result,leftovers=generate_almostsq(ret_array[2],primeslist,ret_array[0])
             #   almostsq_result=generate_almostsq_smallprimes(ret_array[2],primeslist,ret_array[0])
                # psieve: each attempt goes to the usable leftover whose K range is least explored (small K give the most hits)
                budget=max(psieve_min_budget,psieve_share*(default_timer()-siqs_clock))
                pstart=default_timer()
                usable=[]
                for res in almostsq_result:
                    div=res[0]
                    if div!=1 and div%2!=0 and bitlen(div)<(keysize*0.5) and abs(div)<(1<<62) and 2 not in res[3] and len(usable)<400:
                        usable.append(res)
                while len(usable)>0 and default_timer()-pstart < budget:
                    pick=min(usable,key=lambda r:(psieve_cores[r[0]]['k'] if r[0] in psieve_cores else 0, leftover_L(r[0])))   # least explored first, then smallest L(1,chi)
                    div=pick[0]
                    X0=1
                    for j in pick[2]:
                        X0=(X0*math.isqrt(ret_array[1][j]))%n
                    print("[i]psieve2 on a: "+str(div)+" (bits: "+str(bitlen(div))+", primes: "+str(len([p for p in pick[3] if p!=-1]))+") multipliers tried so far: "+str(psieve_cores[div]['k'] if div in psieve_cores else 0)+" | leftovers tried so far: "+str(len(psieve_cores)))
                    factors=psieve2(n,ret_array,primelist_f,div,pick[3],pick[1],X0,budget-(default_timer()-pstart),psieve_cores)
                    if factors is not None:
                        print("[SUCCESS]Factors are: "+str(factors[0])+" and "+str(factors[1]))
                        sys.exit()
                    if div not in psieve_cores:                 # psieve2 refused it
                        usable.remove(pick)
                        continue
                    break
                siqs_clock=default_timer()
           #     if g_debug ==1:
           #         print("seen_factors: ",seen_factors)
           #     print("[i]Performing linear algebra")
                if len(ret_array[0]) > base+10:
                    print("[*]Trying linear algebra")
                    test,test2=QS(n,primelist,ret_array[0],ret_array[2],ret_array[1],ret_array[3])
                    
      #      test,test2=QS(n,primelist,ret_array[0],ret_array[1],ret_array[2]) 
                found=0 


            poly_ind+=1
      
     
        ######To do: Return a list of all odd exponent factors for any smooths found. Now iterate that and call gen_modulus_and_interal for those quadratic coefficients.
    

    return 

def solve_roots(prime,n):
    size=min(k_max,prime)
    sol_array=array.array("i",[0]*2)
    k=0
    while k < size:
        if jacobi(k*n,prime)==1:
            roots=find_roots_poly([1,0,-n*k],prime)
            if len(roots)==0:
                print("fatal ")
                sys.exit()
            if len(roots)>1:
                if roots[0]<roots[1]:
                    sol_array[0]=k
                    sol_array[1]= roots[0]
                    break
                else:
                    sol_array[0]=k
                    sol_array[1]= roots[1]
                    break
            else:
                    sol_array[0]=k
                    sol_array[1]= roots[0]
                    break
        k+=1
    return sol_array

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

def solve_lin_con(a,b,m):
    ##ax=b mod m
    g=gcd(a,m)
    a,b,m = a//g,b//g,m//g
    return pow(a,-1,m)*b%m  
    
#@cython.boundscheck(False)
#@cython.wraparound(False)
#cdef factorise_fast_debug(value,long long [::1] factor_base,max_exp):
#    factors = set()
#    seen_primes=[]
#    if value < 0:
#        factors ^= {-1}
#        value = -value
#    while value % 2 == 0:
#        factors ^= {2}
#        value //= 2
    #    seen_primes.append(2)
#    length=factor_base[0]#len(factor_base)#factor_base[0]
#    cdef Py_ssize_t i=1
#    while i < length:
#        factor=factor_base[i]
#        exp=1
#        while value % factor == 0:
#            factors ^= {factor}
#            value //= factor
#            if exp < (max_exp+1):
#                seen_primes.append(factor)
##            exp+=1
#        i+=1
#    return factors, value,seen_primes



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

def generate_almostsq(facs, fb, values, keep=20):
    #Disclaimer: This functino was generated with claude
    """All single relations plus the lightest combinations, smallest leftover first.

    Returns a list of (leftover, root, idxs) with
        prod(values[i] for i in idxs) == leftover * root**2
    """
    cols = sorted({p for f in facs for p in f})       # -1 sorts first, an ordinary column
    col = {p: j for j, p in enumerate(cols)}

    masks = []
    for f in facs:
        mask = 0
        for p in f:
            mask ^= 1 << col[p]
        masks.append(mask)

    # every relation on its own, then the combinations from the elimination
    candidates = [(mask, 1 << i) for i, mask in enumerate(masks)]
    candidates += almostsq.from_masks(masks, len(cols), keep=keep)

    seen = set()
    results = []
    for prime_mask, history_mask in candidates:
        if history_mask in seen:                      # same set of relations already listed
            continue
        seen.add(history_mask)
        odd  = [cols[j] for j in almostsq.bits(prime_mask)]
        idxs = almostsq.bits(history_mask)
        product = math.prod(values[i] for i in idxs)
        leftover = math.prod(odd)                     # product of the odd primes, sign included
        square, rem = divmod(product, leftover)
        root = math.isqrt(square)
        assert rem == 0 and root * root == square
        results.append((leftover, root, idxs, odd))

    results.sort(key=lambda r: (sum(1 for p in r[3] if p != -1), abs(r[0])))

    #for res in results:
   ##     print(str(bitlen(res[0]))+" amount: "+str(len(res[2]))+" total len: "+str(len(values))+" primes: "+str(res[3]))
   # print(results)
    leftovers=sorted({r[0] for r in results if abs(r[0])!=1},key=abs)   # signed, smallest first
    return results,leftovers

def generate_almostsq_smallprimes(facs, fb, values, limit=50):
    """Single relations and combinations, ranked by the largest prime left in the
    leftover (smaller is better), then by leftover size.

    Returns up to `limit` tuples (leftover, root, idxs, odd) with
        prod(values[i] for i in idxs) == leftover * root**2
    idxs index into the lists as passed in.
    """
    # A repeated value would give the trivial square v*v and rank first: use each value once.
    first = {}
    for i, v in enumerate(values):
        first.setdefault(v, i)
    use = sorted(first.values())

    # Largest prime = column 0.  The elimination pivots on the lowest column, so every
    # reduced row's pivot is its largest prime.
    cols = sorted({p for i in use for p in facs[i]}, reverse=True)   # -1 ends up last
    col = {p: j for j, p in enumerate(cols)}
    masks = []
    for i in use:
        mask = 0
        for p in facs[i]:
            mask ^= 1 << col[p]
        masks.append(mask)

    candidates = [(mask, 1 << k) for k, mask in enumerate(masks)]     # each relation alone
    candidates += almostsq.reduce_with_history(masks, len(cols))      # reduced rows

    seen = set()
    ranked = []
    for prime_mask, history_mask in candidates:
        if history_mask in seen:
            continue
        seen.add(history_mask)
        odd = sorted(cols[j] for j in almostsq.bits(prime_mask))
        leftover = math.prod(odd)
        ranked.append((max((p for p in odd if p > 0), default=1), abs(leftover), leftover, odd, history_mask))
    ranked.sort(key=lambda r: r[:2])

    results = []
    for _, _, leftover, odd, history_mask in ranked[:limit]:          # big products only for what is returned
        idxs = [use[k] for k in almostsq.bits(history_mask)]
        square, rem = divmod(math.prod(values[i] for i in idxs), leftover)
        root = math.isqrt(square) if square >= 0 else 0
        if rem or root * root != square:
            raise RuntimeError("self-check failed for relations %r" % (idxs,))
        results.append((leftover, root, idxs, odd))
    return results

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
                poly_val=(2*new_root)**2-4*n*quad_can
            else: 
                poly_val=(2*new_root)**2+4*n*quad_can
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
            if value != 1:# and isPrime(value,5)==0:#math.isqrt(abs(value))**2 != value:
                k+=1
                continue    
            if 2*new_root not in ret_array[1]:
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
                ret_array[1].append((2*new_root)**2)
                ret_array[0].append(poly_val)
                ret_array[2].append(local_factors)
                ret_array[3].append([])
                div_fac=[]
                                
                if len(ret_array[0])>(base+10):
                    return found
        k+=1

  

    if g_debug ==1:
        if bSeenOnly == 1:
            print("found: "+str(found)+" total: "+str(len(ret_array[0]))+" cmod bits: "+str(bitlen(cmod))+" quad: "+str(quad_can))
        else:
            print("***found: "+str(found)+" total: "+str(len(ret_array[0]))+" cmod bits: "+str(bitlen(cmod))+" quad: "+str(quad_can))
    return found

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

def build_2drootmap(primeslist,hmap,n,qr):
    roots2d=np.zeros([quad_sieve_size+1,base],dtype=np.int32)
 
    pr=np.array(primeslist,dtype=np.int64)
    n_mod=np.array([n%p for p in primeslist],dtype=np.int64)
    if quad_sign != "neg":
        n_mod=(-n_mod)%pr
    qi=0
    while qi < len(qr.quads) and qi < 1:
        quad=qr.quads[qi]
        if quad != 1:
            print("fail")
            sys.exit()
        qr.roots(qi,roots2d[quad])          # fills the whole row for this quad
        x=roots2d[quad,:len(primeslist)].astype(np.int64)
        if np.any((x>=0)&((x*x%pr*quad-n_mod)%pr!=0)):
            print("fatal error, sad face :(")
            sys.exit()
        qi+=1
    roots2d[roots2d<0]=0                    # skipped entries back to 0, as before
    return roots2d

#def find_r(mod,total):
#    mo,i=mod,0
#    while (total%mod)==0:
#        mod=mod*mo
#        i+=1
#    return i

#def brute_force_padic_solutions(prime,k,n,a,sqr):
    ##Add proper hensel later and use this to verify
#    solutions=[]
#    blist=[]
#    exp=1
#    b=0
#    while b < prime**exp:
#       # print("b: "+str(b)+" prime: "+str(prime))
#        poly=[1,-b*sqr,n*k]
#        #poly=[1,-b,n*k]
#        roots=[]
#        x=0
#        while x < prime**exp:
#            if evaluate(poly,x)%prime**exp ==0:
#                roots.append(x)
#                #disc=a*(b**2)-4*n*k*a
#                #disc//=a
#
#            x+=1
#        if len(roots)>0:
#            solutions.extend([b,roots,0])
#        b+=1
#  #  print("solutions: "+str(solutions))
#    if prime == 2:
#        max_lift=12
#    elif prime == 3:
#        max_lift=8
#    elif prime == 5:
#        max_lift=5
#    elif prime == 7:
#        max_lift=3
#    elif prime < 20:
#        max_lift=2
#    else:
#        max_lift=1
#    max_lift=1
#    if max_lift > 1:
#        exp+=1
#        while exp < max_lift:
#            sqr_lifted=lift_root2([1,0,-a],sqr,prime,exp)
#            if evaluate([1,0,-a],sqr_lifted)%prime**exp !=0:
#                print("fatal error")
#                sys.exit()
#            am_to_lift=0
#            new_solutions=[]
#            i=0
#            while i < len(solutions):
#                b=solutions[i]
#                while b < prime**exp:
#                #print("b: "+str(b)+" prime: "+str(prime)+" exp: "+str(exp)+" prime**exp: "+str(prime**exp))
#                    poly=[1,-b*sqr_lifted,n*k]
#                    der=get_derivative(poly)
#                #poly=[1,-b,n*k]
#                    new_roots=[]
#                    hit=0
#                    if solutions[i+2]==0:
#                        roots=solutions[i+1]
#                        j=0
#                        while j < len(roots):
#                            r=roots[j]
#                 #   print("r: "+str(r)+" b: "+str(b))
#                            while r < prime**exp:
#                                if evaluate(poly,r)%prime**exp ==0:
#                                    new_roots.append(r)
#                                    deriv=evaluate(der,r)
#                                    r0=find_r(prime,deriv)
#                                    ceiling=(prime*r0)
#                                    ceiling=prime**ceiling
#                                    if prime**exp > ceiling:
#                                        hit=1
#                              #  print("******************prime**exp: "+str(prime**exp)+" b: "+str(b)+" r: "+str(r)+" "+str(deriv)+" ceiling: "+str(ceiling))
#                           # else:
##                             #   print("prime**exp: "+str(prime**exp)+" b: "+str(b)+" r: "+str(r)+" "+str(deriv)+" ceiling: "+str(ceiling))
 #                               r+=prime**(exp-1)
 #                           j+=1
 #                       if len(new_roots)>0:
 #                   #print("b: "+str(b)+" "+str(len(new_roots)))
 #                           if hit==1:
 #                               new_solutions.extend([b,exp,1])
 #                           else:
 #                               new_solutions.extend([b,new_roots,0])
 #                               am_to_lift+=1
#
#                    else:
#                        new_solutions.extend([solutions[i],solutions[i+1],solutions[i+2]])
#                        break
#                #if solutions[i+2]==1 and len(roots)==0:
#                #    print("fatal error")
#                #    sys.exit()
#                #if solutions[i+2]==1 and hit==0:
#                #    print("fata; error2")
#                #    sys.exit()
#
#
#
#                    b+=prime**(exp-1)
#                i+=3
#            solutions=new_solutions
#       # print("solutions: "+str(prime**exp)+" "+str(len(solutions)//3))
#      #  if am_to_lift==0:
#          #  print("!!!!!!!!!!!!!!!!FULLY LIFTED AT: "+str(exp))
#          #  break
#            if exp+1 == max_lift:
#                break
#            exp+=1
#
#    
#    blist.append(prime**(exp))
#    blist.append([])
#    i=0
#    while i < len(solutions):
#        if solutions[i+2]==1:
#            texp=solutions[i+1]
#            tempb=solutions[i]
#            while tempb < prime**(exp):
#                blist[-1].append(tempb)
#                ##To do: can add a verification loop here later if I'm bored to make sure everything added is correct
#
#                tempb+=prime**texp
#        else:    
#            blist[-1].append(solutions[i])
#
#        i+=3
#    blist[-1].sort()
#    #print("solutions len: "+str(len(blist[-1]))+" k: "+str(k))
#
#    return blist

g_small_primes=None
g_leftover_L={}
def leftover_L(a):
    # L(1,chi) of Q(sqrt(a)) from the primes below 3000. The expected number of psieve hits for a multiplier K is
    #   2^(primes of K) * (number of genera) * window / (sqrt(discriminant) * L(1,chi)),
    # which is the class number formula read backwards; measured hit counts follow it, so a small L is better.
    global g_small_primes
    if g_small_primes is None:
        g_small_primes=get_primes(3,3000)
    if a not in g_leftover_L:
        v=1.0/(1-(1 if a%8 in (1,7) else -1)/2) if a%4==1 else 1.0
        for f in g_small_primes:
            if a%f!=0:
                v/=(1-jacobi(a%f,f)/f)
        g_leftover_L[a]=v
    return g_leftover_L[a]

g_psieve_pool=None
g_psieve_wcores={}
def _psieve_task(args):
    # runs in a worker process: one slice of K for one leftover
    n,a,primes,nsieve,ks,window=args
    key=(n,a)
    if key not in g_psieve_wcores:
        if len(g_psieve_wcores)>8:
            g_psieve_wcores.clear()
        g_psieve_wcores[key]=PSieveCore(n,a,primes,nsieve)
    return g_psieve_wcores[key].run(0,0,window,1<<30,0,0,ks)

def psieve2(n,ret_array,primelist_f,a,odd,root0,X0,budget,state):
    # Looks for X^2 - a*y^2 = 4nK over all multipliers K, where a is the leftover of an SIQS relation (or combination):
    # prod(values) = a*root0^2 and X0^2 = prod(values) mod n. A hit has exactly the odd factors of that relation, so the
    # pair is a dependency and is tried with a gcd straight away. K is not taken in ascending order: KOrder hands out
    # squarefree K whose prime factors all split, best value first (value = 2^(number of prime factors) / K^(1/4)),
    # because such a K carries 2^w ideals and so 2^w chances of a solution for the same sieving cost.
    # Always runs at least one slice. Returns (f1,f2) or None. Progress through K is kept in state[a] between calls.
    global g_psieve_pool
    primes=[p for p in odd if p!=-1]
    if len(primes)==0 or len(primes)>12 or 2 in primes or math.prod(primes)!=abs(a) or abs(a)>=(1<<62) or max(primes)>50_000_000 or (X0*X0-a*root0*root0)%n!=0:
        return None
    if a not in state:
        state[a]={'core':PSieveCore(n,a,primes,psieve_sieve_primes),'order':KOrder(a),'seen':set(),'same':0,'k':0}
    st=state[a]
    start=default_timer()
    first=True
    while first or default_timer()-start < budget:
        first=False
        if workers>1:
            if g_psieve_pool is None:
                g_psieve_pool=multiprocessing.Pool(workers)
                import atexit
                atexit.register(g_psieve_pool.terminate)   # so the worker processes are gone before the interpreter shuts down
            parts=g_psieve_pool.map(_psieve_task,[(n,a,primes,psieve_sieve_primes,st['order'].take(1<<17),psieve_window) for w in range(workers)])
            hits=[h for part in parts for h in part]
        else:
            hits=st['core'].run(0,0,psieve_window,1<<30,0,0,st['order'].take(1<<15))
        st['k']=st['order'].used
        for K,X,y in hits:
            if math.gcd(X,y)>2 or X in st['seen']:           # scaled copy of a hit at a smaller K, or seen before
                continue
            st['seen'].add(X)
            poly_val=X*X-4*n*K
            assert poly_val==a*y*y
            g=math.gcd((X0*y-X*root0)%n,n)
            if g==1 or g==n:                                 # same square root of a mod n as the SIQS relation: nothing new
                st['same']+=1
                continue
            local_factors, value = factorise_fast(poly_val,primelist_f)
            ret_array[1].append(X*X)
            ret_array[0].append(poly_val)
            ret_array[2].append(local_factors)
            ret_array[3].append([])
            print("[i]psieve2 hit: a: "+str(a)+" K: "+str(K)+" square part bits: "+str(bitlen(y*y))+" of "+str(bitlen(poly_val))+" (sqrt(n): "+str(bitlen(math.isqrt(n)))+" bits) odd factors: "+str(sorted(local_factors))+" | multipliers tried: "+str(st['k'])+" hits on the known root (skipped): "+str(st['same']))
            return (g,n//g)
    st['order'].park()                                  # keep memory flat while other leftovers get their turn
    return None

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
