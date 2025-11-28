#!python
#cython: language_level=3
# cython: profile=False
# cython: overflowcheck=False
###Author: Essbee Vanhoutte
###WORK IN PROGRESS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
###Improved QS Variant 

##References: I have borrowed many of the optimizations from here: https://stackoverflow.com/questions/79330304/optimizing-sieving-code-in-the-self-initializing-quadratic-sieve-for-pypy


###To build: python3 setup.py build_ext --inplace


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

key=0                 #Define a custom modulus to factor
build_workers=8
keysize=150           #Generate a random modulus of specified bit length
workers=1 #max amount of parallel processes to use
quad_co_per_worker=1 #Amount of quadratic coefficients to check. Keep as small as possible.
base=1_000
small_base=1000
qbase=1_00
quad_sieve_size=10
g_debug=0 #0 = No debug, 1 = Debug, 2 = A lot of debug
g_lift_lim=0.5
thresvar=40  ##Log value base 2 for when to check smooths with trial factorization. Eventually when we fix all the bugs we should be able to furhter lower this.
lp_multiplier=2
min_prime=1
g_enable_custom_factors=0
g_p=107
g_q=41
mod_mul=0.5
g_max_exp=2
max_smooth_bits=30

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

cdef get_gray_code(n):
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

def formal_deriv(y,x,z):
    result=(z*2*x)-(y)
    return result

def find_r(mod,total):
    mo,i=mod,0
    while (total%mod)==0:
        mod=mod*mo
        i+=1
    return i
        
def QS(n,factor_list,sm,xlist,flist,quad_flist,z_plist):
    g_max_smooths=small_base*1+2+qbase+2
    if len(sm) > g_max_smooths: 
        del sm[g_max_smooths:]
        del xlist[g_max_smooths:]
        del flist[g_max_smooths:]  
    M2 = build_matrix(factor_list, sm, flist,quad_flist,z_plist)
    null_space=solve_bits(M2)
    f1,f2=extract_factors(n, sm, xlist, null_space)
    if f1 != 0:
        print("[SUCCESS]Factors are: "+str(f1)+" and "+str(f2))
        return f1,f2   
    print("[FAILURE]No factors found")
    return 0

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
        sqrt_left = math.isqrt(prod_left)
        ###Debug shit, remove for final version
        sqr1=prod_left%N 
        sqr2=prod_right%N
        if sqrt_right**2 != prod_right:
            print("something fucked up1")
            time.sleep(10000)
        if sqrt_left**2 != prod_left:
            print("something fucked up2")
            time.sleep(10000)
        if sqr1 != sqr2:
            print("ERROR ERROR")
            time.sleep(10000)
        ###End debug shit#########
        sqrt_left = sqrt_left % N
        sqrt_right = sqrt_right % N
        factor_candidate = gcd(N, abs(sqrt_right+sqrt_left))
     #   print(factor_candidate)
        if factor_candidate not in (1, N):
            other_factor = N // factor_candidate
            return factor_candidate, other_factor

    return 0, 0

def solve_bits(matrix):
    n=small_base+2+qbase+2
    lsmap = {lsb: 1 << lsb for lsb in range(n)}
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
        shift = n - col - 1
        val = 1 << shift
        fin = val
        for v in matrix:
            if v & val:
                fin |= v & mark_mask
        nulls.append(fin)
        k += 1
        if k == 10000000000: 
            break
    return nulls

def build_matrix(factor_base, smooth_nums, factors,quad_flist,z_plist):
    fb_map = {val: i for i, val in enumerate(factor_base)}
    fb_map2 = {val: i for i, val in enumerate(z_plist)}
    ind=1

    M2=[0]*(small_base*1+2+qbase+2)
    for i in range(len(smooth_nums)):
        for fac in factors[i]:
            idx = fb_map[fac]
            M2[idx] |= ind
        ind = ind + ind


    offset=(small_base+2)-1
    ind=1
    for i in range(len(quad_flist)):
        for fac in quad_flist[i]:
            idx = fb_map2[fac]
            M2[idx+offset] |= ind
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
    find_comb(n,complete_hmap,primeslist1,primeslist2,small_primeslist)

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


def lift(exp,co,r,n,z,z2,prime):
   # print("r: "+str(r)+" co: "+str(co)+" z: "+str(z))
    ##TO DO: Code duplication with lift_root()
    z_inv=modinv(z,prime**exp)
    temp_r=r*z
    zz=(2*temp_r)%prime
    zz=pow(zz,-1,prime)
    x=(((z*n)-temp_r**2)//prime)%prime
    y=(x*zz)%prime
    new_r=(z_inv*(temp_r+y*prime))%prime**exp
    co2=(formal_deriv(0,new_r,z))%(prime**exp) 
    ret=[]
    ret.extend([co2,new_r])
    new_eq=(z*new_r**2-n)%prime**exp
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
   # print("r: "+str(r)+" prime: "+str(prime)+" z: "+str(z))
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

@cython.profile(False)
def solve_roots(prime,n): 
    s=1  
    while jacobi((s*n)%prime,prime)!=1:
        s+=1
    if s > 50:
        print("big s")
    z_div=modinv(s,prime)  ##To do: Dont need this if we restrict to z=1
    main_root=tonelli((n*z_div)%prime,prime)
    if (s*main_root**2-n)%prime !=0:
        print("fatal error123: ")
        time.sleep(10000000)
    try:
        size=prime*2+1
        if size > 100:#quad_sieve_size*2:
            size= 100#quad_sieve_size*2+1
        temp_hmap = array.array('Q',[0]*size) ##Got to make sure the allocation size doesn't overflow.... 
        temp_hmap[0]=1
      #  while s < prime and s < quad_sieve_size+1:#2: ##To do: Dont need this if we restrict to z=1
        s_inv=modinv(s*z_div,prime)
        if s_inv == None or jacobi(s_inv,prime)!=1:
            print("should this ever happen?")
            return temp_hmap
        root_mult=tonelli(s_inv,prime)
        new_root=((main_root*root_mult))%prime
        if (s*new_root**2-n)%prime !=0:
            print("error2")
        new_co=(2*s*new_root)%prime
        if (s*new_root**2-new_co*new_root+n)%prime !=0: ###To do: For debug delete later
            print("error")
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


@cython.profile(False)
def create_hashmap(n,procnum,return_dict,primeslist):
    i=0
    hmap=[]
    while i < len(primeslist):
        hmap_p=solve_roots(primeslist[i],n)
        hmap.append(hmap_p)
        i+=1
    return_dict[procnum]=[procnum,hmap]
    return 


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

def equation2(y,x,n,mod,z,z2):
    rem=z*(x**2)+y*x-n*z2
    rem2=rem%mod
    return rem2,rem



@cython.boundscheck(False)
@cython.wraparound(False)
cdef factorise_fast(value,long long [::1] factor_base):
    seen_primes=[]
    factors = set()
    if value < 0:
        seen_primes.append(-1)
        factors ^= {-1}
        value = -value
    while value % 2 == 0:
        seen_primes.append(2)
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




cdef get_quads(hmap,indexes,cfact):
    quads=[]
    lin_lists=[]
    i=0
    while i < len(indexes):
        lin_lists.append([])
        index=indexes[i]
        quads.append([])
        length=hmap[index][0]
        j=1
        while j < length:
            lin_lists[-1].append(hmap[index][j+1])
            quads[-1].append(hmap[index][j])
            j+=2
        i+=1
    return quads,lin_lists

cdef create_quad_cos(quads,lin_lists,cfact,local_mod):
    i=0
    parts=[]
    parts2=[]
    length=len(cfact)
    while i < length:
        parts.append([])
        parts2.append([])
        p=cfact[i]
        j=0
        while j < len(quads[i-1]):
            r1 = quads[i-1][j]
            aq = local_mod // p
            invaq = modinv(aq%p, p)
            gamma = r1 * invaq % p
            parts[-1].append(aq*gamma)


            r1 = lin_lists[i-1][j]
            gamma = r1 * invaq % p
            parts2[-1].append(aq*gamma)
            j+=1
        i+=1
    return parts,parts2





cdef get_lin(hmap,cfact,local_mod,indexes,quad_co,n):
    
    all_lin_parts=[]
    j=0
    lin=0
    while j < len(indexes):
        ind=indexes[j]
        prime=math.isqrt(cfact[j])
 

        r1=hmap[ind][2]

        r1=lift_b(prime,n,r1,quad_co,2)
        prime=prime**2
        aq = local_mod // prime
        invaq = modinv(aq%prime, prime)
        gamma = r1 * invaq % prime
        lin+=aq*gamma
        all_lin_parts.append(aq*gamma)
        j+=1
    lin%=local_mod
    return lin,all_lin_parts

def gen_co2(factor_base,hmap,cmod,n):

    length=factor_base[0]
    co2_list=[]
    i=1
    while i < length:
        prime=factor_base[i]

        
        prime_index=i-1
   
 

        co2=hmap[prime_index][2]%prime
        s=hmap[prime_index][1]%prime
        co2_list.append([co2,s])
     
        i+=1
    return co2_list


def filter_quads(qbase):
    valid_quads=[]
    valid_quads_factors=[]
    i=1
    while i < quad_sieve_size+1:
        if i != 1 and math.sqrt(i)%1==0:
            i+=1
            continue
        quad_local_factors, quad_value,seen_primes = factorise_fast(i,qbase)  ##TO DO: move this way up

        if quad_value != 1:
            i+=1
            continue
        valid_quads.append(i)
        valid_quads_factors.append(quad_local_factors)
        i+=1
    return valid_quads,valid_quads_factors


def iterate_mod_co(orig_lin_parts,cfact,grays,new_mod,n,quad):
    lin_co_array=[]
    lin_parts=[]
    q=0
    prev_co=1
    lin3=0
    skip=0
    while q < len(orig_lin_parts):
        prime=cfact[q]
        orig_lin_temp=orig_lin_parts[q]
        z_inv=modinv(quad,prime)
        prime_sqrt=math.isqrt(prime)
        if z_inv == None or jacobi(z_inv,prime_sqrt)!=1:
            skip=1
            break
                    
        root_mult=tonelli(z_inv%prime_sqrt,prime_sqrt)##to do: fix this
        if root_mult**2%prime_sqrt != z_inv%prime_sqrt:
            print("what the fuuk11111111111111111111111111")  


        new_root=((orig_lin_temp*root_mult))%prime
        new_co=(quad*new_root)%prime
        if (new_co**2-n*4*quad)%prime_sqrt !=0:
            print("catastrophic error")
        if new_co > prime // 2:
            new_co=(prime-new_co)%prime  
        new_co=lift_b(prime_sqrt,n,new_co,quad,2)
        aq = new_mod // prime
        invaq = modinv(aq%prime, prime)
        gamma = new_co * invaq % prime
        lin3+=(aq*gamma)
        if (aq*gamma)%prime !=new_co:
            print("something went wrong")
        lin_parts.append(aq*gamma)
        q+=1
     
    if skip == 1:
        return 0
    lin2=lin3%new_mod
    poly_ind=0
    end = 1 << (len(cfact) - 1)
    lin=0
    while poly_ind < end:
        if poly_ind != 0:
            v,e=grays[poly_ind] ##Brilliant trick, whoever come up with this.
            lin=(lin + 2 * e * lin_parts[v])%new_mod

        else:
            lin=lin2
        if (lin**2-n*4*quad)%new_mod !=0:
            print("super big error")
            time.sleep(1000)
        lin_co_array.append(lin)
        poly_ind+=1

    return lin_co_array   

def iterate2(co2_item,prime,valid_quads,orig_lin_parts,cfact,grays,n,new_mod,res):
    i=0
    while i < len(valid_quads):
        quad_co=valid_quads[i]
        if gcd(new_mod,quad_co)!=1:
            i+=1
            continue
        lin_co_array=iterate_mod_co(orig_lin_parts,cfact,grays,new_mod,n,quad_co)
        if lin_co_array == 0:
            i+=1
            continue
        root=get_root(new_mod,lin_co_array[0],quad_co)
        temp_quad=co2_item[1]
        z_div=modinv(temp_quad,prime)
        z_inv=modinv(quad_co*z_div,prime)
        if z_inv == None or jacobi(z_inv,prime)!=1:
            i+=1
            continue
        
        root_mult=tonelli(z_inv,prime)
        temp_co=co2_item[0]
        r=get_root(prime,temp_co,temp_quad)
        r=(r*root_mult)%prime
        r_b=(prime-r)%prime  
        co2=(2*quad_co*r)%prime   
        if co2 > prime // 2:
            co2=(prime-co2)%prime  
        co2_b=(prime-co2)%prime
        root2=r 
        root3=r_b
        #lcmod=new_mod%prime
        exp=1
        while exp < 20:

            if exp > 1:
                root2=lift_root(r,prime**(exp-1),n,quad_co,exp)
            root_dist1=solve_lin_con(new_mod,r-root,prime**exp)
            if root_dist1 > 1000:
                break
            if exp%2 == 0:
                res[i]+=1
      #  modi=modinv(lcmod,prime)
       # root_res=(r-root)%prime
        #root_dist1=(root_res*modi)%prime
       # if root_dist1 < 10:
          #  res[i]+=1
       # root2=lift_root(r,prime,n,quad_co)
       # root_dist1=solve_lin_con(new_mod,r-root,prime)
       # if prime ==f:
                print("r: "+str(r)+" prime: "+str(prime)+" root_dist1: "+str(root_dist1)+" quad_co: "+str(quad_co)+" root_mult: "+str(root_mult)+" lin: "+str(lin_co_array[0])+" mod: "+str(new_mod))
            exp+=1
      #  print("r_b: ",r_b)
#To do: using just one cmod co, loop and see if we can spot a shortcut for root_dist1.
#To do: Can potentially simplify calculations by only using primes at n*4
        i+=1
    return

def iterate(co2_list,valid_quads,new_mod,primelist_f,orig_lin_parts,cfact,grays,n):
    length=primelist_f
    res=[0]*len(valid_quads)
    r=1+small_base
    while r < len(length):
       # if primelist_f[r] ==primelist_f[-1]:
        iterate2(co2_list[r-1],primelist_f[r],valid_quads,orig_lin_parts,cfact,grays,n,new_mod,res)
  


        r+=1
    print("res: ",res)
  #  return

def grab_primes_for_quad(primeslist,quad,n):
    lprimes=[]
    lprimes_indexes=[]
    i=0
    while i < len(primeslist):
        if jacobi(quad*n,primeslist[i])==1:
            lprimes.append(primeslist[i])
            lprimes_indexes.append(i)
        i+=1
    return lprimes,lprimes_indexes
cdef construct_interval(list ret_array,partials,n,primeslist,hmap,large_prime_bound,primeslist2,small_primeslist):
    grays = get_gray_code(20)

    cdef Py_ssize_t i
    cdef Py_ssize_t j
    close_range =30
    too_close = 1   ##TO DO: remove this small prime.. better to not have small primes 
    LOWER_BOUND_SIQS=1

    cdef Py_ssize_t size
    primelist=copy.copy(primeslist)
    primelist.insert(0,2) ##To do: remove when we fix lifting for powers of 2
    primelist.insert(0,-1)
    z_plist=copy.copy(primeslist2)
    z_plist.insert(0,len(primeslist2)+1)
    z_plist=array.array('q',z_plist)
    valid_quads,valid_quads_factors=filter_quads(z_plist)
    primeslist2.insert(0,2)
    primeslist2.insert(0,-1)


    primeslist=array.array('q',primeslist)

    primelist_f=copy.copy(primeslist)
    primelist_f.insert(0,len(primelist_f)+1)
    primelist_f=array.array('q',primelist_f)

    sprimelist_f=copy.copy(small_primeslist)
    sprimelist_f.insert(0,len(sprimelist_f)+1)
    sprimelist_f=array.array('q',sprimelist_f)

    sprimelist=copy.copy(small_primeslist)
    sprimelist.insert(0,2) ##To do: remove when we fix lifting for powers of 2
    sprimelist.insert(0,-1)

    mod_found=0
    mod2_found=0

        



    tnum = int(((2*n)**mod_mul))
    tnum_bit=int(bitlen(tnum)*1)
   
    
    
    
    r=0
    while r < len(valid_quads):
        quad=valid_quads[r]
        seen=[]  ###Moduli seen in generate modulus
        lprimes,lprimes_indexes=grab_primes_for_quad(primeslist,quad,n)
        while 1:
            root_list,cmod_list,cfact_list=generate_modulus(n,lprimes,seen,tnum,close_range,too_close,LOWER_BOUND_SIQS,tnum_bit,hmap,lprimes_indexes,quad)
            if root_list ==0:
                break
           # print("new_mod: "+str(new_mod)+" cfact: "+str(cfact)+" bitlen mod: "+str(bitlen(new_mod))+" bitlen root: "+str(bitlen(root)))
            g=0
            while g < len(root_list):
                root=root_list[g]
                new_mod=cmod_list[g]
                cfact=cfact_list[g]
                test=(quad*root**2-n)%new_mod
                if test !=0:
                    print("something went wrong")
                mod_found+=process_interval(ret_array,n,quad,root,sprimelist_f ,new_mod,valid_quads_factors[r],partials,large_prime_bound,primelist_f)
                if mod_found+1 > 1:  
                    print("", end=f"[i]Smooths: {len(ret_array[0])} / {small_base*1+2+qbase}\r")
                    mod2_found+=mod_found
                    mod_found=0
                    if mod2_found+1 > 100:
                        print("[i]Attempting linear algebra")
                        test=QS(n,sprimelist,ret_array[0],ret_array[1],ret_array[2],ret_array[3],primeslist2)   
                        if test!=0:
                            return 
                        mod2_found=0
                    if len(ret_array[0]) > small_base*1+2+qbase+2:
                        test=QS(n,sprimelist,ret_array[0],ret_array[1],ret_array[2],ret_array[3],primeslist2)  
                        return
                g+=1
        r+=1
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
def find_quads(local_factors,hmap,quad_interval_index,quad_interval,n,poly_val,factor_base):
    local_factors=list(local_factors)
   # print(local_factors)
    length=quad_interval[0]
    primes=[]
    enum_quad=[]
    enum_lin=[]


   
    c=0
    while c < len(local_factors): ###TO DO: THIS SUCKS, USE AN INDEX MAPPING FOR FUCK SAKE
        if local_factors[c]==-1 or local_factors[c] ==2:
            c+=1
            continue
        if poly_val%local_factors[c] == 0:
        
                #print("HIT")
            enum_quad.append(local_factors[c])
            enum_lin.append(local_factors[c])
            length=factor_base[0]
      
            p=1
            while p < length:
                if local_factors[c]==factor_base[p]:
                    p-=1
                    break
                p+=1

            enum_quad.append([])
            enum_lin.append([])
            length2=hmap[p][0]
            j=1
            while j < length2:
  
                enum_quad[-1].append(hmap[p][j]%local_factors[c])
                enum_quad[-1].append(hmap[p][j]%local_factors[c])
                enum_lin[-1].append(hmap[p][j+1]%local_factors[c])
                enum_lin[-1].append(hmap[p][j+1]%local_factors[c])

                j+=2
    
        c+=1


    enum_quad=get_partials(poly_val,enum_quad)
    enum_lin=get_partials(poly_val,enum_lin)

    return enum_quad,enum_lin
def enumerated_product(*args):
    yield from itertools.product(*(range(len(x)) for x in args))


cdef lift_root(r,prime,n,quad_co,exp):
    z_inv=modinv(quad_co,prime**exp)
    c=(quad_co*n)%prime**exp
    temp_r=r*quad_co
    zz=2*temp_r
    zz=pow(zz,-1,prime)
    x=((c-temp_r**2)//prime)%prime
    y=(x*zz)%prime
    new_r=(temp_r+y*prime)%prime**exp
    root2=(new_r*z_inv)%prime**exp
    return root2




#@cython.boundscheck(False)
#@cython.wraparound(False)
cdef process_interval(ret_array,n,quad_co2,root,long long [::1] local_primes,cmod,quad_local_factors,partials,large_prime_bound,primelist_f):
    


    mod_found=0

    
    quad_co=quad_co2
 

  
    poly_val=quad_co*root**2-n
            
    if poly_val%cmod!=0:
        print("very fatal error")
    poly_val2 =poly_val//cmod
   # print("(pos) bitlen polyval/mod: "+str(bitlen(abs(poly_val2)))+" j: "+str(j)+" quad_co: "+str(quad_co))
    co=quad_co*root**2
    local_factors, value,seen_primes = factorise_fast(poly_val2,local_primes)
    quad_local_factors2=copy.deepcopy(quad_local_factors)
    test=math.isqrt(value)
    if test**2 != value:
        if value < large_prime_bound:
            if value in partials:
                rel, lf, pv,ql = partials[value]
                if rel == co:
                    return 0
                co *= rel
                local_factors ^= lf
                poly_val *= pv
                quad_local_factors2 ^=ql
                print("hit large prime")
            else:
                partials[value] = (co, local_factors, poly_val,quad_local_factors2)
                return 0
        else:
            return 0

    if co not in ret_array[1]:
        print("Found: "+str(seen_primes)+" cmod bitlen: "+str(bitlen(cmod))+" bitlen root: "+str(bitlen(root))+" quad: "+str(quad_co)+" Smooth#: "+str(len(ret_array[0])+1)+" bitlen poly_val/cmod: "+str(bitlen(abs(poly_val//cmod)))+" remaining square: "+str(test)+" root: "+str(root)+" mod: "+str(cmod))
        mod_found+=1
        ret_array[0].append(poly_val)
        ret_array[1].append(co)
        ret_array[2].append(local_factors)
        ret_array[3].append(quad_local_factors2)
 
    return mod_found

cdef enum_stuff(enum_quad,enum_lin2,quad,new_mod,enum_quad2):

    indexes=[]
    lin_parts=[]
    i=0
    while i < len(enum_quad2):
        j=0
        index=-1
        while j < len(enum_quad2[i]):
            if enum_quad2[i][j]%enum_quad[i*2]==quad%enum_quad[i*2]:
                index=j
                break
            j+=1
        if index == -1:
            break
        indexes.append(index)
        i+=1
    if len(indexes)!=len(enum_quad2):
        return 0,0
    l=0
    while l < len(enum_quad2):
        lin_parts.append(enum_lin2[l][indexes[l]])
        l+=1
    return lin_parts

def gen_co2(factor_base,hmap,cmod,n):

    length=factor_base[0]
    co2_list=[]
    i=1
    while i < length:
        prime=factor_base[i]

        
        prime_index=i-1
        z=1


        s=1  
        while jacobi((s*n)%prime,prime)!=1:
            s+=1
 

        co2=hmap[prime_index][2]%prime
        co2_list.append([co2,s])
     
        i+=1
    return co2_list


def increase_coefficient(prime,co,temp_quad,cmod,quad,n,total_root):
    target=(cmod*prime**2)
    z_div=modinv(temp_quad,prime)
    z_inv=modinv(quad*z_div,prime)
    if z_inv == None or jacobi(z_inv,prime)!=1:
        return 0
    root_mult=tonelli(z_inv,prime)
    root=get_root(prime,co,temp_quad)
    root=(root*root_mult)%prime
    root=lift_root(root,prime,n,quad,2)
    if (quad*root**2-n)%(prime**2)!=0:
        print("er")
    if cmod !=1:
        root2=(-root)%prime**2
               # print("cmod: "+str(cmod)+" total_root: "+str(total_root)+" new prime: "+str(potential_a_factor))
        dist=solve_lin_con(cmod,root-total_root,prime**2)
        

        dist2=solve_lin_con(cmod,root2-total_root,prime**2)
        if dist2 < dist:
            total_root=total_root+cmod*dist2
        else:
            total_root=total_root+cmod*dist


        if (quad*total_root**2-n)%cmod!=0:
            print("er11")
    else:
        total_root=root
    return total_root


@cython.boundscheck(False)
@cython.wraparound(False)
cdef generate_modulus(n,primeslist,seen,tnum,close_range,too_close,LOWER_BOUND_SIQS,tnum_bit,hmap,primeslist_indexes,quad):
    cdef Py_ssize_t counter 
    cdef Py_ssize_t counter2
    cdef Py_ssize_t counter3
    cdef Py_ssize_t counter4
    cdef Py_ssize_t i
    cdef Py_ssize_t j
    cdef Py_ssize_t const_1=1_00
    cdef Py_ssize_t const_2=1_000

    small_B = len(primeslist)
    lower_polypool_index = 2
    upper_polypool_index = len(primeslist) - 1
    poly_low_found = False
    total_root=0
    for i in range(small_B):  ##To do: Can be moved outside mainloop
        if primeslist[i] > LOWER_BOUND_SIQS and not poly_low_found:
            lower_polypool_index = i
            poly_low_found = True
            break

    small_B=upper_polypool_index
    counter4=0
    while counter4 < const_1:
        counter4+=1
        cmod = 1
        cfact = []
        counter2=0
        while counter2 < const_2:
            found_a_factor = False
            counter=0
            while(found_a_factor == False) and counter < const_2:
                randindex = random.randint(lower_polypool_index, upper_polypool_index)
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


            total_root=increase_coefficient(primeslist[randindex],hmap[primeslist_indexes[randindex]][2],hmap[primeslist_indexes[randindex]][1],cmod,quad,n,total_root)
            if total_root ==0:
                cmod = 1
                s = 0
                cfact = []
                continue
            cmod = cmod * potential_a_factor
            if (quad*total_root**2-n)%cmod !=0:
                print("ERRRRRRRRRRRRRRRRRRROR: ",cmod)
                time.sleep(100000)
            cfact.append(potential_a_factor)
            j = tnum_bit - cmod.bit_length()
            counter2+=1
            if j < too_close:
                cmod = 1
                s = 0
                cfact = []#[0]*base
               # indexes=[]
                continue
            elif j < (too_close + close_range):
                break



 
        #print("intermediate: "+str(bitlen(cmod))+" bitlen root: "+str(bitlen(total_root))+" from goal: "+str(bitlen(quad*total_root**2)-bitlen(int(n))))#+" mindiff: "+str(mindiff))
        root_list=[]
        cmod_list=[]
        cfact_list=[]



        i=0
        while i < len(primeslist):

            if primeslist[i]**2 in cfact:
                i+=1
                continue
            total_root_temp=increase_coefficient(primeslist[i],hmap[primeslist_indexes[i]][2],hmap[primeslist_indexes[i]][1],cmod,quad,n,total_root)
            if total_root_temp ==0:
                continue
            if bitlen(abs((quad*total_root_temp**2-n)//(cmod * primeslist[i]**2)))<max_smooth_bits:
                if (cmod * primeslist[i]**2) in seen:
                    i+=1
                    continue
                cmod_temp = cmod * primeslist[i]**2

                if (quad*total_root_temp**2-n)%cmod_temp !=0:
                    print("ERRRRRRRRRRRRRRRRRRROR:55555555555555555555 ",cmod_temp)
                    time.sleep(100000)

                cfact_temp=copy.copy(cfact)
                cfact_temp.append(primeslist[i]**2)
                root_list.append(total_root_temp)
                cmod_list.append(cmod_temp)
                cfact_list.append(cfact_temp)
            i+=1
        if len(root_list)!=0:
            #print(len(root_list))
            return root_list,cmod_list,cfact_list
        continue


    return 0,0,0









def find_comb(n,hmap,primeslist1,primeslist2,small_primeslist):
    start=default_timer()
    ret_array=[[],[],[],[]]

    
    if g_debug > 0:
        duration = default_timer() - start
        print("Constructing quad interval took: "+str(duration))
    partials={}
    gc.enable()

    large_prime_bound = primeslist1[-1] ** lp_multiplier
    start=default_timer()



    if g_debug > 0:
        duration = default_timer() - start
        print("Processing quad interval took: "+str(duration))
    construct_interval(ret_array,partials,n,primeslist1,hmap,large_prime_bound,primeslist2,small_primeslist)

    return 0

def get_primes(start,stop):
    return list(sympy.sieve.primerange(start,stop))

@cython.profile(False)
def main(l_keysize,l_workers,l_debug,l_base,l_key,l_quad_sieve_size,sbase):
    global key,keysize,workers,g_debug,base,key,quad_sieve_size,small_base
    key,keysize,workers,g_debug,base,quad_sieve_size,small_base=l_key,l_keysize,l_workers,l_debug,l_base,l_quad_sieve_size,sbase
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
        if len(small_primeslist)<small_base:
            small_primeslist.append(primeslist[i])
    
        primeslist1.append(primeslist[i])

        i+=1
    i=0
    while len(primeslist2) < qbase:
        primeslist2.append(primeslist[i])
        i+=1
    launch(n,primeslist1,primeslist2,small_primeslist)     
    duration = default_timer() - start
    print("\nFactorization in total took: "+str(duration))




