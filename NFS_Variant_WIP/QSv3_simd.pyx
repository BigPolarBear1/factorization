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
qbase=10
quad_sieve_size=10
g_debug=0 #0 = No debug, 1 = Debug, 2 = A lot of debug
g_lift_lim=0.5
thresvar=1  ##Log value base 2 for when to check smooths with trial factorization. Eventually when we fix all the bugs we should be able to furhter lower this.
lp_multiplier=2
min_prime=1
g_enable_custom_factors=0
g_p=107
g_q=41
mod_mul=0.05
g_max_exp=10
lin_sieve_size=1000


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
        
def QS(n,factor_list,sm,xlist,flist,quad_flist,z_plist,modk):
    g_max_smooths=small_base*1+2+qbase+2
    if len(sm) > g_max_smooths: 
        del sm[g_max_smooths:]
        del xlist[g_max_smooths:]
        del flist[g_max_smooths:]  
    M2 = build_matrix(factor_list, sm, flist,quad_flist,z_plist)
    null_space=solve_bits(M2)
    f1,f2=extract_factors(n, sm, xlist, null_space,modk)
    if f1 != 0:
        print("[SUCCESS]Factors are: "+str(f1)+" and "+str(f2))
        return f1,f2   
    print("[FAILURE]No factors found")
    return 0

def extract_factors(N, relations, roots, null_space,modk):
    n = len(relations)
    for vector in null_space:
        prod_left = 1
        prod_right = 1
        print("Checking")
        for idx in range(len(relations)):
            bit = vector & 1
            vector = vector >> 1
            if bit == 1:
                prod_left *= roots[idx]
                print("adding left (zx^2): "+str(roots[idx]))
                prod_right *= relations[idx]
                print("adding right (zx^2+yx-n): "+str(relations[idx]))
                print("adding poly: "+str(modk[idx][3])+"*"+str(modk[idx][2])+"^2+"+str(modk[idx][0])+"*"+str(modk[idx][1])+"*"+str(modk[idx][2])+"-"+str(N)+" = "+str(relations[idx]))#str(modk[idx]))
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
            print("ERROR: FIXME... we need to adjust the result... perhaps some type of squar root over a finite field like NFS")
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
    n=small_base+2+qbase+2
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
def solve_roots(prime,n): 
    s=1  
    while jacobi((s*n)%prime,prime)!=1:
        s+=1
    if s > 100:
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

cdef filter_quads(qbase,n):
    valid_quads=[]
    valid_quads_factors=[]
    interval_list=[]
    interval_list_pos=[]
    roots=[]
    i=1
    while i < quad_sieve_size+1:
        if i != 1 and math.sqrt(i)%1==0:
            i+=1
            continue
        if i ==0:
            print("wtf")
        quad_local_factors, quad_value,seen_primes = factorise_fast(i,qbase)  ##TO DO: move this way up

        if quad_value != 1:
            i+=1
            continue

        interval_list.append(array.array('i', [0]*lin_sieve_size))
        interval_list_pos.append(array.array('i', [0]*lin_sieve_size))

        roots.append(math.isqrt(n//i))
        valid_quads.append(i)
        valid_quads_factors.append(quad_local_factors)
        i+=1
    return valid_quads,valid_quads_factors,interval_list,roots,interval_list_pos







cdef sieve(interval_list,valid_quads,roots,n,sprimelist_f,hmap,interval_list_pos):
    length=sprimelist_f[0]
    i=1
    small=1
    while i < length:
        if i > small_base:
            small=0
        prime=sprimelist_f[i]
        prime_index=i-1

      #  print("prime: ",sprimelist_f[i])
        quad=hmap[prime_index][1]%prime
        z_div=modinv(quad,prime)
        
        
        co=hmap[prime_index][2]%prime
        
        r=get_root(prime,co,quad)
        if (quad*r**2-r*co+n)%prime !=0:
            print("error")
       # r=(r*root_mult)%prime
       # r_b=(prime-r)%prime  
        #to do: First calculate the root here
        j=0
        while j < len(valid_quads):
            quad2=valid_quads[j]
            root=roots[j]
            #print("root: ",root)
            if quad2%prime==0:
                j+=1
                continue
            z_inv=modinv(quad2*z_div,prime)
            if z_inv == None or jacobi(z_inv,prime)!=1:
                j+=1
                continue
            root_mult=tonelli(z_inv,prime)
            r2=(r*root_mult)%prime
            r_b=(prime-r2)%prime  
            log=int(math.log2(prime))
            if small==0:
                log*=2
            exp=1
            while exp < g_max_exp+1:

                if exp > 1:
                    r2=lift_root(r2,prime**(exp-1),n,quad2,exp)
                if (small ==0 and exp%2==0) or small ==1:
                    if (quad2*r2**2-n)%prime**exp !=0:
                        print("error2")   
                    dist=((root%prime**exp)-r2)%prime**exp
               # if prime ==3:
                 #   print("dist: "+str(dist)+" exp: "+str(exp)+" prime: "+str(prime)+" quad: "+str(quad2))
                    new_root=root+(-dist%prime**exp)

                    if (quad2*new_root**2-n)%prime**exp !=0:
                        print("error2")  
                        time.sleep(100000) 
                    if dist < lin_sieve_size:
                        miniloop_non_simd(dist,interval_list[j],prime**exp,log,lin_sieve_size)
                    dist2=-dist%prime**exp
                    if dist2 < lin_sieve_size:
                        miniloop_non_simd(dist2,interval_list_pos[j],prime**exp,log,lin_sieve_size)


                    
                    r_b=(-r2)%prime**exp#lift_root(r_b,prime**(exp-1),n,quad2,exp)
                    dist3=((root%prime**exp)-r_b)%prime**exp

                    if dist3 < lin_sieve_size:
                        miniloop_non_simd(dist3,interval_list[j],prime**exp,log,lin_sieve_size)
  
                    dist4=-dist3%prime**exp
                    if dist4 < lin_sieve_size:
                        miniloop_non_simd(dist4,interval_list_pos[j],prime**exp,log,lin_sieve_size)

                    if (quad2*new_root**2-n)%prime**exp !=0:
                        print("error2")  
                        time.sleep(100000) 
                    if dist >lin_sieve_size and dist2 >lin_sieve_size and dist3 >lin_sieve_size and dist4 >lin_sieve_size:
                        break
                exp+=1
            #to do: Then mutate the root here
            j+=1
        i+=1
    print("done")
    return

cdef construct_interval(list ret_array,partials,n,primeslist,hmap,large_prime_bound,primeslist2,small_primeslist):
    





    cdef Py_ssize_t i
    cdef Py_ssize_t j
    close_range =5
    too_close = 1   ##TO DO: remove this small prime.. better to not have small primes 
    LOWER_BOUND_SIQS=1
    primelist_f=copy.copy(primeslist)
    primelist_f.insert(0,len(primelist_f)+1)
    primelist_f=array.array('q',primelist_f)
    cdef Py_ssize_t size
    primelist=copy.copy(primeslist)
    primelist.insert(0,2) ##To do: remove when we fix lifting for powers of 2
    primelist.insert(0,-1)
    z_plist=copy.copy(primeslist2)
    z_plist.insert(0,len(primeslist2)+1)
    z_plist=array.array('q',z_plist)
    sprimelist_f=copy.copy(small_primeslist)
    sprimelist_f.insert(0,len(sprimelist_f)+1)
    sprimelist_f=array.array('q',sprimelist_f)
    valid_quads,valid_quads_factors,interval_list,roots,interval_list_pos=filter_quads(z_plist,n)
    print("[i]Filling in the intervals... this can take a while..")
    sieve(interval_list,valid_quads,roots,n,primelist_f,hmap,interval_list_pos)
    print("[i]Checking intervals for smooths")
    primeslist2.insert(0,2)
    primeslist2.insert(0,-1)


    primeslist=array.array('q',primeslist)





    sprimelist=copy.copy(small_primeslist)
    sprimelist.insert(0,2) ##To do: remove when we fix lifting for powers of 2
    sprimelist.insert(0,-1)

    mod_found=0
    mod2_found=0

    many_primes=[]    


    tnum = int(((2*n)**mod_mul))
    tnum_bit=int(bitlen(tnum)*1)
    rootlist,polylist,flist,quadflist,modk=generate_large_square(n,many_primes,valid_quads,valid_quads_factors,sprimelist_f,interval_list,roots,interval_list_pos,partials,large_prime_bound)
    print("\nSmooths found: "+str(len(rootlist)))#+" rootlist: "+str(rootlist))
    test=QS(n,sprimelist,polylist,rootlist,flist,quadflist,primeslist2,modk)  
    return

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







 

@cython.boundscheck(False)
@cython.wraparound(False)
cdef void miniloop_non_simd(dist1,temp,prime,log,size):
    while dist1 < size :
        temp[dist1]+=log
        dist1+=prime
    return
  

cdef factorise_squares(value,factor_base):
    seen_primes=[]
    total_square=1
    exp=0
    seenp=1
    while value % 2 == 0:
        seenp*=2
        value //= 2
        exp+=1
    if exp>1:
        if exp%2==1:
            exp-=1
        total_square*=2**exp
        seen_primes.append(2**exp)
    #cdef int factor 
    i=0
    while i < len(factor_base):
        factor=factor_base[i]
        exp=0
        seenp=1
        while value % factor == 0:
            seenp*=factor
            value //= factor
            exp+=1

        if exp>1:
            if exp%2==1:
                exp-=1
            seen_primes.append(factor**exp)
            total_square*=factor**exp
        i+=1
    return value,seen_primes,total_square

def reduce_smooth_size(seen_primes,quad,root,n):
    print("seen_primes: ",seen_primes)
    mod=1
    for p in seen_primes:
        mod*=1
    print("mod: ",mod)
    poly_val=quad*root**2-n
    i=0
    while i < 1000:
        y=mod*i
        poly_val=quad*root**2+y*root-n
        print("poly_val:",poly_val)
        i+=1

def generate_large_square(n,many_primes,valid_quads,valid_quads_factors,sprimelist_f,interval_list,roots,interval_list_pos,partials,large_prime_bound):
    root_list=[]
    poly_list=[]
    flist=[]
    quadf_list=[]
    modk=[]
    i=0
    while i <len(valid_quads):
        quad=valid_quads[i]
        root=roots[i]

      #  j=0
      #  while j < lin_sieve_size:
      #      if interval_list[i][j]>thresvar:
      #          
      #          poly_val=quad*(root-j)**2-n
      #          tot=quad*(root-j)**2
              #  print("poly_val: ",poly_val)
      #          local_factors, value,seen_primes = factorise_fast(poly_val,sprimelist_f)

      #          quad_local_factors2=copy.copy(valid_quads_factors[i])
                  #  logged=0
                  #  for p in seen_primes:
                  #      if p != -1 and p != 2:
                  #          logged+=int(math.log2(p))
                 #   print("interval_list[i][j]: "+str(interval_list[i][j])+" seen_primes: "+str(seen_primes)+" logged: "+str(logged)+" quad: "+str(quad))
       #         test=math.isqrt(value)
       #         if test**2 != value:
       #             if value < large_prime_bound:
       #                 if value in partials:
       #                     rel, lf, pv,ql = partials[value]
       #                     if rel == tot:
       #                         j+=1
        #                        continue
        #                    tot *= rel
         #                   local_factors ^= lf
         #                   poly_val *= pv
         #                   quad_local_factors2 ^=ql
         #               else:
         #                   partials[value] = (tot, local_factors, poly_val,quad_local_factors2)
         #                   j+=1
         #                   continue
         #           else:
        #              j+=1 
        #                continue


              
                        
         #       if tot not in root_list:
         #                  # print("adding")
         #           root_list.append(tot)
         #           poly_list.append(poly_val)
         #           flist.append(local_factors)
         #           quadf_list.append(quad_local_factors2)
         #           print("", end=f"[i]Smooths: {len(root_list)} / {small_base*1+2+qbase}\r")
         #           if len(root_list)>small_base+2+qbase+2:
         #               return root_list,poly_list,flist,quadf_list  
         #   j+=1

        j=0
        while j < lin_sieve_size:
            if interval_list_pos[i][j]>10:
                poly_val2=quad*(root+j)**2-n
               # print("initital poly: ",poly_val2)

                tot=quad*(root+j)**2
                quad_local_factors2=copy.copy(valid_quads_factors[i])
                local_factors, value,seen_primes = factorise_fast(poly_val2,sprimelist_f)
                if -1 not in local_factors and 2 not in local_factors:
                    mod=1
                    for p in local_factors:
                        mod*=p
                   # print("\nmod: "+str(mod)+" local: "+str(local_factors))
                    k=1
                    while k < 1000:
                        y=mod*k
                        poly_val2=quad*(root+j)**2+y*(root+j)-n
                        tot=quad*(root+j)**2
                        quad_local_factors2=copy.copy(valid_quads_factors[i])
                        local_factors, value,seen_primes = factorise_fast(poly_val2,sprimelist_f)
                       # print("poly_val:",poly_val)
                       
                   # reduce_smooth_size(local_factors,quad,root+j,n)
                  #  time.sleep(10000)
                   # logged=0
                   # for p in seen_primes:
                    #    if p != -1 and p != 2:
                    #        logged+=int(math.log2(p))
                   # print("interval_list[i][j]: "+str(interval_list_pos[i][j])+" seen_primes: "+str(seen_primes)+" logged: "+str(logged)+" quad: "+str(quad))

                        test=math.isqrt(value)
                        if test**2 == value:
                 #   if value < large_prime_bound:
                 #       if value in partials:
                 #           rel, lf, pv,ql = partials[value]
                 #           if rel == tot:
                 #               j+=1
                  #              continue
                  #          tot *= rel
                  #          local_factors ^= lf
                  #          poly_val2 *= pv
                  #          quad_local_factors2 ^=ql
                  #      else:
                  #          partials[value] = (tot, local_factors, poly_val2,quad_local_factors2)
                  #          j+=1
                  #          continue
                  #  else:
                  #      j+=1 
                  #      continue
                        
                            if tot not in root_list:
                                root_list.append(tot)
                                poly_list.append(poly_val2)
                                flist.append(local_factors)
                                modk.append([mod,k,root+j,quad])
                                quadf_list.append(quad_local_factors2)
                                print("", end=f"[i]Smooths: {len(root_list)} / {small_base*1+2+qbase}\r")
                                if len(root_list)>small_base+2+qbase+2:
                                    return root_list,poly_list,flist,quadf_list,modk
                        k+=1  
            j+=1
            #print("Smooth candidate #1: "+str(poly_val)+" Smooth candidate #2: "+str(poly_val2)+" quad: "+str(quad+i)+" bitlength smooth can #1: "+str(bitlen(poly_val))+" bitlength smooth can #2: "+str(bitlen(poly_val2)))
        i+=1
    return root_list,poly_list,flist,quadf_list,modk










def find_comb(n,hmap,primeslist1,primeslist2,small_primeslist):

    ret_array=[[],[],[],[]]

    

    partials={}


    large_prime_bound = primeslist1[-1] ** lp_multiplier
    construct_interval(ret_array,partials,n,primeslist1,hmap,large_prime_bound,primeslist2,small_primeslist)

    return 0

def get_primes(start,stop):
    return list(sympy.sieve.primerange(start,stop))

@cython.profile(False)
def main(l_keysize,l_workers,l_debug,l_base,l_key,l_lin_sieve_size,l_quad_sieve_size,sbase):
    global key,keysize,workers,g_debug,base,key,quad_sieve_size,small_base,lin_sieve_size
    key,keysize,workers,g_debug,base,quad_sieve_size,small_base,lin_sieve_size=l_key,l_keysize,l_workers,l_debug,l_base,l_quad_sieve_size,sbase,l_lin_sieve_size
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
        if n%primeslist[i]==0:
            i+=1
            continue
        if len(small_primeslist)<small_base:
            small_primeslist.append(primeslist[i])
    
        primeslist1.append(primeslist[i])

        i+=1
    i=0
    while len(primeslist2) < qbase:
        if n%primeslist[i]==0:
            i+=1
            continue
        primeslist2.append(primeslist[i])
        i+=1
    launch(n,primeslist1,primeslist2,small_primeslist)     
    duration = default_timer() - start
    print("\nFactorization in total took: "+str(duration))




