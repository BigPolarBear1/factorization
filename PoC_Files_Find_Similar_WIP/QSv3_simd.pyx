#!python
#cython: language_level=3
# cython: profile=False
# cython: overflowcheck=False
###Author: Essbee Vanhoutte
###WORK IN PROGRESS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
###Factorization_v3 

##References: I have borrowed many of the optimizations from here: https://stackoverflow.com/questions/79330304/optimizing-sieving-code-in-the-self-initializing-quadratic-sieve-for-pypy


###To build: python3 setup.py build_ext --inplace
###To run: python3 run_qs.py -keysize 30 -base 100 -debug 1 -lin_size 100_000 -quad_size 100_000


##Notes: Work in progress, this still needs p-adic lifting for the find_similar_smooth portion to work.

##WORK IN PROGRESS....

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
build_workers=1
keysize=30           #Generate a random modulus of specified bit length
workers=1 #max amount of parallel processes to use
quad_co_per_worker=1 #Amount of quadratic coefficients to check. Keep as small as possible.
base=100
lin_sieve_size=100
lin_sieve_size2=100_000
quad_sieve_size=100_000
g_debug=0 #0 = No debug, 1 = Debug, 2 = A lot of debug
g_lift_lim=0.5
thresvar=40  ##Log value base 2 for when to check smooths with trial factorization. Eventually when we fix all the bugs we should be able to furhter lower this.
thresvar2=20
matrix_mul=1.5 ##1.0 = square.. increase to overshoot min smooths, for this version, gathering double the amount seems best, since we don't concentrate all our smooths at a few quadratic coefficients as other versions do.
lp_multiplier=2
min_prime=1
g_enable_custom_factors=0
g_p=107
g_q=41
mod_mul=0.5



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
        
def QS(n,factor_list,sm,xlist,flist):
    g_max_smooths=round((base+1)*matrix_mul)
    if len(sm) > g_max_smooths: 
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


        sqr1=prod_left**2%N 
        sqr2=prod_right%N
        if sqrt_right**2 != prod_right:
            print("something fucked up")
        if sqr1 != sqr2:
            print("ERROR ERROR")
        ###End debug shit#########
        prod_left = prod_left % N
        sqrt_right = sqrt_right % N
        factor_candidate = gcd(N, abs(sqrt_right-prod_left))
     #   print(factor_candidate)
        if factor_candidate not in (1, N):
            other_factor = N // factor_candidate
            return factor_candidate, other_factor

    return 0, 0

def solve_bits(matrix):
    n=round((base+1)*matrix_mul)
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

def build_matrix(factor_base, smooth_nums, factors):
    fb_len = len(factor_base)
    fb_map = {val: i for i, val in enumerate(factor_base)}
    ind=1
    #factor_base.insert(0, -1)
    M2=[0]*(round((base+1))+1)
    for i in range(len(smooth_nums)):
        for fac in factors[i]:
            idx = fb_map[fac]
            M2[idx] |= ind
        ind = ind + ind
    return M2

@cython.profile(False)
def launch(n,primeslist1):
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
    indexmap=create_hmap2indexmap(complete_hmap,primeslist1)
    
    z=0
    print("[*]Launching attack with "+str(workers)+" workers\n")
    find_comb(n,complete_hmap,primeslist1,indexmap)

    return 


def equation(y,x,n,mod,z,z2):
   # print("z: "+str(z)+" x: "+str(x)+" y: "+str(y)+" z2: "+str(z2)+" n: "+str(n))
    rem=z*(x**2)-y*x+n*z2
    rem2=rem%mod
    return rem2,rem 

def squareRootExists(n,p,b,a):
    b=b%p
    c=n%p
    bdiv = (b*inverse(2,p))%p
    alpha = (pow_mod(bdiv,2,p)-c)%p
    if alpha == 0:
        return 1
    
    if jacobi(alpha,p)==1:
        return 1
    return 0

def inverse(a, m):
    if gcd(a, m) != 1:
        return None
    u1,u2,u3 = 1,0,a
    v1,v2,v3 = 0,1,m
    while v3 != 0:
        q = u3//v3
        v1,v2,v3,u1,u2,u3=(u1-q*v1),(u2-q*v2),(u3-q*v3),v1,v2,v3
    return u1%m

def pow_mod(base, exponent, modulus):
    return pow(base,exponent,modulus)  



@cython.profile(False)
cdef tonelli(long long n, long long p):  # tonelli-shanks to solve modular square root: x^2 = n (mod p)
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
   # i=0
    offset=0
    ret=[]
    while 1:
        root=r+offset
        if root > prime**exp:
            break

        rem,rem2=equation2(0,root,n,prime**exp,z,z2)
        #if rem!=0:
          #  print("error erorr")
       # if g_debug > 1:
           # print("[Debug 2]Prime**exp: "+str(prime**exp)+" root: "+str(root)+" rem2: "+str(rem2)+" rem: "+str(rem))
        if rem ==0:
            co2=(formal_deriv(0,root,z))%(prime**exp)
           # print("co2: "+str(co2)+" orig co: "+str(co)+" prime: "+str(prime)+" z: "+str(z))
           # rem,rem2=equation(co2,root,n,prime**exp,z,z2) 
          #  print("rem: ",rem)
           # if rem == 0:
            ret.extend([co2,root])
                #TO DO: Can we just break here?
          #  else:
              #  print("error error error")
        offset+=prime**(exp-1)
    return ret

def lift_a(prime,n,co,z,totalmod,lin):






    z2=1
    k=0
    ret=[]
    cos=[]
    step_size=[]
    new=[]
    r=get_root(n,prime,co%prime,z) 
    if r==-1:
        return 0



   # rem,rem2=equation2(0,r,n,prime,z,z2)
   # if rem != 0:
      #  print("error") ##If this never triggers we can delete this

    exp=2
    ret=[co,r]
   # print("ret: ",ret)
    cos.append(co)
    step_size.append(0)
    while 1 :
        ret=lift(exp,ret[0],ret[1],n,z,z2,prime)
     #   print("dist1: ",ret[0])
        dist1=solve_lin_con(totalmod,ret[0]-lin,prime**exp)
        #print("dist1: ",dist1)
        test=(lin+totalmod*dist1)**2-n*4*z
        if test % (prime**exp) != 0:
            dist1=-1
           # print("lin: ",lin)
           # print("prime: ",prime)
           # print("co: ",co)
           # print("z: ",z)
           # print("totalmod: ",totalmod)
           # print("error dist1***************************")
           # time.sleep(100000)
        co2=(prime**exp)-ret[0]
        cos.append([ret[0],co2])
        #print("dist2: ",co2)
        dist2=solve_lin_con(totalmod,co2-lin,prime**exp) 
      #  print("dist2: ",dist2)
        test=(lin+totalmod*dist2)**2-n*4*z
        if test % (prime**exp) != 0:
            dist2=-1
        if dist1 > lin_sieve_size2 or dist2 > lin_sieve_size2: ##I don't know man
          #  print("dist1: ",dist1)
           # print("dist2: ",dist2)
          #  print("exp: ",exp)
            break
        step_size.append([dist1,dist2])
        if len(ret) > 2:
            print("SHOULDNT HAPPEN")
        exp+=1
    return cos, step_size   

@cython.profile(False)
def solve_roots(prime,n):
    iN=1      
    try:
        #print(prime)
        size=prime*2
        if size > (quad_sieve_size*2)+1:
            size=(quad_sieve_size*2)+1
        temp_hmap = array.array('Q',[0]*size) ##Got to make sure the allocation size doesn't overflow.... 
        temp_hmap[0]=1
        modi=modinv(4*n,prime)
        
        while iN < prime:
            new_square=(iN*4*n)%prime
            test=jacobi(new_square,prime)
            if test ==1:
                root=tonelli(new_square,prime)
                s=(root**2*modi)%prime
                if s > quad_sieve_size:
                    iN+=1
                    continue
                if root > prime // 2:
                    root=(prime-root)%prime
                end=temp_hmap[0]
                temp_hmap[end]=s
              #  root=lift_a(prime,n,root,s)
                temp_hmap[end+1]=root
                temp_hmap[0]=temp_hmap[0]+2
                
               # time.msleep(100000)
      #  if test == 0: #note: Do we want to use these or not?
           # end=temp_hmap[0]
           # temp_hmap[end]=0
           # temp_hmap[end+1]=0
          #  temp_hmap[0]=temp_hmap[0]+2
            iN+=1   
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

  #  print("hmap: ",hmap)
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





#@cython.boundscheck(False)
#@cython.wraparound(False)
cdef factorise_fast(value,unsigned long long [::1] factor_base):
  
    factors = set()
    if value < 0:
        factors ^= {-1}
        value = -value
    while value % 2 == 0:
        factors ^= {2}
        value //= 2
    #cdef int factor 
    cdef Py_ssize_t length=factor_base[0]#len(factor_base)#factor_base[0]
    cdef Py_ssize_t i=1
    while i < length:
        factor=factor_base[i]
        while value % factor == 0:
            factors ^= {factor}
            value //= factor
        i+=1
    return factors, value


#@cython.boundscheck(False)
#@cython.wraparound(False)
cdef void miniloop_non_simd(unsigned long long dist1,int [::1] temp,unsigned long long prime,int log,Py_ssize_t size):
    while dist1 < size :
        temp[dist1]+=log
        dist1+=prime
    return

#@cython.boundscheck(False)
#@cython.wraparound(False)
cdef construct_interval_2(quad_co,lin_co,cmod,unsigned long long [::1] primeslist,hmap,n,indexmap, int [::1] temp,int [::1] logmap,size,unsigned long long [::1] lprimes,cfact):
    ##Still massively bottlenecking. Its often faster to work with a smaller factor base... that doesn't make much sense.
    cdef Py_ssize_t i=1
    cdef Py_ssize_t k
    cdef Py_ssize_t length = lprimes[0]
    local_primes=array.array('Q',[0]*base)
    local_primes[0]+=1
    for fac in cfact:
        end=local_primes[0]
        local_primes[end]=fac   
        local_primes[0]+=1

    cdef unsigned long long prime_index
    cdef unsigned long long  dist1,dist2
    cdef unsigned long long log
    while i < length:
        prime_index=lprimes[i]
        prime=primeslist[prime_index]
        if prime<min_prime: 
            i+=1
            continue
        end=local_primes[0]
        local_primes[end]=prime
        local_primes[0]+=1
        lcmod=cmod%prime

        log=logmap[prime_index]
        modi=modinv(lcmod,prime)

        k=indexmap[prime_index][quad_co%prime]+1
        co2=hmap[prime_index][k]
        res=(co2-lin_co)%prime
        dist1=(res*modi)%prime
        co2=prime-co2
        res=(co2-lin_co)%prime
        dist2=(res*modi)%prime
    
        miniloop_non_simd(dist1,temp,prime,log,size)  ##Question is it better to do dist1 and dist2 in one function call?
        if dist1 != dist2:
             miniloop_non_simd(dist2,temp,prime,log,size)

        #CAN=(lin_co+dist1*cmod)**2-n*4*quad_co
        #if CAN % (cmod*prime)  != 0:

            #print("EROROREROR",prime)
            #print(CAN%cmod)
            #print(CAN)
            #print("lin_co: ",lin_co)
            #print("quad_co: ",quad_co)
            #print("cmod: ",cmod)
            #print("dist1: ",dist1)
        i+=1
    return temp,local_primes


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


#@cython.boundscheck(False)
#@cython.wraparound(False)
cdef sieve_quads(cfact,n,local_mod,quad_interval,quad_interval_index,quad_co):
    found=0
    local_primes=array.array('Q',[0]*base)
    local_primes[0]+=1
    i=quad_co-1
    j=1
    k=0

    length=quad_interval[i][0]
    hit=0
    while j < length:
        if quad_interval[i][j] in cfact:
            hit+=1
        else:
            end=local_primes[0]
            local_primes[end]=quad_interval_index[i][j]  ###To do: Since we're dealing with primes, we can optimize this using a unique set... 
            local_primes[0]+=1

        j+=1
    if hit == len(cfact):
        found=1
    return found,local_primes


cdef get_lin(hmap,indexmap,cfact,local_mod,indexes,quad_co):
    
    all_lin_parts=[]
    j=0
    lin=0
    while j < len(indexes):
        ind=indexes[j]
        prime=cfact[j]
        ind2=indexmap[ind][quad_co%prime]

        r1=hmap[ind][ind2+1]
        aq = local_mod // prime
        invaq = modinv(aq%prime, prime)
        gamma = r1 * invaq % prime
        lin+=aq*gamma
        all_lin_parts.append(aq*gamma)
        j+=1
    lin%=local_mod
    return lin,all_lin_parts

cdef construct_interval(list ret_array,partials,n,primeslist,hmap,gathered_quad_interval,gathered_ql_interval,rstart,rstop,quad_interval,quad_interval_index,threshold_map,indexmap,seen,large_prime_bound,tmul,tnum_list,lprimes_list,tnum_bit_list,logmap):
    grays = get_gray_code(20)
    target_main = array.array('i', [0]*lin_sieve_size*2)
    target_main2 = array.array('i', [0]*lin_sieve_size2*2)
    cdef Py_ssize_t i
    cdef Py_ssize_t j
    close_range = 5
    too_close = 1
    LOWER_BOUND_SIQS=1#400
    UPPER_BOUND_SIQS=4000
    cdef Py_ssize_t size
    primelist=copy.copy(primeslist)
    primelist.insert(0,2) ##To do: remove when we fix lifting for powers of 2
    primelist.insert(0,-1)
    primeslist=array.array('Q',primeslist)
    mod_found=0
    mod2_found=0
    while len(ret_array[0]) < round((base+1)*matrix_mul):
        i=(rstart-1)
        gathered_quad_interval_len=len(gathered_quad_interval)
        while i < 1: #To do: remove or use for multithreading?
            tnum = tnum_list[i]
            tnum_bit=tnum_bit_list[i]
            primesl=lprimes_list[i]
            local_mod,cfact,indexes=generate_modulus(n,primesl,seen,tnum,close_range,too_close,LOWER_BOUND_SIQS,UPPER_BOUND_SIQS,tnum_bit,quad_interval_index[i])
            print("Generating new modulus: ",local_mod)
            if local_mod>1:
                if len(cfact)<2:
                    break

                c=1
                while c < 2:#quad_sieve_size+1:
                    check,lprimes=sieve_quads(cfact,n,local_mod,quad_interval,quad_interval_index,c)  

                    if check == 0:
                        c+=2
                        continue
                    lin,lin_parts=get_lin(hmap,indexmap,cfact,local_mod,indexes,c)  
                    quad=c
                    if quad == 0 or lin == 0:
                        continue  
  
                    poly_ind=0
                    end = 1 << (len(cfact) - 1)
                    while poly_ind < end:
                        target = target_main[:]
                        if poly_ind != 0:
                            v,e=grays[poly_ind] ##Brilliant trick, whoever come up with this.
                            lin=(lin + 2 * e * lin_parts[v])%local_mod
                        if (lin**2-n*4*quad)%local_mod != 0:
                            print("big error")
                            time.sleep(100000)
                        size=lin_sieve_size*2
                        single_interval,local_primes=construct_interval_2(quad,lin,local_mod,primeslist,hmap,n,indexmap,target,logmap,size,lprimes,cfact)

                        mod_found,sim_found=process_interval(ret_array,single_interval,n,quad,lin,partials,large_prime_bound,local_primes,threshold_map[i],local_mod,size,mod_found,quad_interval,hmap,indexmap,quad_interval_index,grays,target_main2)
                        if mod_found+1 > 10:  ##I don't know what the optimal value is here? 
                            print(str(len(ret_array[0])))#+" / "+str(int((base+1)*matrix_mul)))
                            mod2_found+=mod_found
                            mod_found=0
                        if mod2_found+1 > 500 or sim_found == 1:
                            print("Running linear algebra step with #smooths: ",len(ret_array[0]))
                            test=QS(n,primelist,ret_array[0],ret_array[1],ret_array[2])   
                            if test!=0:
                                return 
                            mod2_found=0
                        if len(ret_array[0]) > round((base+1)*matrix_mul):
                            test=QS(n,primelist,ret_array[0],ret_array[1],ret_array[2])  
                            return

                        

                        poly_ind+=1

                    c+=2
            i+=1
    return 


def get_root(n,p,b,a):
    ##Just factor the quadratic I guess.
    ##Need tonelli shanks only if it is irreducable?
    a_inv=inverse(a,p)
    if a_inv == None:
        return -1
    ba=(b*a_inv)%p 
    c=n%p
    ca=(c*a_inv)%p
    bdiv = (ba*inverse(2,p))%p
    return bdiv%p


def solve_lin_con(a,b,m):
    ##ax=b mod m
    g=gcd(a,m)
    a,b,m = a//g,b//g,m//g
    return pow(a,-1,m)*b%m  

cdef construct_similar(quad_co,lin,totalmod,n, int [::1] temp,local,size,indexmap,indexes1,hmap):
  #  print("quad_co: ",quad_co)
  #  print("lin: ",lin)
 #   print("totalmod: ",totalmod)
 #   print("local: ",local)
    logmod=round(math.log2(totalmod))
    length=local[0]
    i=1
  #  print("totalmod: ",totalmod)
    while i < length:
     #   print("local[i] ",local[i])


        ind=indexes1[i-1]
        ind2=indexmap[ind][quad_co%local[i]]
        r1=hmap[ind][ind2+1]

        ret,step_size=lift_a(local[i],n,r1,quad_co,totalmod,lin) ##To do: We only really need step_size, but keep the coefficients for debugging purposes for now
     #   print("ret: ",ret)
      #  print("step_size: "+str(step_size)+" prime: "+str(local[i]))
        j=1
        while j < len(step_size):
            log=round(math.log2(local[i]))
            v=1
            if totalmod%local[i]==0:
                v=0
            if step_size[j][0] != -1:
                miniloop_non_simd(step_size[j][0],temp,local[i]**(j+v),log,size)
            if step_size[j][1] != -1 and step_size[j][0] != step_size[j][1]:
                miniloop_non_simd(step_size[j][1],temp,local[i]**(j+v),log,size)
            j+=1
        #print("ret: ",ret)
       # print("step_size: ",step_size)
        i+=1
    
    
   # print("temp: ",temp)
   # time.sleep(1000000)
    return

def process_similar(temp,lin,quad_co,ret_array,totalmod,local,n):
    found=0
    i=0
    threshold = int(math.log2(math.sqrt(lin_sieve_size2*(n*4*quad_co))) - thresvar2)
  #  print("temp: ",temp)
   # print("threshold: ",threshold)
    while i < len(temp):
        if temp[i]>threshold:
          #  print("assumed: ",temp[i])
          #  print("lin: "+str(lin)+" totalmod: "+str(totalmod)+" i: "+str(i))
            found+=test_similar_smooth(lin+totalmod*i,quad_co,ret_array,totalmod,local,n)

        i+=1

    return found

def test_similar_smooth(lin,quad_co,ret_array,totalmod,local,n):
  #  print("checking")
    poly_val=(lin**2-n*4*(quad_co))
    local_factors,value=factorise_fast(poly_val,local)
   # print("poly_val: "+str(poly_val)+" lin: "+str(lin))
   # print("local_factors: ",local_factors  )
  #  print("value: ",value)
    #TO DO: we need a better modulus... use generate_modulus for now...
    if poly_val%totalmod != 0: ###Debug, delete later
        print("error")
        time.sleep(10000)
    if value != 1:
        return 0
    if poly_val not in ret_array[0]:
      #  print("found")
        print("Found a similar smooth",local_factors)
        ret_array[0].append(poly_val)
        ret_array[1].append(lin)
        ret_array[2].append(local_factors)
        return 1
    else:
        print("already found")
    return 0

def find_similar_smooth_2(quad_co,primes,n,hmap,indexmap,indexes1,local,ret_array,grays,target_main,size):
    ##To do: We shouldn't have to calculate cfact and totalmod like this... and we definitely should't do list conversions and sorting. this is straight up lazy.
  #  print("primes: ",primes)
    tnum = int(((2*n*4*quad_co)**mod_mul) / (lin_sieve_size2//2))
    totalmod,cfact,indexes=generate_modulus_fast(n,primes,tnum,10,1,1,4000,bitlen(tnum),indexes1)
  #  print("local_mod: ",totalmod)
   # print("indexes: ",indexes)
  #  print("cfact: ",cfact)
    found=0
    lin,lin_parts=get_lin(hmap,indexmap,cfact,totalmod,indexes,quad_co) 
    poly_ind=0
    end = 1 << (len(cfact) - 1)
    while poly_ind < end:
        if poly_ind != 0:
            v,e=grays[poly_ind] 
            lin=(lin + 2 * e * lin_parts[v])%totalmod
        poly_ind+=1
        temp = target_main[:]
      #  print("lin: "+str(lin)+" polyind: "+str(poly_ind)+" local_mod: "+str(totalmod)+" cfact: "+str(cfact))
        construct_similar(quad_co,lin,totalmod,n,temp,local,size,indexmap,indexes1,hmap)
      #  print("temp: ",temp)
      #  time.sleep(10000)
        found+=process_similar(temp,lin,quad_co,ret_array,totalmod,local,n)
   # time.sleep(10000)
    return found

def find_similar_smooth(local_factors1,lin_co,quad_co,quad_interval,n,hmap,indexmap,quad_interval_index,ret_array,grays,target_main):
    local_factors=copy.deepcopy(local_factors1) ###Fix all this shit
    local_factors = list(local_factors)
    local_factors.sort()
    if -1 in local_factors: ##To do: I don't know, this is straight up sloppy
        local_factors.pop(0)
    print("local_factors: ",local_factors)
    size=lin_sieve_size2*2
    ##To do: quad_interval should become a hashmap that we index by prime. That should be faster..
    sim_found=0
    i=0
    while i < len(quad_interval):
        if i+1 == quad_co:
            i+=2
            continue
        length=quad_interval[i][0]
        local=array.array('Q',[0]*base)
        local[0]+=1
        indexes=[]



        j=1
        hit=0


        while j < length:
            if quad_interval[i][j] in local_factors:
                hit+=1
                indexes.append(quad_interval_index[i][j])
                end=local[0]
                local[end]=quad_interval[i][j]
                local[0]+=1
            j+=1

        if hit == len(local_factors) and len(local_factors) > 1:
            
            sim_found+=find_similar_smooth_2(i+1,local_factors,n,hmap,indexmap,indexes,local,ret_array,grays,target_main,size)
            if sim_found > len(local_factors)+1: ##Note: overshooting a little to account for -1 and 2
                print("sim_found: ",sim_found)
                return 1
        i+=2

    return 0

#@cython.boundscheck(False)
#@cython.wraparound(False)
cdef process_interval(ret_array,int [::1] interval,n,quad_co,lin_co,partials, large_prime_bound,unsigned long long [::1] local_primes,int threshold,cmod,Py_ssize_t size,mod_found,quad_interval,hmap,indexmap,quad_interval_index,grays,target_main):
    threshold = int(math.log2((lin_sieve_size//2)*math.sqrt(n*4*quad_co)) - thresvar)
    sim_found=0
    cdef Py_ssize_t j=0
    while j < size:
        if interval[j] > threshold:
            co=abs(lin_co+cmod*j)
            poly_val=(co**2-n*4*(quad_co))
            local_factors, value = factorise_fast(poly_val,local_primes)

            if value != 1:
                if value < large_prime_bound:
                    if value in partials:
                        rel, lf, pv = partials[value]
                        if rel == co:
                            j+=1
                            continue
                        co *= rel
                        local_factors ^= lf
                        poly_val *= pv
                    else:
                        partials[value] = (co, local_factors, poly_val)
                        j+=1
                        continue
                else:
                    j+=1
                    continue
            if poly_val not in ret_array[0]:
                if len(local_factors)==0:
                    sqrt_right = math.isqrt(poly_val)
                    if sqrt_right**2 != poly_val:
                        print("WTF")
                        time.sleep(10000)

                    prod_left = co % n
                    sqrt_right = sqrt_right % n
                    factor_candidate = gcd(n, abs(sqrt_right-prod_left))

                    if factor_candidate not in (1, n):
                        other_factor = n // factor_candidate
                        print("Factors: "+str(factor_candidate)+" "+str(other_factor))
                        sys.exit()
                    else:
                        print("nope")
                        j+=1
                        continue


                mod_found+=1
                ret_array[0].append(poly_val)
                ret_array[1].append(co)
                ret_array[2].append(local_factors)
                if len(local_factors)>1:
                    print("Looking for similar smooths, amount needed: "+str(len(local_factors)+2)+" initial smooth factors: "+str(local_factors))
                    if find_similar_smooth(local_factors,co,quad_co,quad_interval,n,hmap,indexmap,quad_interval_index,ret_array,grays,target_main):
                        sim_found=1
                        return mod_found,sim_found

        j+=1

    return mod_found,sim_found


cdef generate_modulus_fast(n,primeslist,tnum,close_range,too_close,LOWER_BOUND_SIQS,UPPER_BOUND_SIQS,tnum_bit,indexes1):
    i=0
    mod=1
    cfact=[]
    indexes=[]
  #  print("primeslist: "+str(primeslist)+" indexes1: "+str(indexes1))
    while i < len(primeslist):
        mod*=primeslist[i]
        cfact.append(primeslist[i])
        indexes.append(indexes1[i])
        if mod > tnum:
         #   print("early mod break")
            break
        i+=1
    return mod,cfact,indexes


#@cython.boundscheck(False)
#@cython.wraparound(False)
cdef generate_modulus(n,primeslist,seen,tnum,close_range,too_close,LOWER_BOUND_SIQS,UPPER_BOUND_SIQS,tnum_bit,quad_interval_index):
    cdef Py_ssize_t counter 
    cdef Py_ssize_t counter2
    cdef Py_ssize_t counter3
    cdef Py_ssize_t counter4
    cdef Py_ssize_t i
    cdef Py_ssize_t j
    cdef Py_ssize_t const_1=1_000
    cdef Py_ssize_t const_2=100_000

    small_B = len(primeslist)
    lower_polypool_index = 2
    upper_polypool_index = small_B - 1
    poly_low_found = False
    
    for i in range(small_B):  ##To do: Can be moved outside mainloop
        if primeslist[i] > LOWER_BOUND_SIQS and not poly_low_found:
            lower_polypool_index = i
            poly_low_found = True
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
            cmod = cmod * potential_a_factor
            cfact.append(potential_a_factor)
            indexes.append(quad_interval_index[randindex+1])
            j = tnum_bit - cmod.bit_length()
            counter2+=1
            if j < too_close:
                cmod = 1
                s = 0
                cfact = []
                indexes=[]
                continue
            elif j < (too_close + close_range):
                break
        a1 = tnum // cmod
        mindiff = 100000000000000000
        randindex = 0
        for i in range(small_B):
            if abs(a1 - primeslist[i]) < mindiff:
                mindiff = abs(a1 - primeslist[i])
                randindex = i

        found_a_factor = False
        counter3=0
        while not found_a_factor and counter3< const_2:
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

        cmod = cmod * potential_a_factor

        cfact.append(potential_a_factor)
        indexes.append(quad_interval_index[randindex+1])

        diff_bits = (tnum - cmod).bit_length()
        if diff_bits < tnum_bit:
            if cmod in seen:
                continue
            else:
                seen.append(cmod)
                return cmod,cfact,indexes
    return 0,0,0

def construct_quad_interval(hmap,primeslist1,rstart,rstop,n):
    ##To do: only odds
    quad_interval_size=rstop
    quad_interval=[]
    quad_interval_index=[]
    threshold_map=[]

    i=rstart
    while i < rstop:
        quad_interval.append(array.array('Q',[0]*base))
        quad_interval_index.append(array.array('i',[0]*base))

        quad_interval[-1][0]=1
        quad_interval_index[-1][0]=1
        threshold = int(math.log2((lin_sieve_size//2)*math.sqrt(n*4*i)) - thresvar) ###To do: move this into the loop so we can get better estimates
        threshold_map.append(threshold)
        i+=1
    i=0
    gathered_quad_interval=[1]*len(quad_interval)
    gathered_ql_interval=[]
    i=0

    while i < len(quad_interval):
        gathered_ql_interval.append([])
        i+=1

    i=0
    while i < len(hmap):
        prime=primeslist1[i]
        length=hmap[i][0]
        j=1
        while j < length:
            x=hmap[i][j]
            start=rstart%prime
            x-=start
            if x < 0:
                x+=prime
            while x < len(quad_interval):
                gathered_quad_interval[x]*=prime
                gathered_ql_interval[x].append(prime)
                gathered_ql_interval[x].append([hmap[i][j+1]])
                end=quad_interval[x][0]
                quad_interval[x][end]=prime
                quad_interval[x][0]+=1
                end=quad_interval_index[x][0]
                quad_interval_index[x][end]=i
                quad_interval_index[x][0]+=1
                x+=prime

            j+=2
        i+=1

    return quad_interval,threshold_map,quad_interval_index,gathered_quad_interval,gathered_ql_interval



cdef list create_hmap2indexmap(hmap,primeslist1):
    cdef list indexmap=[]
    i=0
    while i < len(hmap):
        indexmap.append(array.array('I',[0]*primeslist1[i]))
        length=hmap[i][0]
        j=1
        while j < length:
 
            indexmap[-1][hmap[i][j]]=j

            j+=2
        i+=1

    return indexmap

cdef create_logmap(primeslist1):
    logmap=[]
    i=0
    while i < len(primeslist1):
        logmap.append(round(math.log2(primeslist1[i])))
        i+=1
    logmap=array.array('i',logmap)
    return logmap

def find_comb(n,hmap,primeslist1,indexmap):
   
    logmap=create_logmap(primeslist1)
    ret_array=[[],[],[]]
    seen=[]
    start=default_timer()
    quad_interval,threshold_map,quad_interval_index,gathered_quad_interval,gathered_ql_interval=construct_quad_interval(hmap,primeslist1,1,(quad_sieve_size)+1,n)
    if g_debug > 0:
        duration = default_timer() - start
        print("Constructing quad interval took: "+str(duration))
    partials={}
    gc.enable()
   # residue_map=[]
    tmul =0.9
    sm=[]
    xl=[]
    fac=[]
    large_prime_bound = primeslist1[-1] ** lp_multiplier
    start=default_timer()
    tnum_list=[]
    tnum_bit_list=[]
    lprimes_list=[]           


    i=0
    while i < len(quad_interval):
        primesl=[]
        g=1
        length=quad_interval[i][0]
        while g < length:
            primesl.append(quad_interval[i][g])
            g+=1
        lprimes_list.append(primesl)

        tnum = int(((2*n*4*(i+1))**mod_mul) / (lin_sieve_size//2))
        tnum_list.append(tnum)
        tnum_bit_list.append(bitlen(tnum))
        i+=1

    if g_debug > 0:
        duration = default_timer() - start
        print("Processing quad interval took: "+str(duration))
    construct_interval(ret_array,partials,n,primeslist1,hmap,gathered_quad_interval,gathered_ql_interval,1,quad_sieve_size+1,quad_interval,quad_interval_index,threshold_map,indexmap,seen,large_prime_bound,tmul,tnum_list,lprimes_list,tnum_bit_list,logmap)

    return 0

def get_primes(start,stop):
    return list(sympy.sieve.primerange(start,stop))

@cython.profile(False)
def main(l_keysize,l_workers,l_debug,l_base,l_key,l_lin_sieve_size,l_quad_sieve_size):
    global key,keysize,workers,g_debug,base,key,lin_sieve_size,quad_sieve_size
    key,keysize,workers,g_debug,base,lin_sieve_size,quad_sieve_size=l_key,l_keysize,l_workers,l_debug,l_base,l_lin_sieve_size,l_quad_sieve_size
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
    print("[i]Modulus length: ",bitlen(n))
    count = 0
    num=n
    while num !=0:
        num//=10
        count+=1
    print("[i]Number of digits: ",count)
    print("[i]Gathering prime numbers..")
    primeslist.extend(get_primes(3,1000000))
    i=0
    while i < base:
        primeslist1.append(primeslist[0])
        i+=1
        primeslist.pop(0)   
    launch(n,primeslist1)     
    duration = default_timer() - start
    print("\nFactorization in total took: "+str(duration))




