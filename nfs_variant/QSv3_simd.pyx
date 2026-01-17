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
qbase=5_00
quad_sieve_size=10
g_debug=0 #0 = No debug, 1 = Debug, 2 = A lot of debug
g_lift_lim=0.5
thresvar=300000  ##Log value base 2 for when to check smooths with trial factorization. Eventually when we fix all the bugs we should be able to furhter lower this.
lp_multiplier=2
min_prime=1
g_enable_custom_factors=0
g_p=107
g_q=41
mod_mul=0.25
g_max_exp=10
lin_sieve_size=1000
max_diff=1_000_000

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
        
def QS(n,factor_list,sm,xlist,flist,disc_sr_list,pval_list,pflist):
    g_max_smooths=base*1+2
    if len(sm) > g_max_smooths*10000: 
        del sm[g_max_smooths:]
        del xlist[g_max_smooths:]
        del flist[g_max_smooths:]  
    M2 = build_matrix(factor_list, sm, flist,pflist)
    null_space=solve_bits(M2,factor_list)
    f1,f2=extract_factors(n, sm, xlist, null_space,disc_sr_list,pval_list,pflist)
    if f1 != 0:
        print("[SUCCESS]Factors are: "+str(f1)+" and "+str(f2))
        return f1,f2   
    print("[FAILURE]No factors found")
    return 0,0

def extract_factors(N, relations, roots, null_space,disc_sr_list,pval_list,pflist):
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
                pval*=pval_list[idx]
                disc_sr*=disc_sr_list[idx]
               # sm: "+str(relations[idx])+" adding root: "+str(roots[idx])+" adding pval: "+str(pval_list[idx])+" adding disc_sr: "+str(disc_sr_list[idx]))
            idx += 1
        pval_sr = math.isqrt(pval)
        if pval_sr**2 != pval:
            print("error pval_sr")
        y1_b=disc_sr-pval_sr
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
      #  sqrt_left = sqrt_left % N
      #  sqrt_right = sqrt_right % N
        factor_candidate = gcd(N, abs(sqrt_right+sqrt_left))
        factor_candidate2 = gcd(N, abs(y1_b+sqrt_left))
        fail=0
        if y1_b**2%N != sqrt_left**2%N:
            fail=1
       # print(factor_candidate)
        if factor_candidate not in (1, N):
            print("y1_b: "+str(factor_candidate2)+" y1_b: "+str(y1_b)+" bfail: "+str(fail)+" sm: "+str(sqrt_right)+" root: "+str(sqrt_left))
            other_factor = N // factor_candidate
            return factor_candidate, other_factor

    return 0, 0

def solve_bits(matrix,factor_base):
    n=len(factor_base)*2#base+2
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

def build_matrix(factor_base, smooth_nums, factors,pflist):
    fb_map = {val: i for i, val in enumerate(factor_base)}

    ind=1

    M2=[0]*(len(factor_base)*2)
    for i in range(len(smooth_nums)):
        for fac in factors[i]:
            idx = fb_map[fac]
            M2[idx] |= ind
        ind = ind + ind



    offset=len(factor_base)#(base+2)-1
    ind=1
    for i in range(len(pflist)):
        for fac in pflist[i]:
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
    print("[i]Creating iN datastructure... this can take a while...")
    primeslist1c=copy.deepcopy(primeslist1)
    plists=[]

    hmap,hmap2=create_hashmap(n,primeslist1)

    #complete_hmap=create_hashmap(n,primeslist1)
    duration = default_timer() - start
    print("[i]Creating iN datastructure in total took: "+str(duration))

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

def solve_roots(prime,n):
    hmap_p={}
    hmap_p2={}
    iN=0
    modi=modinv((4*n)%prime,prime)
    while iN < prime:
        ja= jacobi(iN,prime )
        if ja ==1:
            root=tonelli(iN,prime)
            
            if root > prime // 2:
                root=(prime-root)%prime

            s=(root**2*modi)%prime


            if (root**2-4*n*s)%prime!=0:
                print("error123")
            try:
                c=hmap_p[str(s)]
                c.append(root)
            except Exception as e:
                    c=hmap_p[str(s)]=[root]

            try:
                c=hmap_p2[str(root)]
                c.append(s)
            except Exception as e:
                    c=hmap_p2[str(root)]=[s]
              #  z+=1
     #   if ja ==0:  
      #      roots=[0]
      #      z=0
       #     while z < prime:
           #     for root in roots:
              #      res=(root**2-n*4*z)%prime
                #    if res==0:
                    #    try:
                      #      c=hmap_p[str(z)]
                       #     c.append(root)
                      #  except Exception as e:
                         #   c=hmap_p[str(z)]=[root]
               # z+=1
        iN+=1  
    return hmap_p,hmap_p2

def solve_roots2(prime,n):
    hmap_p={}
    iN=0
    while iN < prime:
        ja= jacobi(iN,prime )

        if ja ==1:
            root=tonelli(iN,prime)

              #  time.sleep(10000)
            if root > prime // 2:
                root=-root%prime
            roots=[root]#,(-root)%prime]
            z=1
            while z < prime:
                for root in roots:
                    res=(root**2-n*4*z)%prime
                    if prime == 11 and z == 8 and roots == [2]:
                        print("hit: "+str(roots)+" res: "+str(res)+" jac: "+str(jacobi(res,prime)))
                    if jacobi(res,prime) != -1:
                        if jacobi(res,prime) == 1:
                            y1=tonelli(res,prime)
                            if y1 > prime // 2:
                                y1=-y1%prime
                        else:
                            y1 =0

                        try:
                            c=hmap_p[str(z)]
                            if (root**2-n*4*z)%prime != y1**2%prime:
                                print("fatal error: ")
                            c.append([root,y1])
                        except Exception as e:
                            if (root**2-n*4*z)%prime != y1**2%prime:
                                print("fatal error: ")
                            c=hmap_p[str(z)]=[[root,y1]]
                z+=1
        if ja ==0:  
            roots=[0]
            z=1
            while z < prime:
                for root in roots:
                    res=(root**2-n*4*z)%prime
                    if jacobi(res,prime) != -1:
                        if jacobi(res,prime) == 1:
                            y1=tonelli(res,prime)
                            if y1 > prime // 2:
                                y1=-y1%prime
                        else:
                            y1 =0
                        try:
                            c=hmap_p[str(z)]
                            c.append([root,y1])
                        except Exception as e:
                            c=hmap_p[str(z)]=[[root,y1]]
                z+=1
        iN+=1  
    return hmap_p

def create_hashmap(n,primeslist):
    i=0
    hmap=[]
    hmap2=[]
    while i < len(primeslist):
        hmap_p,hmap_p2=solve_roots(primeslist[i],n)
        hmap.append(hmap_p)
        hmap2.append(hmap_p2)
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
    #print("done")
    return

def solve_lin_con(a,b,m):
    ##ax=b mod m
    #g=gcd(a,m)
    #a,b,m = a//g,b//g,m//g
    return pow(a,-1,m)*b%m

def solve_quad_integers(a,b,c):
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

def grab_co(co,primeslist,hmap):
    i=0
    modulus=1
    z_list=[]
    while i < len(primeslist):
        prime=primeslist[i]
        try:
            test=co%prime
            if test > prime//2:
                test=(-test)%prime
            z=hmap[i][str(test)]
            z_list.append(prime)
            z_list.append(z)
            modulus*=prime
        except Exception as e:
            i+=1
            continue
        i+=1
    print("co: "+str(co)+" modulus: "+str(modulus)+" z_list: "+str(z_list))
    return modulus,z_list

cdef construct_interval(list ret_array,partials,n,primeslist,hmap,hmap2,large_prime_bound,primeslist2,small_primeslist):
 #   i=0
 #   while i < len(hmap):

  #      print("prime: "+str(primeslist[i])+" "+str(hmap[i]))
  #      i+=1
    cdef Py_ssize_t j

    LOWER_BOUND_SIQS=1
    UPPER_BOUND_SIQS=4000
    primelist_f=copy.copy(primeslist)

    loop_list=copy.copy(primeslist)
    loop_list.insert(0,2)
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

    primeslist2.insert(0,2)
    primeslist2.insert(0,-1)


    





    sprimelist=copy.copy(small_primeslist)
    sprimelist.insert(0,2) ##To do: remove when we fix lifting for powers of 2
    sprimelist.insert(0,-1)

    mod_found=0
    mod2_found=0

    many_primes=[]    

    
    M=lin_sieve_size
  #  tnum = int(math.sqrt(2*n) / M)
    tnum = n#int(((2*n)**mod_mul) / M)
    tnum_bit=int(bitlen(tnum)*1)
    close_range =10
    too_close = 5 
    root_list=[]
    poly_list=[]
    flist=[]
    disc_sr_list=[]
    pval_list=[]
    pflist=[]
 
    x1=1
    y0=1
    o=1
    while o < n:
        y=math.isqrt(n*4)+o
        z=1
        seen=[]
        roots=[]
        seen_y=[]
        seen_mod=[]
        seen_z=[]


        smooths=[]
        coefficients=[]
        while z < 2:
            i=0
    
           # y=66
            while i < n:
                fail=0
                x=x1+i
                poly_val=(z*x**2-y*x)%n
                k=((z*x**2-y*x)-poly_val)//n
                k*=-1
                z2=z*k
                if poly_val==0:
                    i+=1
                    continue

                local_factors, value,seen_primes,seen_primes_indexes = factorise_fast(poly_val,primelist_f)
                if value == 1:
                    j=0
                    total_mod=1
                    prev=1
                    while j < len(seen_primes_indexes):
                        pindex=seen_primes_indexes[j]
                        prime=seen_primes[j]
                        if pindex != -1:
                            try:
                                co=hmap[pindex][str(z2%prime)]
                            except Exception as e:
                                fail=1
                                j+=1
                                continue
                            if co[0]**2%prime == y**2%prime and (z*x*2)%prime == y%prime: 
                                if prime !=prev:
                                    total_mod*=prime
                                    prev=prime
                            else:
                                fail=1

                        j+=1

                    k=((z*x**2-y*x)-poly_val)//n
                    k*=-1
                    if ((z*x**2-y*x)-poly_val)%n !=0:
                        print("super fatal error")
                    if z*x**2-y*x+n*k != poly_val:
                        print("fatal error: ",k)
                    if k > 0 and fail ==0 and poly_val%2==0:# and total_mod > 20: 
                        test=abs(poly_val*z)
                        disc1_squared=y**2-4*(n*k*z-(poly_val*z))
                        disc2_squared=y**2-4*(n*k*z)
                        if disc1_squared < 0:
                            print("should never happen")
                        disc1=math.isqrt(disc1_squared)
                        test_sq=disc1_squared%test
                        if test_sq == 0:
                            i+=1
                            continue
                        test_sqr=math.isqrt(test_sq)
                       # if test_sqr**2 ==test_sq:
                       #     print("poly_val: "+str(poly_val)+" disc1_squared: "+str(disc1_squared)+" disc2_squared: " +str(disc2_squared)+" disc1: "+str(disc1)+" test_sq: "+str(test_sq)+" test_sqr: "+str(test_sqr)+" k: "+str(k)+" y: "+str(y)+" total_mod: "+str(total_mod))
                        if test_sqr**2 ==test_sq and disc2_squared>0:
                            testgcd=gcd(test_sqr+y,n)
                            print(" Factor: "+str(testgcd)+" poly_val: "+str(poly_val)+" disc1_squared: "+str(disc1_squared)+" disc2_squared: "+str(disc2_squared)+" disc2_squared%total_mod: " +str(disc2_squared%total_mod)+" disc1: "+str(disc1)+" test_sq: "+str(test_sq)+" test_sq%total_mod: "+str(test_sq%total_mod)+" test_sqr: "+str(test_sqr)+" k: "+str(k)+" y: "+str(y)+" total_mod: "+str(total_mod))
                            if testgcd != 1 and testgcd !=n:
                                sys.exit()
                i+=1
            z+=1
        o+=1
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

def equation2(y,x,n,mod,z,z2):
    rem=z*(x**2)+y*x-n*z2
    rem2=rem%mod
    return rem2,rem

def formal_deriv(y,x,z):
    result=(z*2*x)-(y)
    return result

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

def lift_b(prime,n,co,z,max_prime):






    z2=1
    k=0
    ret=[]
    cos=[]
    step_size=[]
    new=[]
    r=get_root(prime,co%prime,z) 
  #  print("prime: "+str(prime)+" co: "+str(co)+" z: "+str(z)+" r: "+str(r))
    if r==-1:
        return 0



    rem,rem2=equation2(0,r,n,prime,z,z2)
    if rem != 0:
        print("error") ##If this never triggers we can delete this

    exp=2
    ret=[co,r]
   # print("ret: ",ret)
    cos.append([co,(prime-co)%prime])

    while prime**exp < max_prime+1 :
        ret=lift(exp,ret[0],ret[1],n,z,z2,prime)
       # print("ret: ",ret)
        co2=(prime**exp)-ret[0]
      #  if ret[0] > lin_sieve_size2 and co2 > lin_sieve_size2:
       #     break
        
        cos.append([ret[0],co2])

        if len(ret) > 2:
            print("SHOULDNT HAPPEN")
        exp+=1
    return cos[-1]   

def generate_large_square(n,many_primes,valid_quads,valid_quads_factors,sprimelist_f,interval_list,roots,interval_list_pos,partials,large_prime_bound):
    
#prime,n,co,z,max_prime
   # cos=lift_b(3,n,hmap[index][j+1],hmap[index][j],mod_array[c])
   # print("done")
   # time.sleep(10000)
    root_list=[]
    poly_list=[]
    flist=[]
    quadf_list=[]
    i=0
    while i <len(valid_quads):
        quad=valid_quads[i]
        root=roots[i]

        j=0
        while j < lin_sieve_size:
            if interval_list[i][j]>thresvar:
                
                poly_val=quad*(root-j)**2-n
                tot=quad*(root-j)**2

                local_factors, value,seen_primes = factorise_fast(poly_val,sprimelist_f)
                quad_local_factors2=copy.copy(valid_quads_factors[i])
                  #  logged=0
                  #  for p in seen_primes:
                  #      if p != -1 and p != 2:
                  #          logged+=int(math.log2(p))
                 #   print("interval_list[i][j]: "+str(interval_list[i][j])+" seen_primes: "+str(seen_primes)+" logged: "+str(logged)+" quad: "+str(quad))
                test=math.isqrt(value)
                if test**2 != value:
                    if value < large_prime_bound:
                        if value in partials:
                            rel, lf, pv,ql = partials[value]
                            if rel == tot:
                                j+=1
                                continue
                            tot *= rel
                            local_factors ^= lf
                            poly_val *= pv
                            quad_local_factors2 ^=ql
                        else:
                            partials[value] = (tot, local_factors, poly_val,quad_local_factors2)
                            j+=1
                            continue
                    else:
                        j+=1 
                        continue


              
                        
                if tot not in root_list:
                           # print("adding")
                    root_list.append(tot)
                    poly_list.append(poly_val)
                    flist.append(local_factors)
                    quadf_list.append(quad_local_factors2)
                    print("", end=f"[i]Smooths: {len(root_list)} / {small_base*1+2+qbase}\r")
                    if len(root_list)>small_base+2+qbase+2:
                        return root_list,poly_list,flist,quadf_list  
            j+=1

        j=1
        while j < lin_sieve_size:
            if interval_list_pos[i][j]>thresvar:
                poly_val2=quad*(root+j)**2-n
                tot=quad*(root+j)**2
                quad_local_factors2=copy.copy(valid_quads_factors[i])
                local_factors, value,seen_primes = factorise_fast(poly_val2,sprimelist_f)
                   # logged=0
                   # for p in seen_primes:
                    #    if p != -1 and p != 2:
                    #        logged+=int(math.log2(p))
                   # print("interval_list[i][j]: "+str(interval_list_pos[i][j])+" seen_primes: "+str(seen_primes)+" logged: "+str(logged)+" quad: "+str(quad))

                test=math.isqrt(value)
                if test**2 != value:
                    if value < large_prime_bound:
                        if value in partials:
                            rel, lf, pv,ql = partials[value]
                            if rel == tot:
                                j+=1
                                continue
                            tot *= rel
                            local_factors ^= lf
                            poly_val2 *= pv
                            quad_local_factors2 ^=ql
                        else:
                            partials[value] = (tot, local_factors, poly_val2,quad_local_factors2)
                            j+=1
                            continue
                    else:
                        j+=1 
                        continue
                        
                if tot not in root_list:
                    root_list.append(tot)
                    poly_list.append(poly_val2)
                    flist.append(local_factors)
                    quadf_list.append(quad_local_factors2)
                    print("", end=f"[i]Smooths: {len(root_list)} / {small_base*1+2+qbase}\r")
                    if len(root_list)>small_base+2+qbase+2:
                        return root_list,poly_list,flist,quadf_list  
            j+=1
            #print("Smooth candidate #1: "+str(poly_val)+" Smooth candidate #2: "+str(poly_val2)+" quad: "+str(quad+i)+" bitlength smooth can #1: "+str(bitlen(poly_val))+" bitlength smooth can #2: "+str(bitlen(poly_val2)))
        i+=1
    return root_list,poly_list,flist,quadf_list










def find_comb(n,hmap,hmap2,primeslist1,primeslist2,small_primeslist):

    ret_array=[[],[],[],[]]

    

    partials={}


    large_prime_bound = primeslist1[-1] ** lp_multiplier
    construct_interval(ret_array,partials,n,primeslist1,hmap,hmap2,large_prime_bound,primeslist2,small_primeslist)

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