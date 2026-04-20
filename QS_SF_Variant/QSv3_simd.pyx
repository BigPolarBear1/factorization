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
threshold=2
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
        
def QS(n,factor_list,sm,flist,x_list):#,jsymbols,testl,primeslist2,disc1_squared_list):#,disc_sr_list,pval_list,pflist):
    g_max_smooths=base+2#+qbase
    if len(sm) > g_max_smooths*10000000: 
        del sm[g_max_smooths:]
       # del xlist[g_max_smooths:]
        del flist[g_max_smooths:]  
    M2 = build_matrix(factor_list, sm, flist)#,pflist)
    null_space=solve_bits(M2,factor_list,len(sm))
    f1,f2=extract_factors(n, sm, null_space,x_list)#,disc_sr_list,pval_list,pflist)
    if f1 != 0:
        print("[SUCCESS]Factors are: "+str(f1)+" and "+str(f2))
        sys.exit()
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
        prod_right=x**2
        sqrt_right = x
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

def build_matrix(factor_base, smooth_nums, factors):#,pflist):
    fb_map = {val: i for i, val in enumerate(factor_base)}

    ind=1

    M2=[0]*(len(factor_base)+2)#+qbase)#+2+qbase)
    for i in range(len(smooth_nums)):
        for fac in factors[i]:
            idx = fb_map[fac]
            M2[idx] |= ind
        ind = ind + ind       
    return M2

def get_root(p,b,a):
    a_inv=modinv((a%p),p)
    if a_inv == None:
        return -1
    ba=(b*a_inv)%p 
    bdiv = (ba*modinv(2,p))%p
    return bdiv%p

@cython.profile(False)
def launch(n,primeslist1,primeslist2,small_primeslist):
    found=0
    primelist=copy.copy(primeslist1)
  #  primelist.insert(0,2) ##To do: remove when we fix lifting for powers of 2
    primelist.insert(0,-1)
    smooth_list=[]
    root_list=[]
    factor_list=[]
    primeslist1c=copy.deepcopy(primeslist1)
    plists=[]
    small_base=8
    k=1
    while k < quad_sieve_size+1:
        print("[i]Lifting polynomial coefficients and roots in Z/p^e")
        hmap,root_hmap=create_map(n,primeslist1,k,small_base)
        #for e in hmap[7]:
           # if e[2]**2%primeslist1[7]==5394**2%primeslist1[7]:
           # if e[1]!=0 and e[0]==1 and e[2]==5394%primeslist1[4]:
            #    print(str(e)+" k: "+str(k)+" primeslist1[4]: "+str(primeslist1[4]))
     #   print(hmap[3])
        print("[i]Building Sieve Interval")
        interval=np.zeros([lin_sieve_size,lin_sieve_size],dtype=np.int16)
        i=0
        while i < len(hmap):
            prime=primeslist1[i]
            j=0
            while j < len(hmap[i]):
              #  if hmap[i][j][1]%prime**(hmap[i][j][0]) == 35%prime**(hmap[i][j][0]) and hmap[i][j][2]%prime**(hmap[i][j][0]) == 590%prime**(hmap[i][j][0]):
              #      print("k: "+str(k)+" prime: "+str(primeslist1[i])+" hmap: "+str(hmap[i][j]))
                interval[hmap[i][j][1]::prime**(hmap[i][j][0]),hmap[i][j][2]::prime**(hmap[i][j][0])]+=round(math.log2(prime))
                j+=1
            i+=1

        found+=construct_interval(n,primeslist1,interval,k,smooth_list,root_list,factor_list,root_hmap)
        print("Smooths #: "+str(len(smooth_list)))
        if found >50:
            print("Performing linear algebra")
            QS(n,primelist,smooth_list,factor_list,root_list)
            found=0
        k+=1
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
            return []
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
            return []
        x1 = (y - b_div_2) % p
        x2 = (p - y - b_div_2) % p
        if x1 == 0:
            return [x2]
        elif x2 == 0:
            return [x1]
        else:
            return [x1,x2]

def lift_root(z,y,n, root, p, exp):
    r =(z*root**2+root*y-n)%p**(exp+1)
    deriv = 2*z*root+y
    inv = modinv(deriv, p)
    t = ((-r//p**exp)*inv)%p
    new_root = (root+t*p**exp)%p**(exp+1)
    return new_root

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
   # while value % 2 == 0:
   #     seen_primes.append(2)
   #     seen_primes_indexes.append(-1)
   # #    factors ^= {2}
   #     value //= 2
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

def solve_lin_con(a,b,m):
    ##ax=b mod m
    #g=gcd(a,m)
    #a,b,m = a//g,b//g,m//g
    return pow(a,-1,m)*b%m  

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

def evaluate(f, x):
    res = 0

    for i in range(len(f)-1):
        res += f[i]
        res *= x

    res += f[-1]

    return res

def lift2(n,k):
    ##To do: This needs to be fixed but it's a special case of hensel's lemma... need to refresh how that worked
    ret=[]
    
    solutions=[]
    z=0
    while z < 2:
        y=0
        while y < 2:
            disc=(y**2+4*n*k*z)%2
            if disc == 0:
                solutions.append([1,z,y])
                ret.append([1,z,y,0])
            y+=1
        z+=1

    exp=1
    while exp < 500:
        new_solutions=[]
        for sol in solutions:
            z=sol[1]
            while z < 2**(exp+1):
                y=sol[2]
                while y < 2**(exp+1):
                    disc=(y**2+4*n*k*z)%2**(exp+1)
                    if disc == 0:
                        if z < lin_sieve_size and y < lin_sieve_size:
                            new_solutions.append([exp+1,z,y])
                            ret.append([exp+1,z,y,0])
                    y+=2**exp
                z+=2**exp
        solutions=new_solutions
        if len(solutions)==0:
            break
        exp+=1
  #  print("ret: "+str(ret))
    return ret

def liftp_zero(k,n,prime):
    ##I dont know.. I need to fix this eventually, this is very poor.
    ret=[]
    solutions=[]
    ret.append([1,0,0,0])
    solutions.append([1,0,0])
    

    exp=1
    while exp < 500:
        new_solutions=[]
        for sol in solutions:
            z=sol[1]
            while z < prime**(exp+1):
                y=sol[2]
                while y < prime**(exp+1):
                    disc=y**2+4*n*k*z
                    if disc%prime**(exp+1) ==0:
                        if z < lin_sieve_size and y < lin_sieve_size:
                            ret.append([exp+1,z,y,0])
                            new_solutions.append([exp+1,z,y])
                    y+=prime**exp
                z+=prime**exp
        solutions=new_solutions
        if len(solutions)==0:
            break
        exp+=1

    return ret

def create_map(n,primeslist,k,small_base):
    max_exp=100000
    hmap=[]
    root_hmap=[]
    t=0
    while t < len(primeslist):
        if t < small_base:
            hmap.append([])
        root_hmap.append([])
        prime=primeslist[t]
        print("Calculating Z/"+str(prime)+"^e")
        if prime == 2:
            hmap[-1].extend(lift2(n,k))
            t+=1
            continue

        coeff=[0,0,(-n*k)%prime]
        ranges = [range(start, prime) for start in coeff[:-1]]
        for combo in itertools.product(*ranges):
            cur=list(combo)+[coeff[-1]]
            if cur[0]==0 and cur[1]==0:
               # hmap[-1].extend(liftp_zero(k,n,prime))
                continue

            roots=solve_quadratic_congruence(cur[0], cur[1],cur[2], prime)
            roots2=solve_quadratic_congruence(k, -cur[1],-n*cur[0], prime)
            if len(roots) != 0:
                #print("roots: "+str(roots)+" roots2: "+str(roots2)+" prime: "+str(prime)+" cur: "+str(cur))

                root_hmap[-1].append([cur[0],cur[1],roots,roots2])

            if t > small_base-1:
                continue
            disc=cur[1]**2+4*n*k*cur[0] 
            if disc%prime!=0:
                continue

            
            if len(roots)!=1:
                print("something went wrong")
   
            solutions=[]
            for r in roots:
                solutions.append([cur[0],cur[1],r])
                hmap[-1].append([1,cur[0],cur[1],r])
           # print("roots: "+str(roots)+" solutions: "+str(solutions))
            exp=1

            while exp < max_exp:
                new_solutions=[]
                for sol in solutions:
                    z=sol[0]
                    while z < prime**(exp+1):
                        y=sol[1]
                        while y < prime**(exp+1):
                            root=sol[2]
    

                            root=lift_root(z,0,-n*k, sol[2], prime, exp)       
                            poly=(z*root**2+y*root-n*k)%prime**(exp+1)
                            #if z%prime**(exp+1) == 53%prime**(exp+1) and y%prime**(exp+1) == 346%prime**(exp+1):
                            #    print("orig root: "+str(sol[2])+" new root: "+str(root)+" deriv: "+str(deriv)+" prime: "+str(prime)+" exp : "+str(exp+1)+" poly: "+str(poly))
                            deriv=(2*z*root+y)%prime**(exp+1)

                            if poly == 0 and deriv == 0 and root < prime**(exp+1):            
                              #  if deriv !=0:
                              #      print("something went wrong deriv")
                              #      sys.exit()                   
                                if z < lin_sieve_size and y < lin_sieve_size:# and root < n:
                                #    if z%prime**(exp+1) == 53%prime**(exp+1) and y%prime**(exp+1) == 346%prime**(exp+1):
                                  #      print("ADDDING orig root: "+str(sol[2])+" new root: "+str(root)+" deriv: "+str(deriv)+" prime: "+str(prime)+" exp : "+str(exp+1))

                                    new_solutions.append([z,y,root])
                                    #if (exp+1)%2 == 0:
                                    hmap[-1].append([exp+1,z,y,root])
                            y+=prime**exp
                        z+=prime**exp
                solutions=new_solutions
                if len(solutions) == 0:
                    break
                exp+=1
    #    print("root_hmap: "+str(root_hmap[-1]))
        t+=1
    return hmap,root_hmap

def create_div_interval(root_hmap,primeslist,r1,r2,mod,lin_co):
    ##TO DO: Next step is to just create an interval here. Easy enough.. then after that.. I'll get rid of the interval and just calculate that divisor, but I want to do it the easy way first.
    i=1
    while i < len(primeslist):
        prime=primeslist[i]
        if mod%prime ==  0 or lin_co%prime ==0:
            i+=1
            continue
        found=0
        j=0
        while j < len(root_hmap[i]):
         #   print("Prime: "+str(prime)+" "+str(root_hmap[i][j]))
            if root_hmap[i][j][1] == lin_co%prime:
                h=0
                while h < len(root_hmap[i][j][-2]):
                    if root_hmap[i][j][-2][h] == r1%prime:
                        if r2%prime in root_hmap[i][j][-1]:
                            found =1
                            break
                    h+=1
            if found == 1:
                break
            j+=1
        if found == 0:
           # print("Returning for prime: "+str(prime))
            return 0
        i+=1

    return 1

def ff_square_root(z,y,k,n,seen_primes,root_hmap,primeslist):
    ###To do: Not finished yet, experimental. Once we find a large finite field.. we then need to look at other primefields to figure out what to divide the root by.
    ###This still needs to be implemented. Will implement shortly
    
    
    prev=0
    r=0
    while r < len(seen_primes):
        partials=[]
        lmod=1
        h=0

        prime=seen_primes[r]

        if (z%prime == 0 and y%prime ==0) or prime == 2:
            r+=1
            continue
          #  d_inv=modinv(d,prime)
        roots=solve_quadratic_congruence(1, y,-n*k*z, prime)
        roots2=solve_quadratic_congruence(k*z, -y,-n, prime)
        root=roots[0]
        root2=roots2[0]
        if root == 0:
            r+=1
            continue
        if (root**2+y*root-n*k*z)%prime !=0:
            print("something seems to have gone wrong..")
            sys.exit()
        if (k*z*root2**2-y*root2-n)%prime!=0:
            print("2 something seems to have gone wrong..")
            sys.exit()
        

        prev=prime
        r+=1
        exp=1
        while r < len(seen_primes) and seen_primes[r] == prev:
            root=lift_root(1,0,-n*k*z, root, prime, exp)   
            root2=lift_root(k*z,0,-n, root2, prime, exp)   
            if (root**2+y*root-n*k*z)%prime**(exp+1)!=0:
                print("fatal: "+str(root)+" prime: "+str(prime)+" exp: "+str(exp+1))
                sys.exit()
            if (k*z*root2**2-y*root2-n)%prime**(exp+1)!=0:
                print("2 fatal: "+str(root)+" prime: "+str(prime)+" exp: "+str(exp+1))
                sys.exit()
            exp+=1
            r+=1
        if prime**exp > round(n**0.6):
            ###To do: Need to calculate d in another larger finite field rather then just bruteforce.. or alternatively create a sieve interval for it.
            ###Almost finished now...
            print("Checking ff: "+str(prime**exp))
            d=1
            while d < prime**exp:
                
                if d%prime ==0:
                    d+=1
                    continue
                d_inv=modinv(d,prime**exp)
                new_root=(root*d_inv)%prime**exp
                new_root2=(root2*d)%prime**exp
                result=create_div_interval(root_hmap,primeslist,new_root,new_root2,prime**exp,y)
                if result == 1:
                    print("D: "+str(d))
                else:
                    d+=1
                    continue
                if (d*new_root**2+y*new_root-n*k*z*d_inv)%prime**exp !=0:
                    print("d: ",d)
                    sys.exit()
                if (k*z*d_inv*new_root2**2-y*new_root2-n*d)%prime**exp !=0:
                    print("d: ",d)
                    sys.exit()
                gcdtest=math.gcd(new_root,n)
                gcdtest2=math.gcd(new_root2,n)
                #disc=y**2+4*n*k*z
                if gcdtest != 1 and gcdtest !=n and new_root//gcdtest == 1:
                    print("[!!!!!FOUND BY TAKING SQUARE ROOT OVER FINITE FIELD (Note: This function isn't finished yet, will update soon)!!!!!]Factors of "+str(n)+" : "+str(gcdtest)+" and "+str(n//gcdtest)+" result: "+str(result)+" new_root: "+str(new_root)+" new_root2: "+str(new_root2)+" prime**exp: "+str(prime**exp))#+" d: "+str(d)+" gcdtest2: "+str(gcdtest2)+" new_root//gcdtest: "+str(new_root//gcdtest)+" prime**exp: "+str(prime**exp)+" root: "+str(new_root)+" root2: "+str(new_root2))
                    sys.exit()
                d+=1

  
cdef construct_interval(n,primeslist,interval,k,smooth_list,root_list,factor_list,root_hmap):
    found=0
    primelist_f=copy.copy(primeslist)
    primelist_f.insert(0,len(primelist_f)+1)
    primelist_f=array.array('q',primelist_f)
    #np.set_printoptions(
    #linewidth=500,   # Max characters per line before wrapping
    #precision=1,     # Decimal places for floats
    #suppress=True,   # Avoid scientific notation for small numbers
    #edgeitems=20,#np.inf # Show full array without truncation
    #threshold=10
    #)  
    np.putmask(interval, interval<threshold, 0)
    indexlist=np.nonzero(interval)
    indexlist_x=indexlist[1]#cp.asnumpy(indexlist[1])
    indexlist_y=indexlist[0]#cp.asnumpy(indexlist[0])
    ind=0
    length=len(indexlist_x)
    checked=0
    while ind < length:# length:  
        y=int(indexlist_x[ind])
        z=int(indexlist_y[ind])

        current = [z,y]
        if z == 0:
            ind+=1
            continue
        disc=y**2+4*n*k*z

        factors, value,seen_primes,seen_primes_indexes=factorise_fast(disc,primelist_f)
        if value != 1:
            ind+=1
            continue

        ff_square_root(z,y,k,n,seen_primes,root_hmap,primeslist)
     
        log=0
        g=0
        for prime in seen_primes:
            if prime < 20:
                log+=round(math.log2(prime))
    #    if log != interval[z,y]:
    #        print("fatal error: "+str(seen_primes)+" log: "+str(log)+" inter: "+str(interval[z,y])+" disc: "+str(disc)+" current: "+str(current)+" k: "+str(k))
    #        sys.exit()

        if factors not in factor_list:
            print("found: "+str(disc)+" seen_primes: "+str(seen_primes)+" interval value: "+str(interval[z,y])+" log: "+str(log)+" coeff: "+str(current))
            smooth_list.append(disc)
            root_list.append(current[1])
            factor_list.append(factors)
            found+=1
        ind+=1
    return found






def get_primes(start,stop):
    return list(sympy.sieve.primerange(start,stop))

@cython.profile(False)
def main(l_keysize,l_workers,l_debug,l_base,l_key,l_lin_sieve_size,l_quad_sieve_size,l_threshold):
    global key,keysize,workers,g_debug,base,key,quad_sieve_size,small_base,lin_sieve_size,qbase,threshold
    key,keysize,workers,g_debug,base,quad_sieve_size,lin_sieve_size,threshold=l_key,l_keysize,l_workers,l_debug,l_base,l_quad_sieve_size,l_lin_sieve_size,l_threshold
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

    primeslist.extend(get_primes(2,20000000))
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