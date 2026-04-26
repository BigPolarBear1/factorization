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

def get_derivative(f):
    res = [0]*(len(f)-1)
    for i in range(len(f)-1):
        res[i] = (len(f)-1-i)*f[i]
    return res



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
    primelist_f=copy.copy(primeslist1)
    primelist_f.insert(0,len(primelist_f)+1)
    primelist_f=array.array('q',primelist_f)
    degree=4
  #  threshold=1
    k=1
    hashmap=create_map2(n,primeslist1,k,small_base)
    interval= np.zeros((lin_sieve_size,lin_sieve_size),dtype=np.int16)

    i=0
    while i < len(hashmap):
        prime=primeslist1[i]

        for key, value in hashmap[i].items():
            z=key[0]#hashmap[i][j][0]
            y=key[1]#hashmap[i][j][1]
            exp=value#hashmap[i][j][3]
            ###Should be able to tile it, but its a little bit more complicated with this setup
            interval[z,y]+=round(math.log2(prime**exp))
        i+=1

    threshold=round(math.log2(n))
    
    np.putmask(interval, interval<threshold, 0)
    indexlist=np.nonzero(interval)

    indexlist_x=indexlist[1]
    indexlist_y=indexlist[0]
    ind=0
    length=len(indexlist_x)
    checked=0
    while ind < length:# length:  
        y=int(indexlist_x[ind])
        z=int(indexlist_y[ind])

  
        if z == 0:
            ind+=1
            continue


        if interval[z,y]==0:
            print("Error")

        disc=y**2+4*n*k*z 
        print(str(disc)+" interval[z,y]: "+str(interval[z,y]))
        disc_test=math.isqrt(disc)
        if disc_test**2 == disc:
            gcdtest=math.gcd(disc_test+y,n)
            if gcdtest !=1 and gcdtest !=n:
                print("Factors of "+str(n)+" are: "+str(gcdtest)+" and: "+str(n//gcdtest))
                sys.exit()
        ind+=1
    sys.exit() 
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
            ret=-1
            i=0
            while i < p:
                if (a*i**2+b*i+c)%p ==0:
                    ret=i 
                    break
                i+=1
            if ret == -1:
                return []
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

def power2(poly, f, p, exp):
    if exp == 1: return poly

    tmp = power2(poly, f, p, exp>>1)
    tmp = poly_prod(tmp, tmp)

    if exp&1:
        tmp = poly_prod(tmp, poly)
        return div_poly_mod(tmp, f, p)
    
    else: return div_poly_mod(tmp, f, p)

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

def poly_prod(a, b):
    res = [0]*(max(len(a), len(b))+min(len(a), len(b))-1)

    for i in range(len(a)):
        for j in range(len(b)):
            res[i+j] += a[i]*b[j]

    return res

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

def create_map2(n,primeslist,k,small_base):
    hashmap=[]
    new_map=[]
  #  hashmap.append({})
    max_exp=500
    i=0
    while i < len(primeslist):
        hashmap.append({})
        new_map.append({})
        prime=primeslist[i]
        print("Lifting: "+str(prime))
        if prime == 2:
            ##Skip for now, is a special case of hensel. I need to recheck how it worked.
            i+=1
            continue
        if math.gcd(k,prime)!=1:
            ##This also complicates lifting.. need to take a day to figure this out..
            i+=1
            continue
        coeff=[0,0,(-n*k)%prime]
        ranges = [range(start, prime) for start in coeff[:-1]]
        for combo in itertools.product(*ranges):
            cur=list(combo)+[coeff[-1]]
            if cur[0]==0 and cur[1]==0:
                ##I dont know man
                continue
            roots=solve_quadratic_congruence(cur[0], cur[1],cur[2], prime)
            if len(roots) ==1:
             #   print("prime: "+str(prime)+" poly: "+str(cur)+" roots: "+str(roots))
                solutions=[]
                for r in roots:
                    solutions.append([cur[0],cur[1],r])

                exp=1

                while exp < max_exp:
                    new_solutions=[]
                    verify=[]
                    for sol in solutions:
                        z=sol[0]
                        while z < prime**(exp+1):
                            y=sol[1]
                            while y < prime**(exp+1):
                                root=sol[2]
    

                                root=lift_root(z,0,-n*k, sol[2], prime, exp)       
                                poly=(z*root**2+y*root-n*k)%prime**(exp+1)
                                deriv=(2*z*root+y)%prime**(exp+1)

                                if poly == 0 and deriv ==0 and root < prime**(exp+1):   
                                    if z < lin_sieve_size and y < lin_sieve_size:         
                                        new_solutions.append([z,y,root])
                                     #   hashmap[-1].append([z,y,root,exp])#[tuple([z,y])]=exp+1
                                        hashmap[-1][tuple([z,y])]=exp+1
                                y+=prime**exp
                            z+=prime**exp
                    solutions=new_solutions
                #    print("solutions: "+str(new_solutions))                   
                    if len(solutions) == 0:
                        break
                    exp+=1
        

        for key, value in hashmap[-1].items():
            l=list(key)
            z=int(l[0])
            y=int(l[1])
            if value%2==0:
                zc=z 
                
                while zc < lin_sieve_size:
                    yc=y   
                    while yc < lin_sieve_size:
                        check=(yc**2+4*n*k*zc)//(prime**value)
                        ja=jacobi(check,prime**(value+1))
                        if ja != -1:
                            new_map[-1][tuple([zc,yc])]=value
                        yc+=prime**value
                    zc+=prime**value
        i+=1     

    return new_map
    
def enumerated_product(*args):
    for e in itertools.product(*(range(len(x)) for x in args)):
        yield e

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