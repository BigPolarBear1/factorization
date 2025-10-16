########################!!!!!!!!!!!!!!!!NOTE THIS VERSION IS BROKEN ON PURPOSE. USE V04 FOR A WORKING VERSION THAT WORKS COMPLETELY IN THE INTEGERS. THIS SIMPLY CONVERTS ONE SIDE OF THE CONGRUENCE TO AN ALGEBRAIC SIDE.. IT NEEDS MORE WORK FOR THE LINEAR ALGEBRA STEP TO WORK, WILL FIX IN V06



###Author: Essbee Vanhoutte
###WORK IN PROGRESS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

###notes: use with python3 polarbearalg_v05.py -key 4387, refer to paper for the math.

####I am working on my own custom number field sieve implementation using reducible quadratic polynomials.

###To do:
#Changed zx-yx+N to the shape of zx-yx instead. This will yield correct results mod p. However, we need to add jacobi symbols now.. will fix in v0.6 but got a massive headache today.

###Changelog:
#v0.5: Changed zx-yx+N to the shape of zx-yx.. so smooth candidates are smaller. But needs more work for the linear algebra step to work again
#v0.4: Changed create_hashmap to create all quadratic residues, not just 0 solutions
#v0.3: Added the linear algebra step
#v0.2: Improved sieving for smooths and making sure both sides of the congruence are smooth.

import random
import sympy
import itertools
import sys
import argparse
import multiprocessing
import time
import copy
from timeit import default_timer
import math
import array


g_z=1
key=0                 #Define a custom modulus to factor
keysize=12            #Generate a random modulus of specified bit length
workers=1    #max amount of parallel processes to use
p_amount=11 #amount of primes in factor base
sieve_interval=10000000
g_enable_custom_factors=0
g_p=107
g_q=41



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

def bitlen(int_type):
    length=0
    while(int_type):
        int_type>>=1
        length+=1
    return length   

def gcd(a,b): # Euclid's algorithm
    if b == 0:
        return a
    elif a >= b:
        return gcd(b,a % b)
    else:
        return gcd(b,a)

def launch(n,primeslist1):
    manager=multiprocessing.Manager()
    return_dict=manager.dict()
    jobs=[]
    procnum=0
    start= default_timer()
    print("[i]Creating iN datastructure... this can take a while...")
    hmap,hmap2=create_hashmap(primeslist1,n)
    duration = default_timer() - start
    print("[i]Creating iN datastructure in total took: "+str(duration))
    z=0
    print("[*]Launching attack with "+str(workers)+" workers\n")
    part=(sieve_interval+1)//workers
    rstart=2#round(n**0.5)
    rstop=part#+round(n**0.5)
    if rstart == rstop:
        rstop+=1
    while z < workers:
        p=multiprocessing.Process(target=find_comb, args=(primeslist1,n,procnum,return_dict,rstart,rstop,hmap,hmap2))
        rstart+=part  
        rstop+=part  
        jobs.append(p)
        p.start()
        procnum+=1
        z+=1            
    
    for proc in jobs:
        proc.join(timeout=0)        

    start=default_timer()

    while 1:
        time.sleep(1)
        z=0
        balive=0
        while z < len(jobs):
            if jobs[z].is_alive():
                balive=1
            z+=1
        check=return_dict.values()
        for item in check:
            if len(item)>0:
                factor1=item[0]
                factor2=n//item[0]
                if factor1*factor2 != n:
                    print("some error happened")
                print("\n[i]Factors of " +str(n)+" are: "+str(factor1)+" and "+str(factor2))
                for proc in jobs:
                    proc.terminate()
                return 0
        if balive == 0:
            print("[i]All procs exited")
            return 0    
    return 

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



def get_root(n,p,b,a):
    a_inv=inverse(a,p)
    if a_inv == None:
        return -1
    ba=(b*a_inv)%p 
    c=n%p
    ca=(c*a_inv)%p
    bdiv = (ba*inverse(2,p))%p
    return bdiv%p



def solve_roots(prime,n):
    hmap_p={}
    hmap_p2={}
    iN=0
    while iN < prime:
        ja= jacobi(iN,prime )
        if ja ==1:
            root=tonelli(iN,prime)
            roots=[root,(-root)%prime]
            z=0
            while z < prime:
                for root in roots:
                    res=(root**2-n*4*z)%prime
                    if jacobi(res,prime) != -1:
                        try:
                            c=hmap_p[str(z)]
                            c.append(root)
                        except Exception as e:
                            c=hmap_p[str(z)]=[root]

                        try:
                            c=hmap_p2[str(root)]
                            c.append(z)
                        except Exception as e:
                            c=hmap_p2[str(root)]=[z]

                z+=1
        if ja ==0:  
            roots=[0]
            z=0
            while z < prime:
                for root in roots:
                    res=(root**2-n*4*z)%prime
                    if jacobi(res,prime) != -1:
                        try:
                            c=hmap_p[str(z)]
                            c.append(root)
                        except Exception as e:
                            c=hmap_p[str(z)]=[root]

                        try:
                            c=hmap_p2[str(root)]
                            c.append(z)
                        except Exception as e:
                            c=hmap_p2[str(root)]=[z]
                z+=1
        iN+=1  
    return hmap_p,hmap_p2

def create_hashmap(lists,n):
    i=0
    hmap=[]#index by z
    hmap2=[]#index by y
    while i < len(lists):
        hmap_p,hmap_p2=solve_roots(lists[i],n)
        hmap.append(hmap_p)
        hmap2.append(hmap_p2)
        print(str(lists[i])+" "+str(hmap_p))
        i+=1
    return hmap,hmap2

def isqrt(n): # Newton's method, returns exact int for large squares
    x = n
    y = (x + 1) // 2
    while y < x:
        x = y
        y = (x + n // x) // 2

    return x

def jacobi(a, n):
    #assert(n > a > 0 and n%2 == 1)
    t=1
    while a !=0:
        while a%2==0:
            a /=2
            r=n%8
            if r == 3 or r == 5:
                t = -t
                #return -1
        a, n = n, a
        if a % 4 == n % 4 == 3:
            t = -t
        #   return -1
        a %= n
    if n == 1:
        return t
    else:
        return 0    

def tonelli(n, p):  # tonelli-shanks to solve modular square root: x^2 = n (mod p)
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
            invaq = inverse(aq%prime, prime)
            gamma = r1 * invaq % prime
            new_list[-1].append(aq*gamma)
            k+=1
        i+=2
    return new_list

def formal_deriv(y,x,z):
    result=(z*2*x)-(y)
    return result

def enumerated_product(*args):
    yield from itertools.product(*(range(len(x)) for x in args))


def factorise_fast(value,factor_base):
  
    factors = set()
    if value < 0:
        factors ^= {-1}
        value = -value
    while value % 2 == 0:
        factors ^= {2}
        value //= 2
    length=factor_base[0]
    i=1
    while i < length:
        factor=factor_base[i]
        while value % factor == 0:
            factors ^= {factor}
            value //= factor
        i+=1
    return factors, value

def gen_comb(collected,mod,z,n,factor_base,ret_array):
    enum=[]
    i=0
    while i < len(collected):
        enum.append(collected[i+1])
        i+=2

    for idx in enumerated_product(*enum):
        y0_b=0
        #quad=quad_co
        l=0
        while l < len(idx):
            y0_b+=enum[l][idx[l]]
            l+=1
        y0_b%=mod


            
        y0=y0_b
        if y0 > (mod//2):
            continue
        y1_squared=(y0**2-n*4)%mod ##Should be able to just look this up...
        if jacobi(y1_squared,mod) == 1:
            y1=tonelli(y1_squared,mod)
        else: ###Can use case where its 0 also... fix later
            continue
        x=get_root(n,mod,y0,z) 
        y0d=formal_deriv(y0,x,z)
        
        #disc_y1=y1**2+n*4*z
        eq1_full=(z*x**2+x*y0d)
        eq1=eq1_full%mod
        x2=get_root(n,mod,y1,z) 
        y1d=formal_deriv(y1,x2,z)
        eq1_b_full=(z*x2**2+x2*y1d)
        if eq1_full == 0 or eq1_b_full == 0:
            i+=1
            continue
        factors,rem=factorise_fast(eq1_full,factor_base)
        factors2,rem2=factorise_fast(eq1_b_full,factor_base)#//mod)%mod,factor_base)
        if rem == 1 and rem2 == 1:
            ret_array[0].append(eq1_full)
            ret_array[1].append(y0)
            ret_array[2].append(factors)
            ret_array[3].append(z)
            ret_array[4].append(mod)
            ret_array[5].append(x)
            ret_array[6].append(eq1_b_full)#//mod)%mod)
            ret_array[7].append(factors2)
        eq1_b=eq1_b_full%mod
        #if y0 == 148%mod:
          # print("mod: "+str(mod)+" y0: "+str(y0)+" y1: "+str(y1)+" z: "+str(z)+" x: "+str(x)+" x2: "+str(x2)+" eq1_full: "+str(eq1_full)+" eq1_b_full: "+str(eq1_b_full))

    return

def equation2(y,x,n,mod,z,z2):
    rem=z*(x**2)+y*x-n*z2
    rem2=rem%mod
    return rem2,rem
    
def lift(exp,co,r,n,z,z2,prime):
   # i=0
    offset=0
    ret=[]
    while 1:
        root=r+offset
        if root > prime**exp:
            break
        rem,rem2=equation2(0,root,n,prime**exp,z,z2)
        if rem ==0:
            co2=(formal_deriv(0,root,z))%(prime**exp)
            ret.extend([co2,root])
        offset+=prime**(exp-1)
    return ret

def lift_b(prime,n,co,z,max_exp):
    z2=1
    k=0
    ret=[]
    cos=[]
    step_size=[]
    new=[]
    r=get_root(n,prime,co%prime,z) 
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
    return cos   

def extract_factors(N, relations, null_space,relations2,factors1,factors2):

    n = len(relations)
    
    for vector in null_space:
        prod_left = 1
        prod_right = 1
        bits=[]
        for idx in range(len(relations)):
            bit = vector & 1
            vector = vector >> 1
            bits.append(bit)
          #  print("vector: ",bit)
            if bit == 1:
                #prod_left *= roots[idx]
                prod_right *= relations[idx]
                prod_left *= relations2[idx]
                print("adding to right: "+str(relations[idx]))#+" factors1: "+str(factors1[idx]))
                print("adding to left: "+str(relations2[idx])+" factors2: "+str(factors2[idx]))
            idx += 1


        #print("bits: ",bits)
        sqrt_right = math.isqrt(prod_right)
        sqrt_left = math.isqrt(abs(prod_left))
        print("sqrt_right: "+str(sqrt_right%N)+" sqrt_left: "+str(sqrt_left%N))
 


        if sqrt_right**2%N != sqrt_left**2%N:
            print("something fucked up")

        ###End debug shit#########
        sqrt_left = sqrt_left % N
        sqrt_right = sqrt_right % N
        factor_candidate = gcd(N, abs(sqrt_right-sqrt_left))
        if factor_candidate not in (1, N):
            other_factor = N // factor_candidate
            return factor_candidate, other_factor

    return 0, 0


def solve_bits(matrix):
    n=p_amount*2+4
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

def build_matrix(factor_base, smooth_nums, factors,smooth_nums2,factors2):
    print("smooth_nums: ",smooth_nums)
    print("factors: ",factors)
    print("smooth_nums2: ",smooth_nums2)
    print("factors2: ",factors2)
    fb_len = len(factor_base)
    fb_map = {val: i for i, val in enumerate(factor_base)}
    ind=1
    #factor_base.insert(0, -1)
    M2=[0]*(p_amount*2+4)
    for i in range(len(smooth_nums)):
        for fac in factors[i]:
            idx = fb_map[fac]
            M2[idx] |= ind
        ind = ind + ind
    offset=p_amount+2
    ind=1
    for i in range(len(smooth_nums2)):
        for fac in factors2[i]:
            idx = fb_map[fac]
            M2[idx+offset] |= ind
        ind = ind + ind
  #  for vector in M2:
       # bits=[]
        #for idx in smooth_nums:
       # print("row: ",bin(row))

          #  bit = vector & 1
          #  vector = vector >> 1
          #  bits.append(bit)
      #  print(bits)
    return M2

def QS(n,factor_list,sm,xlist,flist,sm2,flist2):
    factor_list.insert(0,2)
    factor_list.insert(0,-1)
    g_max_smooths=p_amount*2+4
    if len(sm) > g_max_smooths: 
        del sm[g_max_smooths:]
        del xlist[g_max_smooths:]
        del flist[g_max_smooths:] 
        del sm2[g_max_smooths:]
        del flist2[g_max_smooths:]   
    M2 = build_matrix(factor_list, sm, flist,sm2,flist2)
    null_space=solve_bits(M2)
    f1,f2=extract_factors(n, sm, null_space,sm2,flist,flist2)
    if f1 != 0:
        print("[SUCCESS]Factors are: "+str(f1)+" and "+str(f2))
        return f1,f2   
    print("[FAILURE]No factors found")
    return 0
                
def lift_all(old,z,n,max_exp):
    
    new=[]
    i=0
    while i < len(old):
        j=0
        new.append(old[i]**max_exp)
        while j < len(old[i+1]):
            if old[i+1][j] < (old[i]//2+1):

                lifted_co=lift_b(old[i],n,old[i+1][j],z,max_exp)
                new.append(lifted_co)
            j+=1
        i+=2
    return new

def find_comb(lists,n,procnum,return_dict,rstart,rstop,hmap,hmap2):
    ret_array=[[],[],[],[],[],[],[],[]]
    max_exp=1
    z=g_z
    z_max=z+10
    factor_base=copy.deepcopy(lists)
    factor_base.insert(0,len(factor_base)+1)
    factor_base=array.array('Q',factor_base)
    while z < z_max:
        i=0
        while i < len(hmap):
            collected=[]   
            
            mod=1
            #factor_base=[]
            try:
                l=hmap[i][str(z%lists[i])]
                collected.append(lists[i])
                collected.append(l)
                mod*=lists[i]
            except Exception as e:
                i+=1
                continue
            gen_comb(collected,mod**max_exp,z,n,factor_base,ret_array)
            i+=1
        z+=2
    if len(ret_array[0])>0:
        g=0
        while g < len(ret_array[0]):
            y1=formal_deriv(ret_array[1][g],ret_array[5][g],ret_array[3][g])
            print("alg1: "+str(ret_array[0][g])+" y0: "+str(ret_array[1][g])+" y1: "+str(y1)+" factors1: "+str(ret_array[2][g])+" z: "+str(ret_array[3][g])+" mod: "+str(ret_array[4][g])+" root: "+str(ret_array[5][g])+" alg2: "+str(ret_array[6][g])+" factors2: "+str(ret_array[7][g]))
            g+=1

    test=QS(n,lists,ret_array[0],ret_array[1],ret_array[2],ret_array[6],ret_array[7])   
    print("done")
    return 0
     

def get_primes(start,stop):
    return list(sympy.sieve.primerange(start,stop))

def main():
    global key
    global base
    global workers
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
    plist=[]
    print("[i]Modulus length: ",bitlen(n))
    print("[i]Gathering prime numbers..")
    primeslist.extend(get_primes(3,1000000))

    i=0
    while i < p_amount:
        if n%primeslist[i] !=0:
            plist.append(primeslist[i])
        i+=1
   # plist=[11,7]
    launch(n,plist)
    duration = default_timer() - start
    print("\nFactorization in total took: "+str(duration))

def print_banner():
    print("Polar Bear was here       ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀                       ")
    print("⠀         ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀ ⣀⣀⣀⣤⣤⠶⠾⠟⠛⠛⠛⠛⠷⢶⣤⣄⣀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀   ")
    print("⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣀⣤⣴⠶⠾⠛⠛⠛⠛⠉⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⠙⠛⢻⣿⣟ ⠀⠀⠀⠀      ")
    print("⠀⠀⠀⠀⠀⠀⠀⢀⣤⣤⣶⠶⠶⠛⠋⠉⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠙⠳⣦⣄⠀⠀⠀⠀⠀   ")
    print("⠀⠀⠀⠀⠀⣠⡾⠟⠉⢀⣀⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠹⣿⡆⠀⠀⠀   ")
    print("⠀⠀⠀⣠⣾⠟⠀⠀⠀⠈⢉⣿⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⢿⡀⠀⠀   ")
    print("⢀⣠⡾⠋⠀⢾⣧⡀⠀⠀⠈⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣄⠈⣷⠀⠀   ")
    print("⢿⡟⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣀⠀⢹⡆⣿⡆⠀   ")
    print("⠈⢿⣿⣛⣀⣀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠹⣆⣸⠇⣿⡇⠀   ")
    print("⠀⠀⠉⠉⠙⠛⠛⠓⠶⠶⠿⠿⠿⣯⣄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢠⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢠⣿⠟⠀⣿⡇⠀   ")
    print("⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠻⣦⡀⠀⠀⠀⠀⠀⠀⠀⠠⣦⢠⡄⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣿⡞⠁⠀⠀⣿⡇⠀   ")
    print("⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠘⣿⣶⠄⠀⠀⠀⠀⠀⠀⢸⣿⡇⢸⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⣴⠇⣼⠋⠀⠀⠀⠀⣿⡇⠀   ")
    print("⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⡿⣿⣦⠀⠀⠀⠀⠀⠀⠀⣿⣧⣤⣿⡄⠀⠀⠀⠀⠀⠀⠀⠀⣿⣾⠃⠀⠀⠀⠀⠀⣿⠛⠀   ")
    print("⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣾⠀⠘⢿⣦⣀⠀⠀⠀⠀⠀⠸⣇⠀⠉⢻⡄⠀⠀⠀⠀⠀⠀⡘⣿⢿⣄⣠⠀⠀⠀⠀⠸⣧⡀   ")
    print("⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢠⣿⠀⠀⠀⠙⣿⣿⡄⠀⠀⠀⠀⠹⣆⠀⠀⣿⡀⠀⠀⠀⠀⠀⣿⣿⠀⠙⢿⣇⠀⠀⠀⠀⠘⣷   ")
    print("⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣸⡏⠀⠀⢀⣿⡿⠻⢿⣷⣦⠀⠀⠀⠹⠷⣤⣾⡇⠀⠀⠀⠀⣤⣸⡏⠀⠀⠈⢻⣿⠀⠀⠀⠘⢿   ")
    print("⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⣴⠿⠁⠀⠀⢸⡿⠁⠀⠀⠙⢿⣧⠀⠀⠀⠀⠠⣿⠇⠀⠀⠀⠀⣸⣿⠁⠀⠀⢀⣾⠇⠀⠀⠀⠀⣼   ")
    print("⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣠⡾⡁⠀⠀⠀⠀⣸⡇⠀⠀⠀⠀⠈⠿⣷⣤⣴⡶⠛⡋⠀⠀⠀⠀⢀⣿⡟⠀⠀⣴⠟⠁⠀⣀⣀⣀⣠⡿   ")
    print("⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢠⣿⣿⣿⣤⣾⣧⣤⡿⠁⠀⠀⠀⠀⠀⠀⠀⠈⣿⣀⣾⣁⣴⣏⣠⣴⠟⠉⠀⠀⠀⠻⠶⠛⠛⠛⠛⠋⠉⠀   ")
    return

def parse_args():
    global keysize,key,workers,debug,g_z
    parser = argparse.ArgumentParser(description='Factor stuff')
    parser.add_argument('-key',type=int,help='Provide a key instead of generating one') 
    parser.add_argument('-keysize',type=int,help='Generate a key of input size')    
    parser.add_argument('-workers',type=int,help='# of cpu cores to use')
    parser.add_argument('-debug',type=int,help='1 to enable more verbose output')
    parser.add_argument('-z',type=int,help='Z variable')

    args = parser.parse_args()
    if args.keysize != None:    
        keysize = args.keysize
    if args.key != None:    
        key=args.key
    if args.workers != None:  
        workers=args.workers
    if args.debug != None:
        debug=args.debug    
    if args.z != None:
        g_z=args.z 
    return

if __name__ == "__main__":
    parse_args()
    print_banner()
    main()



