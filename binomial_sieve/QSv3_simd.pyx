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
        print(" polyval sqrt: "+str(sqrt_left%N)+" disc constant sqrt: "+str(sqrt_right%N))#+" zx*zxy: "+str(x))
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

#@cython.boundscheck(False)
#@cython.wraparound(False)
def generate_modulus(n,primeslist,seen,tnum,close_range,too_close,LOWER_BOUND_SIQS,UPPER_BOUND_SIQS,tnum_bit,quad):
    const_1=1_000
    const_2=1_000_000

    small_B = len(primeslist)
    lower_polypool_index = 2
    upper_polypool_index = small_B - 1
    poly_low_found = False
    
    for i in range(small_B):  ##To do: Can be moved outside mainloop
        if primeslist[i] > LOWER_BOUND_SIQS and not poly_low_found:
            lower_polypool_index = i
            poly_low_found = True
            break
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
                if  jacobi((-quad*n)%primeslist[randindex],primeslist[randindex])!=1:
                    counter+=1
                    continue
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
            if abs(a1 - primeslist[i]) < mindiff:
                randindex = i
                mindiff = abs(a1 - primeslist[i])
                
        

        found_a_factor = False
        counter3=0
        while not found_a_factor and counter3< const_2 and randindex <len(primeslist):
            if  jacobi((-quad*n)%primeslist[randindex],primeslist[randindex])!=1:
                randindex += 1
                continue
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


def get_lin(cfact,local_mod,indexes,quad,n):
    all_lin_parts=[]
    j=0
    lin=0
    #print("indexes: ",indexes)
    while j < len(indexes):
        ind=indexes[j]
        prime=cfact[j]
        r1=solve_quadratic_congruence(1, 0, n*quad, prime)[0]
        
       # r1=int(roots2d[quad_co,ind])#hmap[ind][2]
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

def construct_interval(primeslist,f,step_size):
    interval=array.array('i',[0]*lin_sieve_size)
    i=0
    while i < len(primeslist):
        p=primeslist[i]
        fc=copy.deepcopy(f)
        roots=find_roots_poly(fc, p)
        for r in roots:
            s=solve_lin_con(step_size,r,p)
            log=round(math.log2(p))
            while s < lin_sieve_size:
                interval[s]+=log
                s+=p
        i+=1
    return interval

def get_mod_root(mod_item,resmap,primeslist1,x):
    mod=mod_item[0]
    cfact=mod_item[1]
    indexes=mod_item[2]
    pl=[]
    i=0
    while i < len(indexes):
        prime=cfact[i]
        if prime != primeslist1[indexes[i]]:
            print("fatal: "+str(prime)+" "+str(primeslist1[indexes[i]]))
        if prime == 2:
            mod//=2
            i+=1
            continue
    
        try:
            xm=x%prime

            res=resmap[indexes[i]][xm]
            pl.append(prime)
            pl.append([])
            for re in res:
                pl[-1].append(re)
        except Exception as e:
            i+=1
            continue
        i+=1
  #  print("pl: "+str(pl)+" x: "+str(x))
    pl=get_partials(mod,pl)
  #  print("pl: "+str(pl))
    kres=0
    i=0
    while i < len(pl):
        kres+=pl[i+1][0]
        i+=2
    kres%=mod
   # sys.exit()
    return kres,mod

def find_similar(y,primeslist1,primelist_f,n,smooth_list,factor_list,root_list,factor_list2):
    grays = get_gray_code(20)
    seen=[]
    facseen=[]
    close_range=10
    too_close=5
    LOWER_BOUND_SIQS=1
    UPPER_BOUND_SIQS=4000
    tnum=int(((n)**0.5) )#/(lin_sieve_size))
    tnum_bit=bitlen(tnum)
    quad=1
    factors=primeslist1[1:]
    factors.sort()
    modlist=[]
    retry=0
    while len(modlist) < 10:
        new_mod,cfact,indexes=generate_modulus(n,factors,seen,tnum,close_range,too_close,LOWER_BOUND_SIQS,UPPER_BOUND_SIQS,tnum_bit,quad)
        if new_mod ==0:
            retry+=1
            if retry == 100:
                break
            continue
        modlist.append([new_mod,cfact,indexes])
        for fac in cfact:
            if fac not in facseen:
                facseen.append(fac)
    print("modlist: "+str(modlist)) 
    if len(modlist)==0:
        print("something went wrong")
        sys.exit()
    ###To do: So lside will have similar factors, but now we need to do residue sieving so pval also has similar factors that match the original B-smooth
    
    
    print("Binomial term: "+str(y))#+" factors to find: "+str(factors))
    found=0
    degree = 2
    while degree < 5:##To do: Worth going even higher?
        sym_x = sympy.symbols("sym_x")
        formula = (sym_x + y)**degree
        poly=formula.expand().as_poly(sym_x).all_coeffs()
        print("poly: "+str(poly))
        resmap=[]
        thresmod=1
        mod=1
        i=0
        while i < len(factors):
            factor=factors[i]
            ###To do: This is dumb... refactor this later

            resmap.append({})
            if factor == -1 or factor == 2:
                i+=1
                continue
            if factor not in facseen:
                i+=1
                continue
            mod*=factor
           # if factor > 20:
            thresmod*=factor
            k=0
            while k < factor:
                
                polyc=copy.deepcopy(poly)
                polyc[-1]=-n*k
                roots=find_roots_poly(polyc, factor)
                if len(roots)>0:
                    for r in roots:
                        try:
                            res=resmap[-1][r]
                            res.append(k)
                        except Exception as e:
                            resmap[-1][r]=[k]
                  #  print("factor: "+str(factor)+" mod: "+str(mod)+" roots: "+str(roots))
                k+=1

           # print("prime: "+str(factor)+" resmap: "+str(resmap[-1]))
            i+=1

       # print("poly: "+str(poly))
       
        print("[i]Done building root map")
        a = 0
        while a < len(modlist):
            mod_item=modlist[a]
            xc=1
            while xc < 10_000:
          #  print("xc: "+str(xc))
                x=y*xc
                poly[-1]=0
                pval=evaluate(poly,x)
                #print("Pval: "+str(pval)+" poly: "+str(poly))
                if pval == 0:
                    xc+=1
                    continue

                kres,mod=get_mod_root(mod_item,resmap,factors,x)
               # print(kres)
                k_start=(pval)//n

                k_start-=(quad_sieve_size//2)*mod
              #  print("kstart: "+str(k_start))
                r=k_start%mod
                d=(r-kres)%mod
                #diff=abs(k_start)-kres
               # diff=diff//mod
                k_start=k_start-d#kres+diff*mod
                if k_start%mod !=kres:
                    print("super fatal")
                    sys.exit()
              #  print("k_start adjusted: "+str(k_start))
                poly[-1]=-n*k_start
                 #   polyc=copy.deepcopy(poly)
                lside=evaluate(poly,x)
                if k_start%mod != kres:
                    print("oops")
                pval=lside
                if pval%mod !=0:
                    print("fatal error")
                    sys.exit()
                k_interval=array.array('i',[0]*quad_sieve_size)
  #              i=0
  #              while i < len(factors):
  #                  factor=factors[i]
  #                  if mod%factor ==0:
  #                      i+=1
  #                      continue
  #                  try:
  #                      xm=x%factor
  #                      res=resmap[i][xm]
  #                      for re in res:
  #                          s=solve_lin_con(mod,re-k_start,factor)
  #                          if (k_start+s*mod)%factor !=re:
  #                              print("fail")
  #                          j=s
  #                          while j < len(k_interval):
  #                              k_interval[j]+=round(math.log2(factor))
  #                              j+=factor
  #                  except Exception as e:
  #                      i+=1
  #                      continue
  #                  i+=1

           
                i = 0
                while i < len(k_interval):
                    if k_interval[i] > -1:
                        k=k_start+i*mod
                        if k == 0:
                            i+=1
                            continue
                        poly[-1]=-n*k
                 #   polyc=copy.deepcopy(poly)
                        lside=evaluate(poly,x)
                    
                        pval=lside
                        lside+=n*k
                        if lside%y**2 !=0:
                            print("error")
                        lside//=y**2
                        factors1, value,seen_primes,seen_primes_indexes=factorise_fast(pval,primelist_f)
                        factors2, value2,seen_primes2,seen_primes_indexes2=factorise_fast(lside,primelist_f)
                        factors1=list(factors1)
                        factors1.sort()
                      #  seen_primes.sort()
                      #  print("pval: "+str(bitlen(pval))+" lside: "+str(lside)+" poly: "+str(poly)+" x: "+str(x)+" seen_primes: "+str(seen_primes)+" value: "+str(value)+" threshold: "+str(threshold)+" indicated threshold: "+str(k_interval[i])+" value2: "+str(value2)+" mod: "+str(mod))
                        if value == 1 and value2 == 1:
                        #print("hit")
                            if factors1 not in factor_list:# and factors2 not in factor_list2:
                                found+=1
                                smooth_list.append(pval)
                                factor_list.append(factors1)
                                root_list.append(lside*y**2)
                                factor_list2.append(factors2)
                                print("Smooth# "+str(len(smooth_list))+"/"+str(len(primeslist1)*2+10)+" Poly: "+str(poly)+" x: "+str(x)+" pval: "+str(pval)+" pval/mod bits: "+str(bitlen(pval//mod))+" seen_primes: "+str(factors1)+" k: "+str(k)+" seen_primes2: "+str(factors2)+" threshold: "+str(threshold)+" indicated threshold: "+str(k_interval[i])+" k center: "+str(k_start+(quad_sieve_size//2)))
                                if len(smooth_list) > len(primeslist1)*2+10:
                                    #print("returning")
                                    return found
                    i+=1
                xc+=1
            a+=1
        degree+=1
   # sys.exit()

    return found

def binomial_sieve(root_hmap,primeslist1,primelist_f,n,smooth_list,factor_list,root_list,factor_list2,primelist):
    found=0
    degree=1
    x_size=1000
    while degree < 2:
        y_exp=0.3
      #  print("exp: "+str(y_exp))
        y_start=round(n**y_exp) ##We can presieve these aswell
        y=y_start
        while y < y_start+lin_sieve_size:
            found+=find_similar(y,primeslist1,primelist_f,n,smooth_list,factor_list,root_list,factor_list2)
            if found > 50:
                print("Performing linear algebra")
                QS(n,primelist,smooth_list,factor_list,root_list,factor_list2)
                found=0

            y+=1
        print("next degree")
        degree+=2
    return found

def mark_interval(interval,pos,prime,log):
    i=pos
    while i < len(interval):
        interval[i]+=log

        i+=prime

def verify_result(pos,y,degree,n,k,prime):
    ##Debug remove later
    sym_x = sympy.symbols("sym_x")
    formula = (sym_x + y)**degree
    poly=formula.expand().as_poly(sym_x).all_coeffs()
    poly[-1]=-n*k
    pval=evaluate(poly,y*pos)
   # print("pval: "+str(pval)+" y: "+str(y)+" pos: "+str(pos)+" prime: "+str(prime)+" poly: "+str(poly))
    if pval%prime !=0:
        print("verify fail...exiting")
        sys.exit()

def process_interval(interval,y,n,k,degree,primelist_f,smooth_list,factor_list,root_list,factor_list2):
    sym_x = sympy.symbols("sym_x")
    formula = (sym_x + y)**degree
    poly=formula.expand().as_poly(sym_x).all_coeffs()
    poly[-1]=-n*k
    i=0
    found=0
    while i < len(interval):
        if interval[i]>threshold:
            pval=evaluate(poly,y*i)
            lside=pval+n*k
            factors1, value,seen_primes,seen_primes_indexes=factorise_fast(pval,primelist_f)
            if lside%y**2 != 0:
                print("fatal")
                sys.exit()
            lside//=y**2
            factors2, value2,seen_primes2,seen_primes_indexes2=factorise_fast(lside,primelist_f)
            if value == 1 and value2 == 1:
                if factors1 not in factor_list:# and factors2 not in factor_list2:
                    print("Smooth# "+str(len(smooth_list))+" pval: "+str(pval)+" lside: "+str(lside)+" seen_primes: "+str(seen_primes)+" seen_primes2: "+str(factors2)+" index: "+str(i)+" poly: "+str(poly))

                    found+=1
                    smooth_list.append(pval)
                    factor_list.append(factors1)
                    root_list.append(lside*y**2)
                    factor_list2.append(factors2)
                    #print("Smooth# "+str(len(smooth_list))+"/"+str(len(primeslist1)*2+10)+" Poly: "+str(poly)+" x: "+str(x)+" pval: "+str(pval)+" pval/mod bits: "+str(bitlen(pval//mod))+" seen_primes: "+str(factors1)+" k: "+str(k)+" seen_primes2: "+str(factors2)+" threshold: "+str(threshold)+" indicated threshold: "+str(k_interval[i])+" k center: "+str(k_start+(quad_sieve_size//2)))
                    if len(smooth_list) > len(primelist_f)*2+10:
                                    #print("returning")
                        return found
        i+=1
    return found



def sieve_loop(n,root_hmap,primeslist,k,degree,primelist_f,smooth_list,factor_list,root_list,factor_list2,primelist):
    interval=array.array("i",lin_sieve_size*[0])
    found=0
    y=1
    while y < 100:
        i=0
        while i < len(root_hmap):
            prime=primeslist[i]
            if y%prime ==0:
                i+=1
                continue
            log=round(math.log2(prime))
            try:
                res=root_hmap[i][y%prime]
           #     print("Prime: "+str(prime)+" res : "+str(res))
                for item in res:
                    for root in item[-1]:
                        pos=solve_lin_con(y,root,prime)

                       
                        #verify_result(pos,y,degree,n,k,prime) ##Debug remove later
                        mark_interval(interval,pos,prime,log)
            except Exception as e:
              #  print(e)
                i+=1
                continue
            i+=1
        found+=process_interval(interval,y,n,k,degree,primelist_f,smooth_list,factor_list,root_list,factor_list2)

        if found > 50:
            print("Performing linear algebra")
            QS(n,primelist,smooth_list,factor_list,root_list,factor_list2)
            found=0
       # print("interval: "+str(interval))
        y+=1


@cython.profile(False)
def launch(n,primeslist1,primeslist2,small_primeslist):
    found=0
    primelist=copy.copy(primeslist1)
  #  primelist.insert(0,2) ##To do: remove when we fix lifting for powers of 2
    primelist.insert(0,-1)
    smooth_list=[]
    root_list=[]
    factor_list=[]
    factor_list2=[]
    primeslist1c=copy.deepcopy(primeslist1)
    plists=[]
    small_base=8
    primelist_f=copy.copy(primeslist1)
    primelist_f.insert(0,len(primelist_f)+1)
    primelist_f=array.array('q',primelist_f)
    root_hmap=create_map(n,primeslist1,1,2)
    sieve_loop(n,root_hmap,primeslist1,1,2,primelist_f,smooth_list,factor_list,root_list,factor_list2,primelist)
    sys.exit()
    ###To do: implement residue sieving with create_map
    root_hmap=[]
    found+=binomial_sieve(root_hmap,primeslist1,primelist_f,n,smooth_list,factor_list,root_list,factor_list2,primelist)
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
    g=gcd(a,m)
    a,b,m = a//g,b//g,m//g
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

def create_map(n,primeslist,k,degree):
    max_exp=100000

    root_hmap=[]
    t=0
    while t < len(primeslist):
        root_hmap.append({})
        prime=primeslist[t]
        coeff=[(-n*k)%prime]
        d=degree
        d_ind=0
        while d_ind < d-1:
            coeff.insert(0,0)
            d_ind+=1


        ranges = [range(start, prime) for start in coeff[:-1]]
        for combo in itertools.product(*ranges):
            cur=[1]+list(combo)+[coeff[-1]]
           # print(cur)
            if cur[0]==0 and cur[1]==0:
                continue
            curc=copy.deepcopy(cur)##I need to fix find_roots_poly.. so we can remove this..
            roots=find_roots_poly(curc, prime)
            #roots=solve_quadratic_congruence(cur[0], cur[1],cur[2], prime)
            
            ##For the second degree.. since the linear coefficient is twice the binomial term, we divide by 2
            ##I need to generalize this to higher degrees
            inv2=modinv(2,prime)
            binomial=(cur[1]*inv2)%prime
            if len(roots) != 0:
                try:
                    res=root_hmap[-1][binomial]
                    res.append([cur,k,roots])
                except Exception as e:
                    root_hmap[-1][binomial]=[[cur,k,roots]]
               # root_hmap[-1].append([cur,k,roots])
        #print("prime: "+str(prime)+" hmap: "+str(root_hmap[-1]))
        t+=1
    return root_hmap

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

def calc_d(root_hmap,primeslist,r1,r2,mod,lin_co):
    ##TO DO: Next step is to just create an interval here. Easy enough.. then after that.. I'll get rid of the interval and just calculate that divisor, but I want to do it the easy way first.
    d_list=[]

    i=1
    while i < len(primeslist):
        d_list.append([])
        prime=primeslist[i]
        if mod%prime ==  0 or lin_co%prime ==0:
            i+=1
            continue

        
        j=0
        while j < len(root_hmap[i]):
         #   print("Prime: "+str(prime)+" "+str(root_hmap[i][j]))
            if root_hmap[i][j][1] == lin_co%prime:
                temp_r1=[]
                h=0
                while h < len(root_hmap[i][j][-2]):
                    r1_d=root_hmap[i][j][-2][h]


                    temp_r1.append(r1_d)
                    h+=1
                temp_r2=[]
                h=0
                while h < len(root_hmap[i][j][-1]):
                    r2_d=root_hmap[i][j][-1][h]


                    temp_r2.append(r2_d)
                    h+=1
                d_list[-1].append([temp_r1,temp_r2])

            j+=1
       # print("prime: "+str(prime)+" d_list: "+str(d_list[-1]))
        i+=1

    return d_list


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