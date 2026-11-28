#!python
#cython: language_level=3
#cython: profile=False
#cython: overflowcheck=False
###Author: Essbee Vanhoutte
###WORK IN PROGRESS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
###Cuda_QS_Variant
###This code is pro-LGBTQ.

##References: I have borrowed many of the optimizations from here: https://stackoverflow.com/questions/79330304/optimizing-sieving-code-in-the-self-initializing-quadratic-sieve-for-pypy


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
from sympy import kronecker_symbol


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
k_max = 12#_000
b_max=10_000
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

def launch(n,primeslist,primeslist2):
    sbase=copy.deepcopy(primeslist[0:15])
  #  print("sbase: "+str(sbase)+" primeslist: "+str(primeslist))
    print("[i]Building psieve Residue Map (to do: some duplication here from merging two algos, fix later)")
    resmaps,resmaps2=build_residues(sbase,n)
    print("[i]Finished building... ")
    psieve(n,resmaps,resmaps2,sbase)
    
    return







def pow_mod(base, exponent, modulus):
    return pow(base,exponent,modulus)  










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




def solve_lin_con(a,b,m):
    ##ax=b mod m
    g=gcd(a,m)
    a,b,m = a//g,b//g,m//g
    return pow(a,-1,m)*b%m  
    
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

def calculate_modulus_root(n,prime_ind,mod_exp,sbase,k_test,resmaps,resmaps2,co,co2):
    prime=sbase[prime_ind]
    colist=resmaps[prime_ind][k_test%prime][1]
    colist2=resmaps2[prime_ind][k_test%prime][1]
    modulus=prime**mod_exp

    exponent=mod_exp
    co2_res=[]
    modr=[]
    col=[]

   # co2=[]
   # r2=[]
    j=0
    while j < len(colist2):
        if colist2[j]==co%(prime):
            if colist[j] not in co2_res:
                co2_res.append(colist[j])
        j+=1
    if len(co2_res) == 0:
       # print("wtf should never happen")
        #sys.exit()
        return -1,col
       # fail_glob=1
       # break
                
    i=0
    while i < len(co2_res):      
        c=co2_res[i]   
        r=-1

        roots=find_roots_poly([1,co,-n*k_test],prime)

        for r1 in roots:
            r1=lift_root2([1,co,-n*k_test],r1,prime,exponent)

            if evaluate([1,-c,n*k_test],r1)%prime==0:
                if r1 != 0:
      
                    col.append(c)
                    modr.append(r1)

                  #  if co2 == 1772:
                  #      print("c: "+str(c)+" roots: "+str(roots)+" co2_res: "+str(co2_res)+" r: "+str(r))
        i+=1
    return modr,col

def psieve(n,resmaps,resmaps2,sbase):
    dupe_map={}
    fb_map = {val: i for i, val in enumerate(sbase)}
    gcdtested_global=0
    gcdtested_global2=0
    debug_seen=[]

    mod_prime=3
    #*3*3*3*3*3#*3#*3#*3#*3*3*3*3
    bound=100
    min_exp=3
    #print("modulus: "+str(modulus)+" note: increase modulus if it fails")
    k_test=1
    while k_test <10_000:
       # print("trying k_test: "+str(k_test))#+" / "+str(modulus))
        if k_test%mod_prime==0:
            k_test+=1
            continue
        roots=find_roots_poly([1,0,n*k_test],mod_prime)
     #   print("k_test roots: "+str(roots))
        for r1b in roots:
            mod_exp=2
            while mod_exp < 1000:
                r1=lift_root2([1,0,n*k_test],r1b,mod_prime,mod_exp)
              #  print("**r1: "+str(r1)+" mod_exp: "+str(mod_exp))
                if r1 > bound and (mod_exp-1) < min_exp:
                    break
                elif r1 > bound and (mod_exp-1) > min_exp-1:
                    mod_exp-=1
                    r1=lift_root2([1,0,n*k_test],r1b,mod_prime,mod_exp)
                    modulus=mod_prime**mod_exp
                    print("Trying discriminant root: "+str(r1)+" finite field exponent: "+str(mod_exp)+" modulus p^e: "+str(modulus)+" k: "+str(k_test))
                
                    disc=r1**2+n*k_test
                    
                    if disc%(modulus)!=0:
                        print("fatal")

                    co_temp=(-r1*2)%(modulus)
                    i=0
                    while i < modulus and i < 500: 
                        co=(co_temp+modulus*i)
                        if co%2!=0:
                            i+=1
                            continue
                        step_dist=((co//2)+r1)//modulus
                        poly1=[1,co,-n*k_test]
                        pval=evaluate(poly1,r1)             
                        j=solve_lin_con(n*modulus,pval,modulus**2)
                        tpval2=pval-n*j*modulus
                        if tpval2%modulus**2!=0:
                            print("linear congruence failed1")
                        tpval2//=modulus**2
                        h=solve_lin_con(n,-step_dist**2+tpval2,modulus)
                        if abs((tpval2-n*h)-step_dist**2)%(mod_prime**(mod_exp)) != 0:
                            print("linear congruence failed2")
                            sys.exit()
                        u2=0
                        while u2 < 50: #To do: this can be turned into a 2d interval
                            interval=array.array("i",[1]*100_000) #to do: bit sieving
                            start_k=(k_test+j*modulus)+(modulus**2)*u2 
                            ind=0
                            while ind < len(sbase):
                                fprime=sbase[ind]

                                if fprime != mod_prime:
                                #marked_mod*=fprime
                                    colist=resmaps2[ind][start_k%fprime][1]
                                    r=0
                                    while r < fprime:
                           #         if r%fprime == 344%fprime and i == 2:
                            #            print("fprime: "+str(fprime)+" colist: "+str(colist)+" start_k: "+str(start_k))
                                        if r not in colist:
                                            dist=solve_lin_con(2*(modulus**2),r-co,fprime)
                                            while dist < len(interval):
                                                interval[dist]=0
                                                dist+=fprime
                                        r+=1
                                ind+=1
                       # print("interval: "+str(interval))
                            i2=0
                            while i2 < len(interval):
                                if interval[i2]==0:
                                    i2+=1
                                    continue
                                #print("Checking i2: "+str(i2))
                                co2=co+((2*modulus**2))*i2





                                new_k=start_k
                                fail=0
                                if g_debug == 1:
                                    ind=0
                                    while ind < len(sbase):
                                        fprime=sbase[ind]

                                        if fprime != mod_prime:
                                      #  if fprime == 5:
                                       #     debugroots=find_roots_poly([1,co,-n*new_k],fprime)
                                       #     if len(debugroots)>0:
                                    #    debug_roots.append(fprime)
                                    #    debug_roots.append([])
                                    #    for tr in debugroots:
                                    #        dis=solve_lin_con(modulus,tr-r1,fprime)
                                     #       debug_roots[-1].append(dis)
                                            #    print("co: "+str(co%fprime)+" debugroots: "+str(debugroots)+" u: "+str(u)+" k_test: "+str(k_test)+" co: "+str(co)+" co%modulus**2: "+str(co%(modulus**2))+" new_k: "+str(new_k)+" r1: "+str(r1%fprime)+" start_k: "+str((k_test+j*modulus))+" i: "+str(i)+" i2: "+str(i2))


                                            colist=resmaps2[ind][new_k%fprime][1]
                                            if co2%fprime not in colist:
                                                print("fatal error should never happen")
                                                sys.exit()
                                        ind+=1
                            
                                fprime=sbase[-1]
                                rr=find_roots_poly([1,co2,-n*new_k],fprime)
                                fexp=15
                                for rr2 in rr:
                                    rr2=lift_root2([1,co2,-n*new_k],rr2,fprime,fexp)
                                    if rr2%modulus == r1: ##To do: This condition is important.. should check how to force it though...
                                        gcdtest=math.gcd(rr2,n)
                                        gcdtested_global+=1
                                        fpval=evaluate([1,co2,-n*new_k],rr2)
                                    #print("gcdtested_global: "+str(gcdtested_global)+" fpval: "+str(fpval)+" rr2: "+str(rr2%modulus)+" r1: "+str(r1)+" (tpval2-n*h)-step_dist**2: "+str((tpval2-n*h)-step_dist**2))
                                        disc=co2**2+4*n*new_k
                                        #fpval=evaluate([1,co,-n*new_k],rr2)
                                        if gcdtest !=1 and gcdtest !=n:
                                            print("Factors of "+str(n)+" are: "+str(gcdtest)+" and "+str(n//gcdtest)+" new_k: "+str(new_k)+" disc: "+str(disc)+" u2: "+str(u2)+" b: "+str(co)+" gcdtested_global: "+str(gcdtested_global)+" fpval: "+str(fpval)+" rr2%modulus: "+str(rr2)+" r1: "+str(r1)+" i: "+str(i)+" i2: "+str(i2)+" interval: "+str(interval[i2]))
                                            sys.exit()

                               
                                i2+=1
                            u2+=1
                        i+=1
                    break
                mod_exp+=1      
        k_test+=1
    return

def find_r(mod,total):
    mo,i=mod,0
    while (total%mod)==0:
        mod=mod*mo
        i+=1
    return i



def enumerated_product(*args):
    yield from itertools.product(*(range(len(x)) for x in args))





def psieve_padic_solutions(prime,exp,n,k):
    blist=[]
    blist2=[]
    blist.append([prime,exp])
    blist2.append([prime,exp])
    rootlist=[] ##To do: precalculate..................
    rootlist2=[]
    p=0
    while p < prime**exp:
        disc=p**2-4*n*k
        nroots=find_roots_poly([1,0,-disc],prime)
        if len(nroots)<2:
            if len(nroots)==1:
                skip=0
                i=0
                while i < len(rootlist):
                    if rootlist[i]==p and rootlist2[i]==nroots[0]:
                        skip=1
                    i+=1
                if skip == 0:
              # if p not in rootlist:
                #if (len(rootlist)==0) or (rootlist[-1]!=p and rootlist2[-1] != nroots[0]): ##To do: Can change up the speed here by skipping these duplicated without needing this check..
                    rootlist.append(p)
               # if nroots[0] not in rootlist2:
                    rootlist2.append(nroots[0])
            p+=1
            continue
        for b in nroots:
            b=lift_root2([1,0,-disc],b,prime,exp)
            if b**2%prime**exp != disc%prime**exp:
                print("fatal: "+str(roots))
                sys.exit()
            roots=find_roots_poly([1,b,-n*k],prime)
            if len(roots)>1:
                for r in roots:
                 
                    r=lift_root2([1,b,-n*k],r,prime,exp)
                    if evaluate([1,b,-n*k],r)%prime**exp == 0:
                        skip=0
                        i=0
                        while i < len(rootlist):
                            if rootlist[i]==p and rootlist2[i]==b:
                                skip=1
                            i+=1
                        if skip == 0:
                        #if p not in rootlist:
                        #if (len(rootlist)==0) or (rootlist[-1]!=p and rootlist2[-1] != b): ##To do: Can change up the speed here by skipping these duplicated without needing this check..
                            rootlist.append(p)
                        #if b not in rootlist2:
                            rootlist2.append(b)
            elif len(roots)==1:
                der=get_derivative([1,b,-n*k])
                nb=evaluate(der,roots[0])
                r=find_roots_poly([1,nb,n*k],prime)
                for r2 in r:
                    r2=lift_root2([1,nb,n*k],r2,prime,exp)
                    if evaluate([1,nb,n*k],r2)%prime**exp == 0:
                        skip=0
                        i=0
                        while i < len(rootlist):
                            if rootlist[i]==p and rootlist2[i]==b:
                                skip=1
                            i+=1
                        if skip == 0:
                        #if p not in rootlist:
                       # if (len(rootlist)==0) or (rootlist[-1]!=p and rootlist2[-1] != b): ##To do: Can change up the speed here by skipping these duplicated without needing this check..
                            rootlist.append(p)
                        #if b not in rootlist2:
                            rootlist2.append(b)
        p+=1    

    
    blist.append(rootlist)
    blist2.append(rootlist2)

    return blist,blist2

def build_residues(sbase,n):
    resmaps=[]
    resmaps2=[]
    for prime in sbase:
        resmaps.append([])
        resmaps2.append([])
        k=0
        while k < prime and k < 1000:
            ##Hensel here isnt doing shit. I need to have a better look here first..........
            exp=1
            temp_plist,temp_plist2=psieve_padic_solutions(prime,exp,n,k)
           # if prime == 5:
           #     print("resmaps: "+str(temp_plist)+" k: "+str(k))
            if len(temp_plist)>0:
                resmaps[-1].append(temp_plist)
                resmaps2[-1].append(temp_plist2)
            k+=1
        #if prime == 5:
           # print("resmaps: "+str(resmaps))
   # sys.exit()
    return resmaps,resmaps2

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