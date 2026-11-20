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
    sbase=copy.deepcopy(primeslist[0:7])
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

def calculate_modulus_root(n,prime_ind,mod_exp,sbase,k_test,resmaps,resmaps2,co):
    prime=sbase[prime_ind]
    colist=resmaps[prime_ind][k_test%prime][1]
    colist2=resmaps2[prime_ind][k_test%prime][1]
    modulus=prime**mod_exp

    exponent=mod_exp
    co2_res=[]
    modr=[]
    j=0
    while j < len(colist2):
        if colist2[j]==co%(prime):
            if colist[j] not in co2_res:
                co2_res.append(colist[j])
        j+=1
    if len(co2_res) == 0:
       # print("wtf should never happen")
        #sys.exit()
        return -1
       # fail_glob=1
       # break
                
                
    for c in co2_res:


        r=-1

        roots=find_roots_poly([1,co,-n*k_test],prime)

        for r1 in roots:
            r1=lift_root2([1,co,-n*k_test],r1,prime,exponent)

            if evaluate([1,-c,n*k_test],r1)%prime==0:
                if r1 != 0:
                    r=r1
                           # debug.append(r1)
                          #  if len(debug)>1:
                           #     print("asddasd: "+str(roots)+" prime: "+str(prime))

        if r == -1:
            continue
     #   s=solve_lin_con(n,r,modulus)
        modr.append(r)

    return modr

def psieve(n,resmaps,resmaps2,sbase):
    dupe_map={}
    fb_map = {val: i for i, val in enumerate(sbase)}
    gcdtested_global=0
    debug_seen=[]
    mod_exp=6
    mod_prime=5
    modulus=mod_prime**mod_exp#*3*3*3*3*3#*3#*3#*3#*3*3*3*3
 #   print(fb_map)
    mod_prime_ind=fb_map[mod_prime]
    if sbase[mod_prime_ind]!=mod_prime:
        print("fatal error")
        sys.exit
    print("modulus: "+str(modulus))
    k_test=1
    while k_test <100:

        marked_mod=1
        clist1=[]
        clist2=[]
        i=0
        while i < len(sbase):
            prime=sbase[i]
            if n%prime ==0 or modulus%prime == 0:
                i+=1
                continue
            clist1.append(resmaps[i][k_test%prime][0][0]**resmaps[i][k_test%prime][0][1])
            colist=resmaps[i][k_test%prime][1]
            clist1.append(colist)
            clist2.append(resmaps2[i][k_test%prime][0][0]**resmaps2[i][k_test%prime][0][1])
            colist2=resmaps2[i][k_test%prime][1] 
            clist2.append(colist2)
            marked_mod*=resmaps[i][k_test%prime][0][0]**resmaps[i][k_test%prime][0][1]
            i+=1
      #  print("clist1: "+str(clist1)+" clist2: "+str(clist2))
        clist1=get_partials(marked_mod,clist1)
        clist2=get_partials(marked_mod,clist2)
    
        cenum=[]
        cenum2=[]
        i=0
        while i < len(clist2):
            cenum.append(clist2[i+1])
            cenum2.append(clist1[i+1])
            i+=2
        for idx in enumerated_product(*cenum):
            co=0
            co2=0
            i=0
            while i < len(idx):
                ind=idx[i]
                co+=cenum[i][ind]
                co2+=cenum2[i][ind]
                i+=1
            co%=marked_mod
            co2%=marked_mod
     #       if co%marked_mod == (672*2)%marked_mod and co2%marked_mod == (886*2)%marked_mod: 
      #          print("!!!!!!!!!!!!!!!!!!!!!!!!![i]Trying k: "+str(k_test)+" co: "+str(co)+" co2: "+str(co2)+" marked_mod: "+str(marked_mod))#+" cenum: "+str(clist1)+" cenum2: "+str(clist2))    

            gcdtested=0           
          #  modrl=calculate_modulus_root(n,mod_prime_ind,mod_exp,sbase,k_test,resmaps,resmaps2,co)
          #  print("modrl: "+str(modrl))
          #  if modrl == -1:
          #      continue
          #  for modr in modrl:
             #   interval=array.array('i',[0]*1)
                #fail_glob=0
               # i=0
               # while i < len(sbase):
                  #  prime=sbase[i]
                  #  if n%prime ==0 or modulus%prime == 0:
                  #      i+=1
                  #      continue





                 #   exponent=1

                #    hit=0

                #    roots=find_roots_poly([1,co,-n*k_test],prime)
                #    for r1 in roots:
                        #r1=lift_root2([1,co,-n*k_test],r1,prime,exponent)

                  #      if evaluate([1,-co2,n*k_test],r1)%prime==0:
                    #        if r1 != 0:
                                #dist=solve_lin_con(modulus,r1-modr,prime)
                                #while dist < len(interval):
                                #    interval[dist]+=1
                                #    dist+=prime
                    #            if r1%(prime**exponent) == ((modr)%modulus)%(prime**exponent):
                    #                hit=1

                    #if hit==0:
                    #    fail_glob=1
                    #    break
                #    i+=1 
             #   if co%marked_mod == (672*2)%marked_mod and co2%marked_mod == (886*2)%marked_mod: 
             #       print("yo yo: "+str(modr)+" marked_mod: "+str(marked_mod))
                #if fail_glob == 1:
                #    continue

             #   print("checking modrl: "+str(modrl))#+" co_res: "+str(co2_res_final))
             #   if (co+co2)%n !=0 and math.gcd(modr,n)==1:# and modr%3 ==0:
             #       print("checking root: "+str(modr)+" co: "+str(co)+" co2: "+str(co2)+" co: "+str(co%n)+" co2: "+str(co2%n)+" marked_mod: "+str(marked_mod))#+" co_res: "+str(co2_res_final))
             #       if modr in debug_seen:
             ##           print("!!!!!!!!!!!hit")
             #       else:
             #           debug_seen.append(modr)


                ##To do: If a root is a multiple of another root.. that would work too.... let me think about this...

            w=0
            while w < 1:
                h = 0
                while h < 1:
                    fco=co+marked_mod*h 
                    fco2=co2+marked_mod*w  
                    disc=fco2**2-fco**2
                    ndiv=modinv(4*n,marked_mod)
                    fk=(disc*ndiv)%marked_mod
                    if fk%marked_mod != k_test%marked_mod:
                        print("hopefully this never hits and we can substitute")
                        sys.exit()
                    if math.gcd(fk,marked_mod)!=1:
                        h+=1
                        continue
                    fmodrl=calculate_modulus_root(n,mod_prime_ind,mod_exp,sbase,fk,resmaps,resmaps2,fco)
                    if fmodrl == -1:
                        h+=1
                        continue
                    poly1=[1,fco,-n*fk]
                    poly2=[1,-fco2,n*fk]
                    rres=[]
                    i=0
                    while i < len(sbase):
                        prime=sbase[i]
                        if n%prime ==0 or modulus%prime == 0:
                            i+=1
                            continue
                            
                        rres.append(prime)
                        rres.append([])
                        roots=find_roots_poly([1,fco,-n*fk],prime)
                        for tr in roots:
                            if evaluate(poly2,tr)%prime ==0:
                                if len(roots)>1 and tr == 0:
                                    continue
                                rres[-1].append(tr)

                        i+=1
                    for fr in fmodrl:
                        dist=[]
                        i=0
                        while i < len(rres):
                            dist.append(rres[i])
                            if len(rres[i+1])>1:
                                print("really shouldn't happen: "+str(rres[i])+" rres[i+1]: "+str(rres[i+1]))
                                #sys.exit()
                            dis=solve_lin_con(modulus,rres[i+1][0]-fr,rres[i])
                            dist.append([dis])
                            i+=2
                        dist=get_partials(marked_mod,dist)
                        fdist=0
                        i=0
                        while i < len(dist):
                            fdist+=dist[i+1][0]
                            i+=2
                        fdist%=marked_mod
                       # if fdist > marked_mod**0.5:
                       #     continue
                        r=fr+modulus*fdist
                        eval1=evaluate(poly1,r)
                        eval2=evaluate(poly2,r)

                        #TO DO: this should evaluate to 0 % modulus as well... which is not always the case... room for improvement

                        if eval1%(marked_mod*modulus)==0 and eval2%(marked_mod*modulus) ==0:
                            gcdtest=math.gcd(r,n)      
                            gcdtested+=1
                            gcdtested_global+=1  
                            if gcdtest != n and gcdtest != 1:# and:# r == 214:# not in debug_seen:
                                print("**Found with psieve, factors of "+str(n)+" are: "+str(gcdtest)+" and "+str(n//gcdtest)+" gcdtested: "+str(gcdtested)+" gcdtested_global: "+str(gcdtested_global)+" fco: "+str(fco)+" fco2: "+str(fco2)+" root: "+str(r)+" k: "+str(k_test)+" marked_mod: "+str(marked_mod)+" h: "+str(h)+" w: "+str(w)+" fdist: "+str(fdist))
                                sys.exit()
                    h+=1
                w+=1



                 #   print("valid x: "+str(modr)+" co: "+str(co)+" co2: "+str(co2))
               # root=(modr)%modulus
               # o=0
               # while o < 1:#len(interval):
                 #   if interval[o]!=len(r_res)//2:
                 #       o+=1
                #        continue
               #     gcdtest=math.gcd(root+o*modulus,n)
               #     gcdtested+=1
               #     gcdtested_global+=1
               #     if gcdtest != n and gcdtest != 1:# and (root+o*modulus) not in debug_seen:
                           # print("interval: "+str(interval[o])+" #: "+str(len(r_res)//2)+" o: "+str(o))
                 #       print("Found with psieve, factors of "+str(n)+" are: "+str(gcdtest)+" and "+str(n//gcdtest)+" gcdtested: "+str(gcdtested)+" gcdtested_global: "+str(gcdtested_global)+" co: "+str(co)+" co2: "+str(co2)+" root: "+str(root+o*modulus)+" o: "+str(o)+" k: "+str(k_test)+" marked_mod: "+str(marked_mod))
                     #   debug_seen.append(root+o*modulus)
                      #  sys.exit()
                #    o+=1
 

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