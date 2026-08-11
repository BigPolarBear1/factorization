import pstats, cProfile
import sympy
import random
import itertools
import sys
import argparse
import time
import copy
from timeit import default_timer
import math
import gc
import array
import numpy as np
import os

min_lin_sieve_size=10_000
max_bound=10_000_000
key=0                 #Define a custom modulus to factor
build_workers=8
keysize=150           #Generate a random modulus of specified bit length
workers=1 #max amount of parallel processes to use
quad_co_per_worker=1 #Amount of quadratic coefficients to check. Keep as small as possible.
base=1_000
qbase=10
lin_size=1
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
    int_type=abs(int_type)
    length=0
    while(int_type):
        int_type>>=1
        length+=1
    return length 

def get_primes(start,stop):
    return list(sympy.sieve.primerange(start,stop))

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
    global keysize,key,workers,debug,base,lin_size,quad_size
    parser = argparse.ArgumentParser(description='Factor stuff')
    parser.add_argument('-key',type=int,help='Provide a key instead of generating one') 
    parser.add_argument('-keysize',type=int,help='Generate a key of input size')    
    parser.add_argument('-workers',type=int,help='# of cpu cores to use')
    parser.add_argument('-debug',type=int,help='1 to enable more verbose output')
    parser.add_argument('-base',type=int,help='Size of the factor base')
    parser.add_argument('-lin_size',type=int,help='Size of the factor base')
    parser.add_argument('-quad_size',type=int,help='Size of the factor base')
    args = parser.parse_args()
    if args.keysize != None:    
        keysize = args.keysize
    if args.key != None:    
        key=args.key
    if args.workers != None:  
        workers=args.workers
    if args.debug != None:
        debug=args.debug  
    if args.base != None:
        base=args.base  
    if args.lin_size != None:
        lin_size=args.lin_size  
    if args.quad_size != None:
        quad_size=args.quad_size   
    return

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

def power2(poly, f, p, exp):
    if exp == 1: return poly

    tmp = power2(poly, f, p, exp>>1)
    tmp = poly_prod(tmp, tmp)

    if exp&1:
        tmp = poly_prod(tmp, poly)
        return div_poly_mod(tmp, f, p)
    
    else: return div_poly_mod(tmp, f, p)

def gcd(a,b): # Euclid's algorithm ##
    if b == 0:
        return a
    elif a >= b:
        return gcd(b,a % b)
    else:
        return gcd(b,a)

def modinv(n,p):
    p2=p
    n = n % p
    x =0
    u = 1
    while n:
        x, u = u, x - (p // n) * u
        p, n = n, p % n
    return x%p2

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

def binomial_coeffs_fast(y, n):
    coeffs=[1]                
    c=1
    yk=1
    for k in range(1, n):
        c=c*(n-k+1)//k 
        yk*=y                
        coeffs.append(c*yk)
    return coeffs

def polygen(n,d):
    bx=math.ceil(n**0.5)
    bx=74

  #  by=math.ceil((bx**2-n)**0.5)
    f_x=binomial_coeffs_fast(-bx, d)
 #   g_x=binomial_coeffs_fast(by, 2)
    print("bx: "+str(bx)+" f_x: "+str(f_x))
    return f_x,bx#,g_x

def g_x_residues(fbase,n,d):
    exp=2
    g_x_res=[]
    i=0
    while i < len(fbase):
        j=0
        g_x_res.append({})
        while j <1:
            k=0
            while k < fbase[i]**exp:
                co=binomial_coeffs_fast(j, d)
                g_x=co+[-n*k]
                roots=find_roots_poly(g_x,fbase[i])
                roots.sort()
                if len(roots)>1:
                    new_roots=[]
                    for r in roots:
                        new_r=lift_root2(g_x, r, fbase[i], exp)
                        test=evaluate(g_x,new_r)
                        if test%fbase[i]**exp !=0:
                            print("fatal error: "+str(roots))
                            sys.exit()
                        new_roots.append(new_r)
                    new_roots.sort()
                    #print("prime: "+str(fbase[i])+" co: "+str(g_x)+" k: "+str(k)+" roots: "+str(roots))
                    try:
                        res=g_x_res[-1][k]
                        res.append([new_roots,j])
                    except Exception as e:
                        g_x_res[-1][k]=[[new_roots,j]]
                k+=1
            j+=1
        print("g_x: "+str(fbase[i])+" "+str(g_x_res[-1]))
        i+=1
    return g_x_res

def evaluate(f, x):
    res = 0

    for i in range(len(f)-1):
        res += f[i]
        res *= x

    res += f[-1]

    return res

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

def get_quadratic_residues(bin,n,fbase):
  #  print("fbase: "+str(fbase))
    d=2
    f_x_res=[]
    i=0
    while i < len(fbase):
        prime=fbase[i]
        f_x_res.append(prime**2)
        f_x_res.append([])
        co=binomial_coeffs_fast(-(bin), d)
        co+=[0]
       # der=get_derivative(co)
        k=0
        while k < prime:
            f_x=co
            f_x[-1]=n*k
           # print("f_x: "+str(f_x))
            f_xt=copy.deepcopy(f_x)
            roots=find_roots_poly(f_xt,prime)#to do: remove deepcopy.. just make sure it never gets altered
            if len(roots)==1:
                
                if roots[0]!=bin%prime:
                    print("something fucked up")

                ###Left to degree 2
                k2=k
                while k2 <  prime**2:
                    r=bin%prime**2
                    f_x[-1]=n*k2
                    if evaluate(f_x,r)%prime**2 ==0:
                        f_x_res[-1].append(k2)
                    k2+=prime

            
            k+=1

    #    print("f_x: "+str(f_x)+" prime: "+str(fbase[i])+" "+str(f_x_res[-1]))
        i+=1

    return f_x_res

def get_quartic_residues(bin,n,fbase):
    d=4
    f_x_res=[]
    i=0
    while i < len(fbase):
        prime=fbase[i]
        f_x_res.append(prime)
        f_x_res.append([])
        co=binomial_coeffs_fast(-(bin), d)
        co+=[0]
       # der=get_derivative(co)
        k=0
        while k < fbase[i]:
            f_x=co
            f_x[-1]=n*k
           # print("f_x: "+str(f_x))
            f_xt=copy.deepcopy(f_x)
            roots=find_roots_poly(f_xt,prime)#to do: remove deepcopy.. just make sure it never gets altered
            if len(roots)==1:
               # f_x_res[-1].append([k,roots,1,0])
                if roots[0]!=bin%prime:
                    print("something fucked up")
                ###Left to degree 4
                k2=k
                while k2 <  prime**4:
                    r=bin%prime**4
                    f_x[-1]=n*k2
                    if evaluate(f_x,r)%prime**4 ==0:
                        f_x_res[-1].append(k2)
                    k2+=prime
            k+=1

     #   print("f_x: "+str(f_x)+" prime: "+str(fbase[i])+" "+str(f_x_res[-1]))
        i+=1
    return f_x_res

def compute_result(quad_res,quar_res,bin,fbase,n):
    i=0
    while i < len(fbase):
        prime=fbase[i]

        lquad=quad_res[i]
        lquar=quar_res[i]
        if len(lquad)!=len(lquar):
            print("fatal")
            sys.exit()
        j=0
        while j < len(lquad):
            k1=lquad[j][0]
          #  k1=19*23
            k2=lquar[j][0]
            k2=(k2-(bin**2*k1))%prime**4
            mod_inv_k1=modinv(k1,prime**4)
            print(k2)
            k2=(k2*mod_inv_k1)%prime**4
            print("mod_inv_k1: "+str(mod_inv_k1))
            print("prime: "+str(prime)+" quad_res: "+str(lquad[j])+" k: "+str(lquad[j][0]))
            print("prime: "+str(prime)+" quar_res: "+str(lquar[j])+" k: "+str(lquar[j][0]))
            print("solution: "+str(k2)+" k1: "+str(k1)+" k2: "+str(lquar[j][0])+" bin: "+str(bin))
            k2_sqr=math.isqrt(k2)
            if k2_sqr**2==k2 and k2%n == bin**2%n:
                print("k2_sqr: "+str(k2_sqr)+" bin: "+str(bin))
                gcdtest=math.gcd(k2_sqr+bin,n)
                if gcdtest != 1 and gcdtest != n:
                    print("factors of N are: "+str(gcdtest)+" and "+str(n//gcdtest))
                    sys.exit()
            j+=1
        i+=1
    return


def generate_modulus(n,primeslist,seen,tnum,close_range,too_close,LOWER_BOUND_SIQS,UPPER_BOUND_SIQS,tnum_bit):
    const_1=1_000
    const_2=1_000_000

    small_B = base#len(primeslist)
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
        while counter2 < const_1:
            counter2+=1
            found_a_factor = False
            counter=0
            while(found_a_factor == False) and counter < const_2:
                randindex = random.randint(lower_polypool_index, upper_polypool_index)
               # if  jacobi((-quad*n)%primeslist[randindex],primeslist[randindex])!=1:
                  #  counter+=1
                  #  continue
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
            cfact.append(math.isqrt(potential_a_factor))
        #    if  jacobi((-quad*n)%primeslist[randindex],primeslist[randindex])!=1:#hmap[randindex][1]!=quad%primeslist[randindex]:
       #         print("THE FUC")
      #          time.sleep(1000000)
            indexes.append(randindex)
            j = tnum_bit - cmod.bit_length()
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
        while not found_a_factor and counter3< const_1 and randindex <base:
     #       if  jacobi((-quad*n)%primeslist[randindex],primeslist[randindex])!=1:
     #           randindex += 1
    #            counter3+=1
   #             continue
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
      #  if  jacobi((-quad*n)%primeslist[randindex],primeslist[randindex])!=1:
      #      print("THE FUC: ",randindex)
       #     time.sleep(1000000)
        cfact.append(math.isqrt(potential_a_factor))
        indexes.append(randindex)

        diff_bits = (tnum - cmod).bit_length()
        if diff_bits < tnum_bit:
            if cmod in seen:
                continue
            else:
                seen.append(cmod)
                return cmod,cfact,indexes
    return 0,0,0

def generate_sieve_interval():



    return

def factorise_fast(value,factor_base):
    if value == 0:
        print("blah")
        return [],-1
    factors = set()
    if value < 0:
        factors ^= {-1}
        value = -value
    while value % 2 == 0:
        factors ^= {2}
        value //= 2

    length=factor_base[0]#len(factor_base)#factor_base[0]
    i=1
    while i < length:
        factor=factor_base[i]
        while value % factor == 0:
            factors ^= {factor}
            value //= factor
        i+=1
    return factors, value


def solve_lin_con(a,b,m):
    ##ax=b mod m
    #g=gcd(a,m)
    #a,b,m = a//g,b//g,m//g
    return pow(a,-1,m)*b%m  

def build_interval(fbase,init_k,odd_mod,bin,n,smoothcan_org):
    f_x=binomial_coeffs_fast(bin, 2)
    f_x+=[0]
 #   print("f_x: "+str(f_x))
    interval=np.zeros(100_000,dtype=np.uint16)   
    i=1
    length=fbase[0]
    while i < length:
        k=0
        prime=fbase[i]
     #   print(prime)
        if odd_mod%prime==0:
            i+=1
            continue
        odd_mod_inv=modinv(odd_mod,prime)
        while k < prime:
            #f_x[-1]=k*n
            #f_x[-1]=(f_x[-1]*odd_mod_inv)%prime
            #f_x[0]=odd_mod
            #f_xc=copy.deepcopy(f_x)
            test=bin**2-n*k
            test*=odd_mod_inv
            test%=prime
            if compute_legendre_character((test)%prime,prime)==-1:
         #   roots=find_roots_poly(f_xc,prime)
          #  if len(roots)==0:
                count=solve_lin_con(odd_mod,k-init_k,prime)
                if (init_k+odd_mod*count)%prime != k:
                    print("euhm??")
                    sys.exit()
                test=bin**2-n*(init_k+odd_mod*count)
                test*=odd_mod_inv
                test%=prime
                if compute_legendre_character((test)%prime,prime)!=-1:
                    print("failed: "+str(f_x)+" prime: "+str(prime)+" k: "+str(k)+" bin: "+str(bin))
                    sys.exit()
             #   else:
             #       print("ok")
                interval[count::prime]=1
            k+=1
        i+=1


    return interval

def check_higher_degrees(smoothcan_org,n,fbase,ret_array,bin,odd_mod,local_factors_org,mod):
    ##Very rudimentary.. needs to be improved now
    degree=4
    while degree < 20:
        ccount=0
        init_k=bin**degree-smoothcan_org**(degree//2)
        if init_k%n!=0:
            print("fatal")
            sys.exit()
        init_k=init_k//n
        newcan=bin**degree-n*(init_k+odd_mod)
        count_start=newcan//(n*odd_mod)
    #   newcan2=bin**4-n*(init_k+odd_mod*count_start)
    #   print(str(newcan)+" "+str(newcan2))
        interval=build_interval(fbase,init_k,odd_mod,bin**(degree//2),n,smoothcan_org)
        w=0
        while w < len(interval):
            if interval[w]==0:
                test=bin**degree-n*(init_k+odd_mod*w)
                if test%odd_mod !=0:
                    print("fatal")
                    sys.exit()
             #   print("test: "+str(test)+" odd_mod: "+str(odd_mod))
                test=test//odd_mod
                length=fbase[0]
                e=1
                while e < length:
                    l=compute_legendre_character(test,fbase[e])
                    if l == -1:
                        print("l: "+str(l)+" prime: "+str(fbase[e]))
                        sys.exit()
                    e+=1
                test2=math.isqrt(test)
                if test2**2 == test:

                #    print("found "+str(math.gcd(test2+bin**(degree//2),n))+" test2: "+str(test2**2%n)+" "+str((bin**(degree))%n))
                    newcan=bin**degree-n*(init_k+odd_mod*w)
                    newcan=newcan*smoothcan_org
                    root=bin**(degree//2)*bin
                    testf=math.isqrt(abs(newcan))
                    if testf**2 == newcan and root%n != testf%n and (root+testf)%n !=0:
                        f1=math.gcd(testf+root,n)
                        f2=n//f1
                        print("succeeded early: "+str(f1)+" and "+str(f2))
                        sys.exit()
                else:
                    print("not found")
            w+=1
       # while ccount < 5:
        #    count=count_start+ccount
       #     ccount+=1
       #    # root=bin**(degree//2)
      #      newcan=bin**degree-n*(init_k+odd_mod*count)
    #    if bitlen(newcan//odd_mod) < 32:
    #        print(str(bitlen(newcan//odd_mod))+" odd_Mod: "+str(bitlen(odd_mod)))
      #      if newcan%odd_mod !=0:
      #          print("fatal")
      #          sys.exit()
       
     #   local_factors, value = factorise_fast(newcan,fbase)
      #      root=bin**(degree//2)*bin
      #      newcan=newcan*smoothcan_org
      #      if root**2%n!=newcan%n:
      #          print("super fatal")
      #          sys.exit()
      #      local_factors2, value2 = factorise_fast(newcan,fbase)
      #      local_factors3, value3 = factorise_fast(newcan//(smoothcan_org*odd_mod),fbase)
      #      test=math.isqrt(value2)
      #  print("newcan: "+str(newcan))
        #    if test**2 == value2 and local_factors2 not in ret_array[2]:
        #        print("**Smooths: "+str(len(ret_array[0]))+" local_factors: "+str(local_factors2)+" degree: "+str(degree)+" "+str(newcan/(smoothcan_org*odd_mod))+" odd_mod: "+str(odd_mod)+" "+str(value3))#+" local2: "+str(local_factors2)+" value2: "+str(value2)+" value1: "+str(value)+" local_org: "+str(local_factors_org))
        #        ret_array[1].append(root**2)
        #        ret_array[0].append(newcan)
        #        ret_array[2].append(local_factors2)
        #        ret_array[3].append([])
        degree+=2
        
    return

def process_sieve_interval(k,n,bin,mod,fbase,ret_array):
    ##Barebones, improve later
    if (bin**2-k*n)%mod !=0:
        print("fatal error123: "+str(mod)+" "+str(bin)+" "+str(math.gcd(bin,mod))+" k: "+str(k))
        return
        sys.exit()

    i=0
    while i < lin_size:
        smoothcan=(bin+mod*i)**2-k*n
        
        local_factors, value = factorise_fast(smoothcan,fbase)
        odd_mod=1
        for odd in local_factors:
            if odd != -1 and odd != 2:
                odd_mod*=odd
        odd_mod*=value
       # print("smoothcan "+str(smoothcan)+" "+str(local_factors))

        
        test=math.isqrt(value)
        if value**2 == test and local_factors not in ret_array[2]:
            print("local_factors: "+str(local_factors)+" bin: "+str(bin)+" smoothcan: "+str(smoothcan)+" mod: "+str(mod))
            check_higher_degrees(smoothcan,n,fbase,ret_array,bin+mod*i,odd_mod,local_factors,mod)    
            print("Smooths: "+str(len(ret_array[0]))+" local_factors: "+str(local_factors)+" value: "+str(value))
           
            ret_array[1].append((bin+mod*i)**2)
            ret_array[0].append(smoothcan)
            ret_array[2].append(local_factors)
            ret_array[3].append([])
        i+=1

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
       # print(" polyval sqrt: "+str(sqrt_left%N)+" disc constant sqrt: "+str(sqrt_right%N))#+" zx*zxy: "+str(x))
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

def solve_roots(prime,n):
    hmap_p=[]
    hmap_p2=[]
    sq=[1,0,-37]
    sqr=find_roots_poly(sq, prime)
    
    if len(sqr)==0:
        return [],[]
    print("SQR: "+str(sqr))
    div_37=modinv(37,prime)
    k=0
    while k < prime:# and k < quad_sieve_size+1:
        hmap_p2.append({})
        hmap_p.append({})
        coeff=[(-n)%prime]
     #   coeff2=[(-n)%prime]


        d=2
        d_ind=0
        while d_ind < d:
            coeff.insert(0,0)

            d_ind+=1

      #  ranges = [range(start, prime) for start in coeff[:-1]]
        predef=[]
        i=0
        while i < prime:

            combo=[k,(sqr[0]*i)%prime]
         #   combo2=[k,-i]
        #ranges = [(1,1),(0,prime)]
        #for combo in itertools.product(*ranges):
            current = list(combo) +  [coeff[-1]]
         #   current2= list(combo2) +  [coeff2[-1]]
        #    current[0]=k
          #  current2[0]=k
          #  if current[0]==0:
           #     if current[1]==0:
          #          continue
           #     x=solve_lin_con(current[1],n*k,prime)
          #      if (x*current[1]-n*k)%prime != 0:
         ##           print("fatal")
         #       root=[x]
       #     else:
              #  ##to do: Use tonelli instead if we arn't going to use higher degrees...

            if (27)%prime == i and k == 37%prime:
                print("eval: "+str(current))
            root = find_roots_poly(current, prime)
         #   root2 = find_roots_poly(current2, prime)
            if len(root)>0:
                    f=[-i]
                    f.append(root)
                    try:
                        res=hmap_p2[-1][i]
                        print("shouldn't happen")
                        sys.exit()
                    except Exception as e:
                        hmap_p2[-1][i]=[root]
               #     hmap_p2[-1].append(f)
          #  if len(root2)>0:
            #        f=[i]
           #         f.append(root2)
            #      #  hmap_p[-1].append(f)
            #        try:
            #            res=hmap_p[-1][i]
                #        print("shouldn't happen")
               #         sys.exit()
               #     except Exception as e:
               #         hmap_p[-1][i]=[root2]
            i+=1
        if len(hmap_p2[-1]) > 0 and k == 1%prime:
            print("- prime: "+str(prime)+" k: "+str(k)+" "+str(hmap_p2[-1]))
            try:
                res=hmap_p2[-1][27%prime]
                print("- prime: "+str(prime)+" k: "+str(k)+" "+str(res))
            except Exception as e:
                print("FAILED")
                sys.exit()
                disc=37*(27)**2+4*n
                leg=compute_legendre_character(disc%prime,prime)
                if leg != -1:
                    print("asdjasd")
                    sys.exit()
                print("failed: ",leg)
       # if len(hmap_p[-1]) > 0 and k == 23%prime:
        #    print("+ prime: "+str(prime)+" k: "+str(k)+" "+str(hmap_p[-1]))
        k+=1
  #  print("\n")
    return hmap_p,hmap_p2

def create_hashmap(n,primeslist):
    i=0
    hmap=[]
    hmap2=[]

    while i < len(primeslist):
        hmap_p,hmap_p2=solve_roots(primeslist[i],n)
        hmap.append(hmap_p)
        hmap2.append(hmap_p2)
       # print("prime: "+str(primeslist[i])+" hmap_p: "+str(hmap_p))

        i+=1

    return hmap,hmap2

def psieve(n,fbase):
    ret_array=[[],[],[],[]]
    fbase_opt=copy.copy(fbase)
    fbase_opt.insert(0,len(fbase_opt)+1)
    fbase_opt=array.array('q',fbase_opt)
    fbase_fin=copy.copy(fbase)
    fbase_fin.insert(0,2) ##To do: remove when we fix lifting for powers of 2
    fbase_fin.insert(0,-1)
 #   print("[i]Building residue map")
 #   hmap,hmap2=create_hashmap(n,fbase) ##Fix this later
   # sys.exit()
  #  print("[i]Building interval")
    seen=[]
    close_range=20
    too_close=1
    LOWER_BOUND_SIQS=1
    UPPER_BOUND_SIQS=4000
    tnum=int(((n)**0.2) /(1))
    found=0
    while 1:


        mod,cfact,indexes=generate_modulus(n,fbase,seen,tnum,close_range,too_close,LOWER_BOUND_SIQS,UPPER_BOUND_SIQS,bitlen(tnum))
        #mod=13
        if mod == 0:
            print("failed to generate modulus")
            sys.exit()
        div=mod
    #found=0
   # div=1
    #while div <100000:
       # local_factors, value = factorise_fast(div,fbase_opt)
     #   if value !=1:
          #  div+=1
         #   continue
        print("[i]Building interval with divisor: "+str(div))
        interval=np.ones([lin_size,lin_size],dtype=np.int16)

      #  modi=modinv(div,n)
        t=0
        while t < len(fbase):
            prime=fbase[t]
            sq=[1,0,-div]
            sqr=find_roots_poly(sq, prime)
    
            if len(sqr)==0:
                t+=1
                continue

            k=0
            while k < prime:
                i=0 
                while i < prime:
                    if math.gcd(div,prime)!=1:
                        i+=1
                        continue
                    combo=binomial_coeffs_fast(i,2)
                  #  combo[0]*=div
                    combo[1]*=sqr[0]
                    combo[1]%=prime
                    current = list(combo)+[-n*k]

                 ##   print("current: "+str(current)+" prime: "+str(prime))
                    curc=copy.deepcopy(current)
                    root = find_roots_poly(curc, prime)

                    if len(root)==0:

                        if (k*current[0])%prime < lin_size and i < lin_size:
                            interval[(k*current[0])%prime::prime,i::prime]=0
                   # else:
                    #    if len(sqr_div)>0:
                    #        divinv=modinv(div,prime)
                    #        disc=(sqr_div[0]*current[1])**2+n*k*current[0]
                   #         leg=compute_legendre_character((disc*modular_div)%prime,prime)
                    #        print(" disc: "+str(disc)+" divinv: "+str(divinv)+" Prime: "+str(prime)+" current: "+str(current)+" root: "+str(root)+" sqdiv: "+str(sqr_div)+" legendre: "+str(leg))

                   
                    i+=1
                k+=1
            t+=1
       # sys.exit()
        print("[i]Checking interval")

        indexlist=np.nonzero(interval)
        indexlist_x=indexlist[1]
        indexlist_y=indexlist[0]
        ind=0
        length=len(indexlist_x)

        while ind < length:# length:  
            y=int(indexlist_x[ind])
            k=int(indexlist_y[ind])

  
            if k == 0:
                ind+=1
                continue
        #    f_x=binomial_coeffs_fast(y, 2)
      #  f_x+=[n*k]
            disc=div*(2*y)**2+4*n*k
            disc_test=math.isqrt(disc)

            bsmooth=(div*2*y)**2+4*n*k*div
          #  print("disc: "+str(disc)+" (bsmooth//div)**0.5: "+str((bsmooth//div)**0.5)+" bsmooth: "+str(bsmooth))
            if disc_test**2 !=disc:
                ind+=1
                continue

            local_factors, value = factorise_fast(div,fbase_opt)
            print("value: "+str(value))

            test=math.isqrt(value)

            if test**2 == value:# and y**2 not in ret_array[1]:
                print("**Smooths: "+str(len(ret_array[0])))
        #        print("**Smooths: "+str(len(ret_array[0]))+" local_factors: "+str(local_factors2)+" degree: "+str(degree)+" "+str(newcan/(smoothcan_org*odd_mod))+" odd_mod: "+str(odd_mod)+" "+str(value3))#+" local2: "+str(local_factors2)+" value2: "+str(value2)+" value1: "+str(value)+" local_org: "+str(local_factors_org))
                ret_array[1].append((div*2*y)**2)
                ret_array[0].append(bsmooth)
                ret_array[2].append(local_factors)
                ret_array[3].append([])
                found+=1
            if found > 0:
                test,test2=QS(n,fbase_fin,ret_array[0],ret_array[2],ret_array[1],ret_array[3]) 
                if test !=0:
                    print("\n\n\n\nFound at: ",len(ret_array[0]))
                    return 
                found=0
          #  sys.exit()
           # u=math.isqrt(disc)
          #  gcdtest=math.gcd(y+u,n)
          #  if gcdtest != 1 and gcdtest != n:
           #     print("factors are: "+str(gcdtest)+" and "+str(n//gcdtest))
           #     sys.exit()
            ind+=1
    return

def main():
    global key
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
    psieve(n,primeslist1)
    duration = default_timer() - start
    print("\nFactorization in total took: "+str(duration))

if __name__ == "__main__":
    parse_args()
    print_banner()
    main()