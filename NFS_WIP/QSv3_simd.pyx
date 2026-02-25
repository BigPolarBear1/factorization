##NFS related code is borrowed from: https://github.com/basilegithub/General-number-field-sieve-Python (note: Very impressively written, helped me big time, thanks)



import argparse
import cProfile
from numpy.polynomial.legendre import leggauss
import os
from datetime import datetime
import random
import itertools
import multiprocessing
import time
import copy
from timeit import default_timer
import numpy as np
from scipy.linalg import det
import math, sys
import sympy
import gc
from cpython cimport array
import array
cimport cython
from sympy import symbols, Poly, Matrix


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
## Set path to find config file

CONFIG_PATH = os.path.abspath("../config/config.ini")


##Key gen function ##
def power2(x, y, p):
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
    x = power2(a, d, n);
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
    if math.gcd(a, m) != 1:
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
        if math.gcd(e, (p - 1) * (q - 1)) == 1:
            break

    phi=(p - 1) * (q - 1)
    d = findModInverse(e, (p - 1) * (q - 1))
    publicKey = (n, e)
    privateKey = (n, d)
    print('[i]Public key - modulus: '+str(publicKey[0])+' public exponent: '+str(publicKey[1]))
    print('[i]Private key - modulus: '+str(privateKey[0])+' private exponent: '+str(privateKey[1]))
    return (publicKey, privateKey,phi,p,q)
##END KEY GEN##

# Naive smoothness test: trial division by every prime in the base  
# Note that in this function, we test for smoothness (candidate) and not abs(candidate), see QS.sieve_and_smooth
def trial(pair, primes, const1, const2, div):
    large1 = 1
    large2 = 1
    
    result = abs(pair[3])
    for p in primes:
        while not result%p: 
            result //= p
        if result == 1: 
            break

    if result > 1:
        if result < const1: large1, large2 = 1, result
        elif result < const2 and fermat_primality(result):
            result = pollard_rho(result)
            large1, large2 = min(result), max(result)
        else: 
            return False, 1, 1, 1, 1
    else: 
       # print("wohooo1")
        large1, large2 = 1, 1
   # print("hit")
    result = abs(pair[1])//div
   # print("result: ",result)
    for p in primes:
        while not result%p: 
            result //= p
        if result == 1: 
          #  print("WOHOOO!!")
            return True, 1, 1, large1, large2

    if result > 1:
        if large2 > 1: return False, 1, 1, 1, 1
        elif result < const1: return True, 1, result, 1, 1
        elif result < const2 and fermat_primality(result):
            result = pollard_rho(result)
            return True, min(result), max(result), 1, 1
        else: return False, 1, 1, 1, 1
        
    return True, 1, 1, 1, 1
    
    
# Pollard rho algorithm to factor small primes
# This is used to find the two large primes when needed
def pollard_rho(n):
    a = int(random.randint(1, n-3))
    s = int(random.randint(0, n-1))
    x = s
    y = s
    d = 1
    while d == 1:
        e = 1
        X = x
        Y = y
        for k in range(0, 100):
            x = (x*x+a)%n
            y = (y*y+a)%n
            y = (y*y+a)%n
            e = e*(x-y)%n
        d = math.gcd(e, n)
    if d == n:
        x = X
        y = Y
        d = 1
        while d == 1:
            x = (x*x+a)%n
            y = (y*y+a)%n
            y = (y*y+a)%n
            d = math.gcd(x-y, n)
    if d == n:
        return pollard_rho(n)
    return d, n//d
# Efficiently performs d[i] ^= 1 for all i
def switch_indices(d):
    return ~d
    
# d selects the coefficients of A
def multiply_d(A, d):
    return [A[i]&d for i in range(len(A))]
    
# Compute W_inv and the indices for d
# Basically performs gaussian elimination
def block(T, N):     
    M = concatenate(T, identity(N), N)
    S = []
    
    for j in range(N):
        for k in range(j, N):
            if (M[k] >> 2*N - j - 1)&1 != 0:
                M[k], M[j] = M[j], M[k]
                break
        
        if (M[j] >> 2*N - j -1)&1 != 0:
            for k in range(N):
                if k != j:
                    if (M[k] >> 2*N - j - 1)&1 != 0:
                        M[k] ^= M[j]
            S.append(j)
        else:
            for k in range(j, N):
                if (M[k] >> N - j - 1)&1 != 0:
                    M[k], M[j] = M[j], M[k]
                    break
                
            if (M[j] >> N - j - 1)&1 == 0:
                for e in M:
                    print(e)
                return False
            for k in range(N):
                if k != j:
                    if (M[k] >> N - j - 1)&1 == 1:
                        M[k] ^= M[j]
                M[j] = 0
                
    new = [M[z]&((1<<N)-1) for z in range(N)]

    d = [0]*N
    for index in S: d[index] = 1
    d_new = 0
    for i in range(len(d)-1):
        d_new ^= d[i]
        d_new <<= 1
        
    d_new ^= d[-1]
    
    return new, d_new
    
# The whole block lanczos algorithm
# See the README for sources
def block_lanczos(B, nb_relations, BLOCK_SIZE):
    Y = [random.randint(0, (1<<BLOCK_SIZE)-1) for _ in range(nb_relations)]

    X = [0]*nb_relations
    b = transpose_sparse(B, nb_relations)
    Vo = sparse_multiply(b, sparse_multiply(B, Y))
    i = 0
    
    P = [0 for _ in range(nb_relations)]
    V = Vo
    d = 1
    while d and i <= int(len(B)/(BLOCK_SIZE-0.764))+10:
        Z = sparse_multiply(b, sparse_multiply(B, V))
        vAv = dense_multiply(transpose_dense(V, BLOCK_SIZE), Z)
        vAAv = dense_multiply(transpose_dense(Z, BLOCK_SIZE), Z)
        
        W_inv, d = block(vAv, BLOCK_SIZE)
        
        X = add_vector(X, dense_multiply(V, dense_multiply(W_inv, dense_multiply(transpose_dense(V, BLOCK_SIZE), Vo))))
        
        neg_d = switch_indices(d)
        
        c = dense_multiply(W_inv, add_vector(multiply_d(vAAv, d), multiply_d(vAv, neg_d)))
        
        tmp1 = multiply_d(Z, d)
        tmp2 = multiply_d(V, neg_d)
        tmp3 = dense_multiply(V, c)
        tmp4 = multiply_d(vAv, d)
        tmp5 = dense_multiply(P, tmp4)
        tmp6 = dense_multiply(V, W_inv)
        tmp7 = multiply_d(P, neg_d)
            
        V, P = add_vector(add_vector(add_vector(tmp1, tmp2), tmp3), tmp5), add_vector(tmp6, tmp7)
        
        i += 1
 
    print("lanczos halted after "+str(i)+" iterations")
    x = add_vector(X, Y)
    Z = concatenate(x, V, BLOCK_SIZE)
    matrix = transpose_dense(sparse_multiply(B, Z), BLOCK_SIZE<<1)
    Z = transpose_dense(Z, BLOCK_SIZE<<1)
    matrix, Z = solve(matrix, Z, len(B))

    solutions = []
    for i in range(len(matrix)):
        if matrix[i] == 0 and Z[i] != 0 and Z[i] not in solutions:
            solutions.append(Z[i])

    if len(solutions) == 0:
        solutions = block_lanczos(B, nb_relations, BLOCK_SIZE<<1)

    return solutions
    
# Performs gaussian elimination
def solve(matrix, block, nb_relations):
    k = 0
    
    for l in range(nb_relations):
        for i in range(k, len(matrix)):
            if (matrix[i] >> nb_relations - i - 1)&1 != 0:
                matrix[k], matrix[i] = matrix[i], matrix[k]
                block[k], block[i] = block[i], block[k]
                k += 1
                break
                
        for z in range(i+1, len(matrix)):
            if (matrix[z] >> nb_relations - l - 1)&1 == 1:
                matrix[z] ^= matrix[k-1]
                block[z] ^= block[k-1]

    return matrix, block

# Computes a^(-1) (mod m)
def invmod(a, m):
    (r, u, R, U) = (a, 1, m, 0)
    while R:
        q = r//R
        (r, u, R, U) = (R, U, r - q *R, u - q*U)
    return u%m
    
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
        
# Compute x such that x² = n (mod p)
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
        
# Use chinese remainder theorem to generate some sieve polynomial coefficients
def CRT(moduli, a ,second_part):
    sum = 0
    for i in range(len(moduli)):
        sum = (sum+moduli[i]*second_part[i])%a
    return sum      
    
def is_prime(n):
    if n == 2 or n == 3:
        return True
    if not n&1 or n==1:
        return False
    r, s = 0, n - 1
    while not s&1:
        r += 1
        s >>= 1
    for _ in range(0, 50):
        a = random.randrange(2, n - 1)
        x = pow(a, s, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(0, r - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                break
        else:
            return False
    return True
    
# Fast but not always true Fermat primality test
# There are very few known Fermat primes, so it is most likely almost always correct
# Returns true is the number is not prime
def fermat_primality(n): return pow(2, n-1, n)-1


def fac(n):
    res = 1
    for i in range(2, n+1): res *= i
    return res
    
def binom(k, n):
    res = 1

    for i in range(n-k+1, n+1): res *= i
    for i in range(2, k+1): res //= i

    return res
        
## polynomial operations functions

def poly_div_bin(poly_a, poly_b):
    remainder = poly_a
    quotient = 0
    deg_a = poly_a.bit_length() - 1
    deg_b = poly_b.bit_length() - 1

    for j in reversed(range(deg_a - deg_b + 1)):
        if remainder & (1 << (deg_b + j)):
            quotient |= (1 << j)
            remainder ^= (poly_b << j)

    return quotient, remainder

def poly_prod_bin(poly_a, poly_b):
    res = 0

    for i in range(poly_a.bit_length()-1, 0, -1):
        if poly_a >> i & 1: res ^= poly_b
        res <<= 1

    if poly_a & 1: res ^= poly_b

    return res
    
def poly_gcd_bin(poly_a, poly_b):
    A = poly_a
    B = poly_b
    
    while B != 0:
        (Q, R) = poly_div_bin(A, B)
        A, B = B, R
        
    return A
    

## linear algebra functions
    
def add_vector(a, b):
    new = [a[i]^b[i] for i in range(len(a))]
    return new
    
def sparse_dot_prod(lbd, x):
    res = 0
    for i in lbd: res ^= x[i]
    
    return res
    
def sparse_multiply(A, b):
    new = [0]*len(b)
    for i in range(len(A)):
        tmp = 0
        for j in A[i]: tmp ^= b[j]
        new[i] = tmp
    return new
    
def dense_multiply(A, B):
    C = [0]*len(A)
    for i in range(len(A)):
        tmp = 0
        for j in range(len(B)):
            if (A[i] >> len(B)-j-1)&1:
                    tmp ^= B[j]
        C[i] = tmp
    return C
    
def transpose_dense(A, N):
    dim1_len = N
    B = [0]*dim1_len
    
    for i in range(len(A)-1):
        tmp = A[i]
        for j in range(dim1_len-1, -1, -1):
            if tmp&1: B[j] ^= 1
            B[j] <<= 1
            tmp >>= 1
            
    tmp = A[-1]
    for j in range(dim1_len-1, -1, -1):
        if tmp&1: B[j] ^= 1
        tmp >>= 1
        
    return B
    
def transpose_sparse(A, n):
    A_transpose = [set() for _ in range(n)]
    
    for i in range(len(A)):
        for index in A[i]:
            A_transpose[index].add(i)
        
    return A_transpose
    
def identity(n):
    new = [(1<<i) for i in range(n-1, -1, -1)]
    return new
    
def concatenate(A, B, N):
    return [(A[i]<<N)^B[i] for i in range(len(A))]


## Statistics function

def central(k,x): return pow(x, k/2-1)/(math.exp(x/2)*pow(2, k/2)*math.gamma(k/2))

def non_central(k, l ,x):
    if x <= 0: return 0
    res = 0
    for i in range(100): res += central(k+2*i, x)*pow(l/2, i)/fac(i)
    
    return res*math.exp(-l/2)
        
        
## time formating functions

def format_duration(delta):
    hours = delta.seconds//3600
    minutes = delta.seconds//60 - hours*60
    seconds = delta.seconds%60
    return str(delta.days)+" days, "+str(hours)+" hours, "+str(minutes)+" minutes, "+str(seconds)+" seconds"

def sieve(length, f_x, rational_base, algebraic_base, b, logs, leading,div,x,mod,n):
    ###TO DO: FIX LEADING COEFF AND ALL THAT
    ###WIP
    tmp_poly = new_coeffs(f_x, b)
    sieve_len=length
    pairs=[]
    x=(x*b)%mod
    sieve_array = [0]*sieve_len
    sieve_array_neg = [0]*sieve_len
    y=f_x[1]
    for q, p in enumerate(rational_base):
        if mod%p ==0:
            continue
        log = logs[q]

        for r in algebraic_base[q]:
            r=(r*b)%p
            
            s=solve_lin_con(mod,r-x,p)
            root=x+s*mod
            if (root**2+(y*b)*root-n*b**2)%(mod*p)!=0:
                print("fataaaaal")

            i=s
            while i < len(sieve_array):
                sieve_array[s]+=log

                i+=p
            i=-s%p
            while i < len(sieve_array):
                sieve_array_neg[s]+=log

                i+=p
    i=0
    while i < len(sieve_array):
        if sieve_array[i] > keysize*0.7: ##TO DO: FIX THRESHOLD
            a=x+i*mod
            eval1=a**2+(y*b)*a-n*(b**2)
            if eval1%mod !=0:
                print("error")
            eval2=a+(y*b)


            if eval1!=0 and eval2!=0  and math.gcd(a, b) == 1:
                pairs.append([[-b, a*leading], eval1, [[1 ,1]], eval2, [1], [-b,a], 1])


        i+=1

    i=0
    while i < len(sieve_array_neg):
        if sieve_array_neg[i] > keysize*0.7: ##TO DO: FIX THRESHOLD
            a=x-i*mod
            eval1=a**2+(y*b)*a-n*(b**2)
            if eval1%mod !=0:
                print("error")
            eval2=a+(y*b)



            if eval1!=0 and eval2!=0  and math.gcd(a, b) == 1:
                pairs.append([[-b, a*leading], eval1, [[1 ,1]], eval2, [1], [-b,a], 1])


        i+=1   
    
    return pairs,tmp_poly

def square_root(f, N, p, m0, m1, leading, bound):
 #   print("Square_root info f: "+str(f)+" N: "+str(N)+" p: "+str(p)+" m0: "+str(m0)+" m1: "+str(m1)+" leading: "+str(leading)+" bound: "+str(bound))
    d = len(f)-1
    root = compute_root_mod_Newton(N, f, p)
    print("Initial root computed, lifting...")
    root = Newton(root, N, f, p, bound)
    print("Lift done")
    return eval_F(leading*m0, m1, root, d-1)
                
def Newton(root, gamma, f, p, bound):
    modulo = p
    old = [1]
    while modulo < bound<<1:
        modulo *= modulo
        tmp = div_poly_mod(poly_prod(root, root), f, modulo)
        tmp = div_poly_mod(poly_prod(tmp, gamma), f, modulo)
        tmp = [-i%modulo for i in tmp]
        tmp[-1] = (tmp[-1]+3)%modulo
        tmp = div_poly_mod(poly_prod(root, tmp), f, modulo)
        tmp2 = invmod(2, modulo)
        root = [tmp2*i%modulo for i in tmp]
        
        tmp = div_poly_mod(poly_prod(root, gamma), f, modulo)
        for i in range(len(tmp)):
            if tmp[i] > modulo>>1: tmp[i] -= modulo

    return tmp
        
                
def compute_root_mod_Newton(number, f, p):
    n = div_poly_mod(number, f, p)
    for k in range(len(n)):
        if n[k]:
            n = n[k:]
            break
            
    s = pow(p, len(f)-1)-1
    r = 0
    while not s&1:
        s >>= 1
        r += 1
        
    z = [random.randint(0, p-1) for i in range(len(f)-1)]
    while quadratic_residue(z, f, p) != p-1:
        z = [random.randint(0, p-1) for i in range(len(f)-1)]
    
    lbd = power(n, f, p, s)
    omega = power(n, f, p, (s+1)>>1)
    zeta = power(z, f, p, s)
    while True:
        print("lbd: ",lbd)
        if lbd == [0]*len(lbd): return [0]
        if lbd == [1]: return power(omega, f, p,pow(p, len(f)-1)-2)
        k = 0
        for m in range(1, r):
            k = m
            if div_poly_mod(power(lbd, f, p, 1<<m), f, p)==[1]: break
            
        tmp = 1<<(r-k-1)
        tmp = power(zeta, f, p, tmp)
        omega = div_poly_mod(poly_prod(omega, tmp), f, p)
        
        tmp = 1<<(r-k)
        tmp = power(zeta, f, p, tmp)
        lbd = div_poly_mod(poly_prod(lbd, tmp), f, p)
        
def compute_root_mod(number, f, p):
    n = div_poly_mod(number, f, p)
    for k in range(len(n)):
        if n[k]:
            n = n[k:]
            break
            
    s = pow(p, len(f)-1)-1
    r = 0
    while not s&1:
        s >>= 1
        r += 1
        
    z = [random.randint(0, p-1) for i in range(len(f)-1)]
    while quadratic_residue(z, f, p) != p-1:
        z = [random.randint(0, p-1) for i in range(len(f)-1)]
    
    lbd = power(n, f, p, s)
    omega = power(n, f, p, (s+1)>>1)
    zeta = power(z, f, p, s)
    while True:
        if lbd == [0]*len(lbd): return [0]
        if lbd == [1]: return omega
        k = 0
        for m in range(1, r):
            k = m
            if div_poly_mod(power(lbd, f, p, 1<<m), f, p)==[1]: break
            
        tmp = 1<<(r-k-1)
        tmp = power(zeta, f, p, tmp)
        omega = div_poly_mod(poly_prod(omega, tmp), f, p)
        
        tmp = 1<<(r-k)
        tmp = power(zeta, f, p, tmp)
        lbd = div_poly_mod(poly_prod(lbd, tmp), f, p)

# Used when the gaussian pivot is used, convert sparse vector of indices into dense binary vector
def convert_to_binary(z, n):
    res = [0]*n
    for index in z: res[index] ^= 1
    return res
    
def convert_to_binary_lanczos(z, n):
    res = [0]*n
    for i in range(n):
        if (z >> n - i - 1)&1: res[i] = 1
    return res
    
def create_solution(pairs, null_space, n, len_primes, primes, f_x, m0, m1, inert, f_prime_sq, leading, f_prime_eval, u,qua):
    #print("n: "+str(n)+" len_primes: "+str(len_primes)+" f_x: "+str(f_x)+" m0: "+str(m0)+" m1:"+str(m1)+" inert: "+str(inert)+" f_prime_sq: "+str(f_prime_sq)+" leading: "+str(leading)+" f_prime_eval: "+str(f_prime_eval)+" u: "+str(u))

    ###Pairs used: index 6, index 4, index 3,index 0
   # print("len nspace: "+str(len(null_space))+" len pairs: "+str(len(pairs)))
   # i=0
   # while i < len(pairs):
        #[i][1]=1

       # i+=1
    xtotal=1
    int_square=1
    f_norm = 0
    tmp = 1
    for x in f_x:
        f_norm += x*x*tmp
        tmp *= leading
    f_norm = int(math.sqrt(f_norm))+1
    fd = int(pow(len(f_x)-1, 1.5))+1
    S = 0
    
    x = f_prime_eval
    rat=create_rational(null_space, n, len_primes, primes, pairs)
    x = x*rat%n
    hit=0
    rational_square = [i for i in f_prime_sq]


    for k in range(len(null_space)):
        if null_space[k]:
           

            
            ###pair info: [[-b, a*leading], eval1, [[1, 1]], eval2, [1], [-b,a], 1]
            ###If free relation: [[leading_coeff*p], leading_coeff*p, [[1,1]], p*m1, [1], [0,p], 1]
         #   print("adding pair[k] [-b, a*leading]: "+str(pairs[k][0])+" eval1: "+str(pairs[k][1])+" eval2: "+str(pairs[k][3])+" [-b,a]: "+str(pairs[k][5]))



            prod=poly_prod(rational_square, pairs[k][0])
            rational_square = div_poly(prod, f_x)
            S += pairs[k][6] ##exp?
    




    lead=pow(leading, S>>1, n)  
    x = x*lead%n
   # print("Rational info f_prime_eval: "+str(f_prime_eval)+" rat: "+str(rat)+" lead: "+str(lead)+" final: "+str(x)+" nullspace len: "+str(hit)+" S: "+str(S)+" inert: "+str(inert)+" u: "+str(u))          
    coeff_bound = [fd*pow(f_norm, len(f_x)-1-i)*pow(2*(leading*u)*f_norm, S>>1) for i in range(len(f_x)-1)]
    #print(" f_prime_sq: "+str(f_prime_sq)+" prod: "+str(prod))
    y = square_root(f_x, rational_square, inert, m0, m1, leading, max(coeff_bound))
    y = y*pow(m1, S>>1, n)%n
  
    
    return x, y
    
def create_rational(null_space, n, len_primes, primes, pairs):
    x = 1
    vec = [0]*(len_primes)

    for z in range(len(null_space)):
        if null_space[z]:
            for large_prime in pairs[z][4]: x = x*large_prime%n
            cmp = pairs[z][3]
            for p in range(len_primes):
                tmp = primes[p]
                tmp2 = 0
                while not cmp%tmp:
                    tmp *= primes[p]
                    tmp2 += 1
                vec[p] += tmp2

    for i in range(len_primes): 
        x = x*pow(primes[i], vec[i]>>1, n)%n

    return x
    


def compute_factors(pairs_used, vec, n, primes, g, g_prime, g_prime_sq, g_prime_eval, m0, m1, leading_coeff,
                    inert_set, M, time_1,Q):

    x, y = create_solution(pairs_used, vec, n, len(primes), primes, g, m0, m1, inert_set[-1], g_prime_sq,
                              leading_coeff, g_prime_eval, M<<1,Q)
    print("x**2%n: "+str(x**2%n)+" y**2%n: "+str(y**2%n))
    if x != y and math.gcd(x-y, n) != 1 and math.gcd(x+y, n) != 1:
        print_final_message(x, y, n, time_1)
        return True, 1

    return False, 0
    
def print_final_message(x, y, n, time_1):
    time_2 = datetime.now()

    print("found factor : "+str(math.gcd(x-y,n)))
    print("found factor : "+str(math.gcd(x+y,n)))
    print("null space found in "+format_duration(time_2-time_1)+".\n")
    sys.exit()
# Create the sparse matrix:
# Consider the dense binary matrix A[i][j]
# The sparse matrix M is defined by M[i] contains j if and only if A[i][j] = 1
# For large numbers to factor, this reduces that memory requirement for the linear algebra part
# Each prime in the factor base corresponds to one line
# Each relation corresponds to one column
def build_sparse_matrix(pairs, rat, alg, qua, div_lead):
    matrix = []
    rc=0
    while rc < len(pairs):
        r=pairs[rc]
  #  for r in pairs:
        index = 1
        
        if r[3] < 0: 
            line = [0]
        else: line = []

        for p in rat:
            tmp = p
            tmp2 = 0
            while not r[3]%tmp:
                tmp *= p
                tmp2 ^= 1
            if tmp2: 
                line.append(index)
            index += 1
            
        if r[1] < 0: 
            line.append(index)
        index += 1

        for p in range(len(alg)):
            tmp = rat[p]
            tmp2 = 0
            while not r[1]%tmp:
                tmp *= rat[p]
                tmp2 ^= 1
           # if rc == 0:
             #   print("r[1]: "+str(r[1])+" tmp: "+str(tmp)+" tmp2: "+str(tmp2))
            for i in range(len(alg[p])):
                if not eval_mod(r[5], alg[p][i], rat[p]):
                 #   if rc == 0:
                     #   print("eval_mod(r[5]: "+str(r[5])+" alg[p][i]: "+str(alg[p][i])+" rat[p]: "+str(rat[p])+" eval_mod(r[5], alg[p][i], rat[p]): "+str(eval_mod(r[5], alg[p][i], rat[p])))
 
                    if tmp2: 
                        line.append(index)
                index += 1
            
        for p in range(len(qua)):
            if compute_legendre_character(eval_mod(r[5], qua[p][1], qua[p][0]), qua[p][0]) == -1: 
                line.append(index)
            index += 1
       # print("div_lead: "+str(div_lead))
        for p in range(len(div_lead)):
            if r[7+p]:
                tmp = div_lead[p]
                tmp2 = 0

                while not r[1]%tmp:
                    tmp *= div_lead[p]
                    tmp2 ^= 1
                if tmp2: 
                    line.append(index)

            index += 1
            
        if r[6]&1: 
            line.append(index)
        matrix.append(set(line))
        rc+=1

    return matrix
    
# Reduce the sparse matrix, according to various rules:
# 1. If a line contains only one non-zero value, then it is deleted as well as the column containing the corresponding non-zero value
# 2. If a line contrains no non-zero value, then it is deleted
# 3. We only keep 10 more columns than lines, to ensure we still have solutions while reducing the matrix size
def reduce_sparse_matrix(matrix, pairs):
    flag = True
    while flag:
        flag = False

        singleton_queue = [index for index, row in enumerate(matrix) if len(row) == 1]
        active_cols = set(range(len(pairs)))

        flag = len(singleton_queue)

        while singleton_queue:
            index = singleton_queue.pop()

            if len(matrix[index]) != 1: continue
            coeff = next(iter(matrix[index]))

            for j, row in enumerate(matrix):
                if j != index and coeff in row:
                        row.remove(coeff)
                        if len(row) == 1: singleton_queue.append(j)

            matrix[index].clear()
            active_cols.discard(coeff)

        matrix = [row for row in matrix if row]

        pairs = [pairs[index] for index in sorted(active_cols)]

        mapping = {old: new for new, old in enumerate(sorted(active_cols))}
        matrix = [{mapping[c] for c in row if c in mapping} for row in matrix]

        active_cols = set(range(len(pairs)))
        target_n_cols = len(matrix)+10
        to_delete = len(pairs) - target_n_cols
        current_len_row = 2

        while to_delete >= current_len_row-1:

            index = -1
            for i in range(len(matrix)):
                if len(matrix[i]) == current_len_row:
                    flag = True
                    index = i
                    break

            if index != -1:
                for coeff in matrix[index]:
                    active_cols.discard(coeff)

                    for j, row in enumerate(matrix):
                        if j != index and coeff in row:
                            row.discard(coeff)

                to_delete -= current_len_row - 1

                matrix[index].clear()

            else:
                current_len_row += 1

        matrix = [row for row in matrix if row]

        pairs = [pairs[index] for index in sorted(active_cols)]

        mapping = {old: new for new, old in enumerate(sorted(active_cols))}
        matrix = [{mapping[c] for c in row if c in mapping} for row in matrix]

    return matrix, pairs
                        
## Creation of polynomials
def find_relations(f_x, leading_coeff, g, primes, R_p, Q, B_prime, divide_leading, prod_primes, pow_div, pairs_used,const1, const2, logs, M,mult,x,mod,n):

    full_found = 0 
    b = mult
  
    while b < 20:
        
        div = 1
        for q in range(len(divide_leading)):
            p = divide_leading[q]
            if not b%p:
                tmp = p
                while not b%tmp and tmp <= pow_div[q]:
                    div *= p
                    tmp *= p
        
        pairs,tmp_poly = sieve(M, f_x, primes, R_p, b, logs, leading_coeff, div,x,mod,n) ##to do: returning tmp_poly for debug purposes.Remove later
        for i, pair in enumerate(pairs):
            z = trial(pair, primes, const1, const2, div)
            if z[0]:
                tmp = [u for u in pair]
                for p in divide_leading: 
                    tmp.append(not pair[5][0]%p)
                if z[2] == 1 and z[4] == 1:
                  #  print("adding: "+str(tmp))
                    pairs_used.append(tmp)
                    full_found += 1
        b+=1
                                                                                        #cycle_len, full_found, partial_found_pf)                  
    return pairs_used

def coefficients(n, m, d):
    tmp1, tmp2 = n, 1
    coeff = [0]*(d+1)
    i = 0
    for k in range(d+1):
        tmp3 = tmp1%(tmp2*m)
        coeff_k = (tmp3//tmp2)
        coeff[d-k] = coeff_k
        tmp1 -= tmp3
        tmp2 *= m
    return coeff

def get_derivative(f):
    res = [0]*(len(f)-1)
    for i in range(len(f)-1):
        res[i] = (len(f)-1-i)*f[i]
    return res
    
def new_coeffs(f, x):
    b = 1
    tmp = [i for i in f]
    for i in range(len(f)-1):
        tmp[i] *= b
        b *= x
    tmp[-1] *= b
    return tmp
    
def power(poly, f, p, exp):
    if exp == 1: return poly

    tmp = power(poly, f, p, exp>>1)
    tmp = poly_prod(tmp, tmp)

    if exp&1:
        tmp = poly_prod(tmp, poly)
        return div_poly_mod(tmp, f, p)
    
    else: return div_poly_mod(tmp, f, p)
    
def shift(poly,k):
    res = [i for i in poly]
    for i in range(len(poly)-1):
        for j in range(i+1, len(poly)):
            res[j] += binom(j-i, len(poly)-1-i)*pow(k, j-i)*poly[i]
    return res
    
## Characteristics of one polynomial
    
def irreducibility(f, p):
    g = [1,0]
  #  print("start: "+str(p)+" poly: "+str(f))
    for i in range((len(f)-1)//2):
        
        g = power(g, f, p, p)
        tmp2 = [u for u in g]
        if len(tmp2) == 1: 
            tmp2 = [-1, tmp2[0]]
        else: 
            tmp2[-2] -= 1
        tmp = gcd_mod(f, tmp2, p)
      #  print("tmp: "+str(tmp)+" range: "+str((len(f)-1)//2+1))
        if len(tmp) > 1: 
            return False
    return True
    


def fast_roots(f, p):
    tmp_f = [i%p for i in f]
    for k in range(len(tmp_f)):
        if tmp_f[k]:
            tmp_f = tmp_f[k:]
            break

    roots = []
    tmp = [0]*len(f)
    for i in range(len(f)): tmp[i] = eval_mod(tmp_f, i, p)
    for q in range(1, len(tmp_f)):
        for k in range(len(tmp_f)-1, q-1, -1): tmp[k] -= tmp[k-1]

    res = tmp[0]
    if not res: roots.append(0)
    for q in range(1, p):
        res = (res+tmp[1])%p
        if not res:
            roots.append(q)
            if len(roots) == len(f)-1: return roots
        for k in range(1, len(tmp)-1): tmp[k] += tmp[k+1]
    return roots
    
    
def gcd_mod(f, poly, p):
    while poly != [0]*len(poly):
        (f,poly) = (poly, div_poly_mod(f, poly, p))

    return f
    
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
    
def roots(g, p):
    if len(g) == 1: return []
    if len(g) == 2: return [-g[1]*invmod(g[0], p)%p]
    if len(g) == 3:
        tmp = (g[1]*g[1]-4*g[0]*g[2])%p
        if tmp == 0: return [-g[1]*invmod(2*g[0], p)%p]
        if compute_legendre_character(tmp, p) == -1: return []
        tmp = compute_sqrt_mod_p(tmp, p)*invmod(2*g[0], p)%p
        return [(-g[1]*invmod(g[0]<<1, p)+tmp)%p, (-g[1]*invmod(g[0]<<1, p)-tmp)%p]
    
    h = [1]
    while len(h) == 1 or h == g:
        a = random.randint(0, p-1)
        h = power([1, a], g, p, (p-1)>>1)
        for k in range(len(h)):
            if h[k]:
                h = h[k:]
                break
        h[-1] -= 1
        h = gcd_mod(h, g, p)
    r = roots(h, p)
    h = quotient_poly_mod(g, h, p)
    return r+roots(h, p)
    
def quadratic_residue(poly, f, p):
    return power(poly, f, p, (pow(p, len(f)-1)-1)>>1)[0]
    

## Operations between two polynomials
    
def poly_prod(a, b):
    res = [0]*(max(len(a), len(b))+min(len(a), len(b))-1)

    for i in range(len(a)):
        for j in range(len(b)):
            res[i+j] += a[i]*b[j]

    return res
    
def poly_prod_mod(a, b, p):
    res = [0]*(max(len(a), len(b)) + min(len(a), len(b))-1)

    for i in range(len(a)):
        for j in range(len(b)):
            res[i+j] = (res[i+j]+a[i]*b[j])%p

    return res
    
def div_poly(a,b):
    remainder = [i for i in a]
    difference = len(a)-len(b)+1
    leading_coeff = b[0]

    for j in range(difference):
        quotient = -remainder[j]//leading_coeff
        for k in range(len(b)):
            remainder[j+k] += quotient*b[k]
            
    for k in range(len(remainder)):
        if remainder[k]: return remainder[k:]
        
    return [0]
    
def div_poly_mod(a, tmp_b, p):
    remainder = [i%p for i in a]
    b = [i%p for i in tmp_b]
    
    #print(remainder, b)
    while not b[0]: del b[0]
    
    difference = len(a)-len(b)+1
    coeff = invmod(-b[0], p)
    for j in range(difference):
        if remainder[j]:
            quotient = remainder[j]*coeff%p
            remainder[j] = 0
            for k in range(1,len(b)): remainder[j+k] = (remainder[j+k]+quotient*b[k]%p)%p
            
    for k in range(len(remainder)):
        if remainder[k]: return remainder[k:]
        
    return [0]
    
def quotient_poly_mod(a, b, p):
    remainder = [i%p for i in a]
    b = [i%p for i in b]
    
    while not b[0]: del b[0]
    
    difference = len(a)-len(b)+1
    coeff = invmod(-b[0], p)
    res = [0]*difference

    for j in range(difference):
        quotient = remainder[j]*coeff%p
        res[j] = -quotient
        for k in range(len(b)):
            remainder[j+k] = (remainder[j+k]+quotient*b[k]%p)%p
            
    for k in range(len(res)):
        if res[k]: return res[k:]
        
    return [0]
    
def gcd_mod(f, poly, p):
    while poly != [0]*len(poly):
        (f, poly) = (poly, div_poly_mod(f, poly, p))
    return f
    
## Evaluation of polynomials

def eval_mod(f, x, n):
    res = 0

    for i in range(len(f)-1):
        res += f[i]
        res *= x
        res %= n

    res += f[-1]
    res %= n

    return res
    
def evaluate(f, x):
    res = 0

    for i in range(len(f)-1):
        res += f[i]
        res *= x

    res += f[-1]

    return res
    
def eval_F(x, y, f, d):
    tmp = 0
    tmp2 = 1

    for k in range(d):
        tmp += f[k]*tmp2
        tmp *= x
        tmp2 *= y

    tmp += f[d]*tmp2

    return tmp

## NFS functions

    

    
def get_sieve_region(f, B):
    F = poly_prod(f, f)
    d = len(f)-1
    ratios = [math.log(1e-7+abs(f[i+1]/(1+abs(f[i])))) for i in range(d)]
    s = 2*int(math.exp(sum(ratios)/d)) # skew factor
    k = 1 # shift
    best_norm = None

    while k > 0:
        updated = False

        if s-k > 0:
            norm = get_Lnorm(F, s-k, B)
            if best_norm == None or norm < best_norm:
                s = s-k
                k <<= 1
                best_norm = norm
                updated = True

        if s+k < B:
            norm = get_Lnorm(F, s+k, B)
            if best_norm == None or norm < best_norm:
                s = s+k
                k <<= 1
                best_norm = norm
                updated = True

        if not updated: k >>= 1

    x = round(B*math.sqrt(s))

    return x, s
    
def get_Lnorm(F, s, B):
    sqrt = math.sqrt(s)
    n = len(F)
    
    base_X, base_Y = B*sqrt/2, 1+B/sqrt
    current_X, current_Y = pow(base_X, n), base_Y
    base_X, base_Y = base_X*base_X, base_Y*base_Y

    res = 0.0

    for i in range(0, n, 2):
        res += (2*current_X)*F[i]*(current_Y-1)/((i+1)*(n-i))
        current_X /= base_X
        current_Y *= base_Y

    return math.log(abs(res))/2
        
def initialize_2(f_x, d, primes, leading_coeff):
    m1=f_x[0]
    pairs_used, R_p, logs, tmp = [], [], [], 10*d
    divide_leading, pow_div = [], []
    for p in primes:
        if not leading_coeff%p:
            divide_leading.append(p)
            u = p
            while not leading_coeff%u: 
                u *= p
            pow_div.append(u//p)

        logs.append(round(math.log2(p)))
        if p > tmp: 
            R_p.append(find_roots_poly(f_x, p))
        else: 
            R_p.append(fast_roots(f_x, p))
        if len(R_p[-1]) == d:
            if d&1: 
                pairs_used.append([[leading_coeff*p], leading_coeff*pow(p, d), [[1,1]], p*m1, [1], [0,p], 1])
            else: 
                pairs_used.append([[leading_coeff*p], leading_coeff*p, [[1,1]], p*m1, [1], [0,p], 1])
               # print("free: ",pairs_used[-1])
    for p in divide_leading:
        for i in range(len(pairs_used)): 
            pairs_used[i].append(True)
        
    B_prime = sum([len(i) for i in R_p]) + 1
    
    return pairs_used, R_p, logs, divide_leading, pow_div, B_prime
    
def initialize_3(n, f_x, f_prime, const1, leading_coeff):
    k = 3*n.bit_length()
    Q = []
    q = const1+1
    if not q%2: q += 1
    while len(Q) < k:
        if is_prime(q) and leading_coeff%q:
            if not n%q: 
                q+=2
                continue
            else:
                tmp = find_roots_poly(f_x, q)
                for r in tmp:
                    if eval_mod(f_prime, r, q): Q.append([q, r])
        q += 2
        
    return Q, k
    
def coefficients(P):
    Pc = []
    for i in range(P.degree(), -1, -1):
        Pc.append(P.nth(i))
    return tuple(Pc)
def sylvester(P, Q):
    rows = []
    m = P.degree()
    n = Q.degree()
    size = m + n
    CP = coefficients(P)
    CQ = coefficients(Q)
    for i in range(size):
        tail = []
        row = []
        if i in range(0, n):
            row = list(CP)
            row.extend((n-1-i)*[0])
            row[:0] = [0]*i
            rows.append(row)
        if i in range(n, size):
            row = list(CQ)
            row.extend((size-1-i)*[0])
            row[:0] = [0]*(size-len(row))
            rows.append(row)
    return Matrix(rows)
def resultant(P, Q):
    return sylvester(P,Q).det()

def NFS_sieve(n,f_x,primeslist,mult,R_p,pairs_used,logs,pow_div,divide_leading,B_prime,x,mod):

    LARGE_PRIME_CONST=10000
    n = int(n)
    d=2
    B=primeslist[-1]
    primes =primeslist
    prod_primes = math.prod(primes)
    const1, const2 = LARGE_PRIME_CONST*primes[-1], LARGE_PRIME_CONST*primes[-1]*primes[-1]
    
    f_prime = get_derivative(f_x)
    M=lin_sieve_size
    leading_coeff = f_x[0]
    g = [1, f_x[1]]
    for i in range(2, len(f_x)): 
       g.append(f_x[i]*pow(leading_coeff, i-1))
    Q, k = initialize_3(n, f_x, f_prime, B, leading_coeff)
    pairs_used = find_relations(f_x, leading_coeff, g, primes, R_p, Q, B_prime, divide_leading,prod_primes, pow_div, pairs_used, const1, const2, logs,M,mult,x,mod,n)
    return pairs_used,R_p,Q,divide_leading,g,M

def NFS_solve(n,f_x,primeslist,pairs_used,R_p,Q,divide_leading,V,g,M,inert_set):
    m0=-f_x[1]
    m1=f_x[0]
    d=2
    leading_coeff = f_x[0]

    f_prime = get_derivative(f_x)
    BLOCK_SIZE=8


    g_prime = get_derivative(g)
    g_prime_sq,g_prime_eval = div_poly(poly_prod(g_prime, g_prime), g), pow(leading_coeff, d-2, n)*eval_F(m0, m1, f_prime, d-1)%n
    print("sieving complete, building matrix...")
    matrix = build_sparse_matrix(pairs_used, primeslist, R_p, Q, divide_leading)
    matrix = transpose_sparse(matrix, V)
    print("matrix built "+str(len(matrix))+"x"+str(len(pairs_used))+" reducing...")
    print("pairs before: "+str(len(pairs_used)))
    matrix, pairs_used = reduce_sparse_matrix(matrix, pairs_used)
    print("pairs after: "+str(len(pairs_used)))
    print("matrix built "+str(len(matrix))+"x"+str(len(pairs_used))+" finding kernel...")
    time_1 = datetime.now()
    while True:
        null_space = block_lanczos(matrix, len(pairs_used), BLOCK_SIZE)
                
        print(str(len(null_space))+" kernel vectors found")
        for vec in null_space:
            vec = convert_to_binary_lanczos(vec, len(pairs_used))
            flag, res = compute_factors(pairs_used, vec, n, primeslist, g, g_prime, g_prime_sq, g_prime_eval,
                                                             m0, m1, leading_coeff, inert_set, M,
                                                             time_1,Q)
            if flag: return res

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
            return -1
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
            return -1
        x1 = (y - b_div_2) % p
        x2 = (p - y - b_div_2) % p
        return [x1,x2]

def lift_root(z,y,n, root, p, exp):
    r =(z*root**2+root*y-n)%p**(exp+1)
    deriv = 2*z*root+y
    inv = modinv(deriv, p)
    t = ((-r//p**exp)*inv)%p
    new_root = (root+t*p**exp)%p**(exp+1)
    return new_root

def solve_roots(prime,n):
    hmap_p=[]
    k=0
    while k < prime and k < quad_sieve_size:
        hmap_p.append({})
        k+=1
   

    k=1
    while k < prime and k < quad_sieve_size+1:
        y=0
        while y < prime:
            roots=solve_quadratic_congruence(k, y, -n, prime)
            if roots == -1:
                y+=1
                continue
            y2=y
            while y2 < prime:
                for x in roots:
                    if (k*x**2+y2*x-n)%prime!=0: 
                        print("error")   
                    try:
                        ret=hmap_p[k-1][y2]
                        ret.append(x)

                    except Exception as e:
                        hmap_p[k-1][y2]=[x]
             

   
                y2+=prime
            y+=1
        k+=1
    return hmap_p

def create_hashmap(n,primeslist):
    i=0
    hmap=[]
    hmap2=[]
    while i < len(primeslist):
        hmap_p=solve_roots(primeslist[i],n)
        hmap.append(hmap_p)
        hmap2.append([])
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

@cython.boundscheck(False)
@cython.wraparound(False)
cdef factorise_fast_quads(value,long long [::1] factor_base):
    factors = set()
    if value % 2 == 0:
        factors ^= {2}
        value //= 2
        if value % 2 == 0:
            return -1, -1
    length=factor_base[0]#len(factor_base)#factor_base[0]
    cdef Py_ssize_t i=1
    while i < length:
        factor=factor_base[i]
        if value % factor == 0:
            factors ^= {factor}
            value //= factor
            if value % factor == 0:
                return -1, -1
        i+=1
    return factors, value

def filter(qbase,n,start,size):
    ###Note: We look for quadratic coefficients that factor over the factor base but have no even exponents. This garantuees unique results
    interval= np.zeros(size+1,dtype=np.int32)
    interval2= np.zeros(size+1,dtype=np.int32)

    length=qbase[0]
    i=1
    while i < length:
        prime=qbase[i]
        log=round(math.log2(prime))
        dist=start%prime
        dist=-dist%prime
        interval[dist::prime]+=log
        if (start+dist)%prime !=0:
            print("fatal error")
        i+=1
    
    
    valid_quads=[]
    valid_quads_factors=[]

    roots=[]
    i=0
    while i < size:
        if start+i > 999:
            threshold=round(math.log2((start+i)*0.50))
        if start+i < 1000 or interval[i]>threshold:
            quad_local_factors, quad_value = factorise_fast_quads(start+i,qbase) 
         #   print("checing")
            if quad_value != 1:
                i+=1
                continue

            interval2[i]=1
         #   print("found")
            valid_quads.append(start+i)
            valid_quads_factors.append(quad_local_factors)
        
        i+=1
    return valid_quads,valid_quads_factors,interval2

def solve_quad_integers(a,b,c):
    ##To do: Could I use this for the sieving process?
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
#@cython.boundscheck(False)
#@cython.wraparound(False)
def generate_modulus(n,primeslist,seen,tnum,close_range,too_close,LOWER_BOUND_SIQS,UPPER_BOUND_SIQS,tnum_bit,quad):
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
        while not found_a_factor and counter3< const_2 and randindex <base:
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

def solve_lin_con(a,b,m):
    ##ax=b mod m
    #g=gcd(a,m)
    #a,b,m = a//g,b//g,m//g
    return pow(a,-1,m)*b%m  

def find_xy(primeslist,n,x,collected,z,lin):
    xy=x+lin

    i=0
    while i < len(collected):
        
        prime=collected[i][0]
        i+=1
    return xy

def retrieve(hmap,primeslist,x,lin,n):

    collected=[1]*(quad_sieve_size+1)
    mult_list=[]
    interval=array.array("i",[0]*lin_sieve_size)
    i=0
    while i < quad_sieve_size:
        mult_list.append(array.array("i",[0]*lin_sieve_size))
        i+=1
    i=0
    while i < len(hmap):
        
        prime=primeslist[i]#**2
        log=round(math.log2(prime))
        c=1
        while c <  len(hmap[i])+1:
            z=c
            xy=(z*x)+lin
            if (z*x)%prime !=0 and xy%prime !=0:
                mul=solve_lin_con(xy,hmap[i][c-1][(x)%prime],prime)
                t=mul   
                y=(xy*mul)-x*z
                if (z*x**2+y*x-n)%prime !=0:
                    print("fatal error: "+str(z)+" p: "+str(prime)+" x: "+str(x))
                    sys.exit()
                while t < lin_sieve_size:
                    mult_list[c-1][t]+=log
                    t+=prime
       
              #  mult_list[c].append([prime,mul])

             #   collected.extend([prime,c])
            c+=1
        i+=1

    return mult_list

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

def coefficients(P):
    Pc = []
    for i in range(P.degree(), -1, -1):
        Pc.append(P.nth(i))
    return tuple(Pc)
def sylvester(P, Q):
    rows = []
    m = P.degree()
    n = Q.degree()
    size = m + n
    CP = coefficients(P)
    CQ = coefficients(Q)
    for i in range(size):
        tail = []
        row = []
        if i in range(0, n):
            row = list(CP)
            row.extend((n-1-i)*[0])
            row[:0] = [0]*i
            rows.append(row)
        if i in range(n, size):
            row = list(CQ)
            row.extend((size-1-i)*[0])
            row[:0] = [0]*(size-len(row))
            rows.append(row)
    return Matrix(rows)
def resultant(P, Q):
    return sylvester(P,Q).det()

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

cdef construct_interval(list ret_array,partials,n,primeslist,hmap,hmap2,large_prime_bound,primeslist2,small_primeslist):
    i=0
  #  while i < len(hmap):
      #  print("prime: "+str(primeslist[i])+" "+str(hmap[i]))
      #  i+=1

    grays = get_gray_code(20)
    primelist_f=copy.copy(primeslist)
    primelist_f.insert(0,len(primelist_f)+1)
    primelist_f=array.array('q',primelist_f)

    primelist2_f=copy.copy(primeslist2)
    primelist2_f.insert(0,len(primelist2_f)+1)
    primelist2_f=array.array('q',primelist2_f)
    cdef Py_ssize_t size
    primelist=copy.copy(primeslist)
    primelist.insert(0,2) ##To do: remove when we fix lifting for powers of 2
    #primelist.insert(0,-1)
    print("[i]Filtering Quadratic Coefficients (quad_size) (to do: can be saved to disk for re-use)")
    valid_quads,valid_quads_factors,qival=filter(primelist_f,n,1,quad_sieve_size+1)
    print("[i]Filtering interval indices (lin_size) (to do: can be saved to disk for re-use)")
    start=1#round(n**0.25)

    valid_ind,valid_ind_factors,lival=filter(primelist_f,n,1,lin_sieve_size+1)
    print("[i]Entering attack loop")
    smooths=[]
    coefficients=[]
    factors=[]
    jsymbols=[]
    testl=[]
    disc1_squared_list=[]
    xy_list=[]
    xy_f_list=[]
    x_list=[]
    x_f_list=[]
    x1=1
    y0=1
    o=1
    close_range=20
    too_close=1
    LOWER_BOUND_SIQS=1
    UPPER_BOUND_SIQS=4000
    tnum=int(((n)**0.3) /1)
    
    threshold = int(math.log2((lin_sieve_size)*math.sqrt(abs(n))) - thresvar)
    if threshold < 0:
        threshold = 1
    all_pairs_used=[]

    divide_leading=[]


    z=1
    y_start=1#round(n**0.2)
    y_ind=0
    while y_ind<100_000_000:  ##This is shit... fix this.. 
        y=y_start+y_ind
        f_x=[1,y,-n]
        seen=[]

        d=2
        pairs_used, R_p, logs, divide_leading, pow_div, B_prime = initialize_2(f_x, d, primelist, z)

        while 1:
            new_mod,cfact,indexes=generate_modulus(n,primeslist,seen,tnum,close_range,too_close,LOWER_BOUND_SIQS,UPPER_BOUND_SIQS,bitlen(tnum),z)
            if new_mod ==0:
                print("exiting, no new modulus found")
                break
            collected=[]
            e=0
            while e < len(indexes):
            
                ind=indexes[e]
                prime=primeslist[ind]
                l=hmap[ind][z-1]
                collected.append(l)
         #   print("prime: "+str(prime)+" l: "+str(l))
                e+=1

            t=0
            hit=0
            x_list=[]
 
            mod=1
            while t < len(collected):
                prime=primeslist[indexes[t]]
                try:
                    ret=collected[t][y%prime]
                    x_list.extend([prime,ret])

                    mod*=prime
                except Exception as e:
                    hit=1
                    break

                t+=1
            if hit == 1:
                continue  

        #    print(" y: "+str(y)+" x_list: "+str(x_list))
            enum=[]
            x_list=get_partials(mod,x_list)
            t=0
            while t < len(x_list):
                enum.append(x_list[t+1])
                t+=2
            root_list=[]
            mult=1
            for comb in itertools.product(*enum):
                root=0
                for l in comb:
                    root+=l 
                root%=mod
                if (z*root**2+y*root-n)%mod !=0:
                    print("fatal error")
                root_list.append(root)

            
            for x in root_list:
                pairs_used,R_p,Q,divide_leading,g,M=NFS_sieve(n,f_x,primelist,mult,R_p,pairs_used,logs,pow_div,divide_leading,B_prime,x,mod)
                find = 3 + len(primeslist) + B_prime + len(Q) + len(divide_leading)
                print("found: "+str(len(pairs_used))+" / "+str(find))     
                if len(pairs_used)>find+10:
                    inert_set = []
                    m1=f_x[0]
                    leading_coeff=f_x[0]
                    for p in primeslist:
                        if m1%p and leading_coeff%p and irreducibility(g, p): 
            
                            inert_set.append(p)
                    if len(inert_set)==0:
                        continue
                    NFS_solve(n,f_x,primelist,pairs_used,R_p,Q,divide_leading,find,g,M,inert_set)



        y_ind+=1




   
    return      

def find_comb(n,hmap,hmap2,primeslist1,primeslist2,small_primeslist):

    ret_array=[[],[],[],[]]

    

    partials={}


    large_prime_bound = primeslist1[-1] ** lp_multiplier
    construct_interval(ret_array,partials,n,primeslist1,hmap,hmap2,large_prime_bound,primeslist2,small_primeslist)

    return 0

def get_primes(start,stop):
    return list(sympy.sieve.primerange(start,stop))

@cython.profile(False)
def main(l_keysize,l_workers,l_debug,l_base,l_key,l_lin_sieve_size,l_quad_sieve_size):
    global key,keysize,workers,g_debug,base,key,quad_sieve_size,small_base,lin_sieve_size,qbase
    key,keysize,workers,g_debug,base,quad_sieve_size,lin_sieve_size=l_key,l_keysize,l_workers,l_debug,l_base,l_quad_sieve_size,l_lin_sieve_size
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
    
    while len(primeslist2) < qbase:
        if n%primeslist[i]==0:
            i+=1
            continue
        primeslist2.append(primeslist[i])
        i+=1
    launch(n,primeslist1,primeslist2,small_primeslist)     
    duration = default_timer() - start
    print("\nFactorization in total took: "+str(duration))