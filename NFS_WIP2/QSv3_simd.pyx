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




def square_root(f, N, p, m0, m1, leading, bound):
    d = len(f)-1
    root = compute_root_mod_Newton(N, f, p)
    if root == -1:
        return -1
    root = Newton(root, N, f, p, bound)
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
    count=0
    while True:
        count+=1
        if count > 50:
            print("square root stuck.. something probably went wrong")
            return -1
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




    
def create_solution(pairs, null_space, n, len_primes, primes, m0, m1, leading, u,qua,f_x_init):
    #total_z=1
    #for k in range(len(null_space)):
    #    if null_space[k]:
    #        total_z*=pairs[k][-1]
    #print("total_z: "+str(total_z))

    g = [1, f_x_init[1]]
    for i in range(2, len(f_x_init)): 
       g.append(f_x_init[i]*pow(f_x_init[0], i-1))
    inert_set = []
    if m1 != f_x_init[0]:
        print("oops")
        sys.exit()

    for p in primes:
    #    print(" irreducible: "+str(irreducibility(g, p))+" g: "+str(g)+" p: "+str(p))
        if m1%p and f_x_init[0]%p and irreducibility(g, p): 
            
            inert_set.append(p)
    if len(inert_set)==0:
        print("fatal error empty inert: "+str(g)+" f_x_init: "+str(f_x_init))
        sys.exit()
        return 0,0

    f_x=g
    inert=inert_set[-1]
    #print("g: "+str(f_x)+" f_x: "+str(f_x_init))
    d=2
    f_prime = get_derivative(f_x_init)
    g_prime = get_derivative(f_x)
    f_prime_sq,f_prime_eval = div_poly(poly_prod(g_prime, g_prime), f_x), pow(f_x[0], d-2, n)*eval_F(m0, m1, f_prime, d-1)%n


    S = 0
    rat_square=1
    alg_square=1

    rational_square = [i for i in f_prime_sq]
   # print("rational_square: "+str(rational_square))
    count=0
    for k in range(len(null_space)):
        if null_space[k]:
            rat_square*=pairs[k][3]
            alg_square*=pairs[k][1]
         #   print("using pair: "+str(pairs[k]))
            count+=1
            if pairs[k][1]==0:
                continue
            rational_square=poly_prod(rational_square, pairs[k][0])
    if count < 3:
        return 0,0
    rat_test=math.isqrt(rat_square)
    alg_test=math.isqrt(alg_square)
    if rat_test**2 != rat_square or alg_test**2 != alg_square:
        print("something went wrong lol")
        sys.exit()
  #  print("rational_square: "+str(rational_square))

    for k in range(len(null_space)):
        if null_space[k]:    
            if pairs[k][1]==0:
                continue
            rational_square = div_poly(rational_square, f_x)
            S += pairs[k][6] 
    
  #  print("rational_square: "+str(rational_square))

    x = f_prime_eval
    rat=create_rational(null_space, n, len_primes, primes, pairs)
    x = x*rat%n
    lead=pow(leading, S>>1, n)  
    x = x*lead%n


    f_norm = 0
    tmp = 1
    for x1 in f_x:
        f_norm += x1*x1*tmp
        tmp *= leading
    f_norm = int(math.sqrt(f_norm))+1
    fd = int(pow(len(f_x)-1, 1.5))+1


    coeff_bound = [fd*pow(f_norm, len(f_x)-1-i)*pow(2*(leading*u)*f_norm, S>>1) for i in range(len(f_x)-1)]
    y = square_root(f_x, rational_square, inert, m0, m1, leading, max(coeff_bound))
    if y == -1:
        return 0,0
    y = y*pow(m1, S>>1, n)%n
  
    
    return x, y
    
def create_rational(null_space, n, len_primes, primes, pairs):
    x = 1
    vec = [0]*(len_primes)

    for z in range(len(null_space)):
        if null_space[z]:
            if pairs[z][1]==0:
                continue
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
    


def compute_factors(pairs_used, vec, n, primes, m0, m1, leading_coeff, M, time_1,Q,f_x):

    x, y = create_solution(pairs_used,vec,n,len(primes),primes,m0,m1,leading_coeff,M<<1,Q,f_x)
    print("x**2%n: "+str(x**2%n)+" y**2%n: "+str(y**2%n))
    if x**2%n != y**2%n:
        print("FATAL")
        sys.exit()
    if x != y and math.gcd(x-y, n) != 1 and math.gcd(x+y, n) != 1:
        print_final_message(x, y, n, time_1)
        return True, 1

    return False, 0
    
def print_final_message(x, y, n, time_1):
    time_2 = datetime.now()

    print("found factor : "+str(math.gcd(x-y,n)))
    print("found factor : "+str(math.gcd(x+y,n)))

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
            if r[3]==0:
                break
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
            if r[1]==0:
                break
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
            if r[1]==0:
                break
            if compute_legendre_character(eval_mod(r[5], qua[p][1], qua[p][0]), qua[p][0]) == -1: 
                line.append(index)
            index += 1
       # print("div_lead: "+str(div_lead))
        for p in range(len(div_lead)):
            if r[1]==0:
                break
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
    

    
def power(poly, f, p, exp):
    if exp == 1: return poly

    tmp = power(poly, f, p, exp>>1)
    tmp = poly_prod(tmp, tmp)

    if exp&1:
        tmp = poly_prod(tmp, poly)
        return div_poly_mod(tmp, f, p)
    
    else: return div_poly_mod(tmp, f, p)
    

    
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
def initialize_2(f_x, d, primes, leading_coeff):
    m1=f_x[0]
    pairs_used, R_p, logs, tmp = [], [], [], 10*d
    divide_leading, pow_div = [], []
    for p in primes:
        if not leading_coeff%p:
            divide_leading.append(p)
            u = p
            while not leading_coeff%u: u *= p
            pow_div.append(u//p)

        logs.append(round(math.log2(p)))
        if p > tmp: R_p.append(find_roots_poly(f_x, p))
        else: R_p.append(fast_roots(f_x, p))
        if len(R_p[-1]) == d:
            if d&1: pairs_used.append([[leading_coeff*p], leading_coeff*pow(p, d), [[1,1]], p*m1, [1], [0,p], 1])
            else: pairs_used.append([[leading_coeff*p], leading_coeff*p, [[1,1]], p*m1, [1], [0,p], 1])

    for p in divide_leading:
        for i in range(len(pairs_used)): pairs_used[i].append(True)
        
    B_prime = sum([len(i) for i in R_p]) + 1
    
    return pairs_used,R_p, logs, divide_leading, pow_div, B_prime
    
    
def initialize_3(n, f_x, f_prime, const1, leading_coeff):
   # print("hit")
    k = 3*n.bit_length()
    Q = []
    q = const1+1
    if not q%2: q += 1
    while len(Q) < k:
      #  print(q)
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
            roots=solve_quadratic_congruence(1, y, -n*k, prime)
            if roots == -1:
                y+=1
                continue
            y2=y
            while y2 < prime:
                for x in roots:
                    if (x**2+y2*x-n*k)%prime!=0: 
                        print("error22")
                        sys.exit()   
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

def build_matrix(factor_base, algebraic_factors,rational_factors,symbols):#,pflist):
    fb_map = {val: i for i, val in enumerate(factor_base)}
   # print("factor_base: "+str(factor_base)+" fb_map: "+str(fb_map))
    qbase=len(symbols[0])
    algbase=len(factor_base)
    ratbase=len(factor_base)
    ind=1
   # print("total length: "+str(algbase+ratbase+qbase))
    M2=[0]*(algbase+ratbase+qbase)
    for i in range(len(algebraic_factors)):
        for fac in algebraic_factors[i]:
            idx = fb_map[fac]
      #      print("idx: "+str(idx))
            M2[idx] |= ind
        ind = ind + ind       


    offset=algbase
    ind=1
    for i in range(len(rational_factors)):
        for fac in rational_factors[i]:
            idx = fb_map[fac]
            M2[idx+offset] |= ind
        ind = ind + ind

    offset=algbase+ratbase
    ind=1
    
    for i in range(len(symbols)):
        j=0
        while j < len(symbols[i]):
            fac=symbols[i][j]
            if fac == 1:
                M2[j+offset] |= ind
            j+=1
        ind = ind + ind
    return M2

def solve_bits(matrix,factor_base,length):
    n=length#len(factor_base)*1#base+2
    lsmap = {lsb: 1 << lsb for lsb in range(n+10000)}
    m = len(matrix)
    marks = []
    cur = -1
    mark_mask = 0
    for row in matrix:
       # if cur % 100 == 0:
         #   print("", end=f"{cur, m}\r")
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

def check_poly_irreducible(poly_coeffs):
    """
    Checks if a polynomial defined by coefficients [an..., a0]
    is irreducible over the Rational field (Q).
    """
    x = symbols('x')
    # Construct polynomial: [4, 7, 21, 28] -> 4x^3 + 7x^2 + 21x + 28
    p = Poly(poly_coeffs, x)
    return p.is_irreducible

def try_square_root(ret_pairs,f_x,n,M,primelist_solve,primelist_f,B,primelist):
    pairs_used=ret_pairs
    rational_factors=[]
    algebraic_factors=[]

    symbols=[]

    z=1
    d=2
    free_pairs,R_p, logs, divide_leading, pow_div, B_prime = initialize_2(f_x, d, primelist, z)
    pairs_used.extend(free_pairs)
    f_prime = get_derivative(f_x)
    Q, k1 = initialize_3(n, f_x, f_prime, B, z)

    k=0
    while k < len(pairs_used):
        local_factors2, value2,seen_primes2,seen_primes_indexes2 = factorise_fast(pairs_used[k][1],primelist_f)
        if value2 !=1:
            print("fatal1")
            sys.exit()

        algebraic_factors.append(local_factors2)
        local_factors2, value2,seen_primes2,seen_primes_indexes2 = factorise_fast(pairs_used[k][3],primelist_f)
        if value2 !=1:
            print("fatal2: ",pairs_used[k])
            sys.exit()

        rational_factors.append(local_factors2)
        sym_line=[]
        for p in range(len(Q)):
            if compute_legendre_character(eval_mod(pairs_used[k][5], Q[p][1], Q[p][0]), Q[p][0]) == -1: 
                sym_line.append(1)
            else:
                sym_line.append(0)
        symbols.append(sym_line)
        k+=1

    M2 = build_matrix(primelist_solve, algebraic_factors,rational_factors,symbols)#,pflist)
    null_space=solve_bits(M2,primelist_solve,len(pairs_used))

    for vector in null_space:
        fvec=[]
        test_poly=1
        test_zxy=1
        for idx in range(len(pairs_used)):
            bit = vector & 1
            vector = vector >> 1
            if bit == 1:
              #  print("pairs_used: "+str(pairs_used[idx])+" co_list: "+str(co_list[idx]))
                test_poly*=pairs_used[idx][1]
                test_zxy*=pairs_used[idx][3]
                fvec.append(1)
            else:
                fvec.append(0)
              #  print("polyval:  "+str(relations[idx])+" disc constant "+str(x_list[idx]))
            idx += 1

        test_poly_sqr=math.isqrt(test_poly)
        test_zxy_sqr=math.isqrt(test_zxy)
        if test_poly_sqr**2 != test_poly:
            print("something went wrong test_poly")
            sys.exit()
        if test_zxy_sqr**2 != test_zxy:
            print("something went wrong test_zxy")
            sys.exit()
        x, y = create_solution(pairs_used,fvec,n,len(primelist),primelist,-f_x[1],f_x[0],f_x[0],M<<1,Q,f_x)
        if x ==0 and y==0:
            continue
        gcdres=gcd(x+y,n)
        print("gcd: "+str(gcdres))
        if gcdres != 1 and gcdres != n:
            sys.exit()

def find_common_coefficient(pairs_used,pairs_flist,pair_unfactored,n,M,primelist_solve,primelist_f,B,primelist):
    ##TO DO: We need some kind of formula to quickly find a common coefficient. Everything I know so far indicates that it should be possible. but it is a little complicated.
    ##For now this is just doing bruteforcing.


    k=1
    while k < 100:
        local_factors, value,seen_primes,seen_primes_indexes = factorise_fast(k,primelist_f)
        if value !=1:
            k+=1
            continue        
        hmp_pairs={}
        i=0
        while i < len(pairs_used):
            pair=pairs_used[i]
            poly_val=pair[1]*k
            a=pair[3]*k
            value2=pair_unfactored[i]
            a=a//value2
            x_start=value2
            y_start=a-x_start
            z_start=1
            seen_primes2=copy.deepcopy(pairs_flist[i])
            seen_primes2.extend(seen_primes)
            seen_primes2.append(1)
            seen_primes2.sort()
            u=0
            prev=1
            exp=1
            while u < len(seen_primes2):
                z=1
                factor=seen_primes2[u]
                if factor == prev:

                    exp+=1
                    factor=factor**exp
                else:
                    prev=factor
                    exp=1
                x=x_start*factor 
                if (a*value2) % x !=0:
                    print("fatal: "+str(x)+" a: "+str(a)+" value2: "+str(value2)+" k: "+str(k)+" x_start: "+str(x_start)+" factor: "+str(factor)+" pair[3]: "+str(pair[3])+" seen_primes: "+str(seen_primes)+" seen_primes2: "+str(seen_primes2))
                    sys.exit()
                y=((a*value2)//x)-x

                
                if z*x**2+(y)*x-n*k != poly_val:
                    print("something went wrong")
                    sys.exit()

                zxy=z*x+(y)
             #   pair=[[-1, x*z], poly_val, [[1 ,1]], zxy, [1], [-1,x], 1]
                if poly_val!=0 and zxy !=0:
                    #print("zxy: "+str(zxy)+" y: "+str(y)+" x: "+str(x)+" z: "+str(z)+" pair: "+str(pair))
                    try:
                        test=hmp_pairs[y]
                        test.append([i,x,y,z,zxy,poly_val])

                    except Exception as e:
                   #     print(e)

                        hmp_pairs[y]=[[i,x,y,z,zxy,poly_val]]
                u+=1

            #pair=pairs_used[i]
            #print("pairs_used: "+str(pair))#+" co: "+str(co_list[i]))
            #temp=(pair[1]+n*2)#*2#len(pairs_used)
            #x=zxy_total
            #y=temp-x#temp
            #poly=x**2+y*x-n*x*4#len(pairs_used)
            #print("poly: "+str(poly)+" x: "+str(x)+" y: "+str(y))
            i+=1

       
        for y, test in hmp_pairs.items():

            if len(test)>(base*2):
                ret_pairs=[]
                for ind in test:
                                

                    new_pair=[[-1, ind[1]*ind[3]], ind[5], [[1 ,1]], ind[4], [1], [-1,ind[1]], 1]
                    ret_pairs.append(new_pair)
                #print("found one: ",y)
                f_x=[1,y,-n*k]
                if check_poly_irreducible(f_x) != 1:
                    continue

                try_square_root(ret_pairs,f_x,n,M,primelist_solve,primelist_f,B,primelist)
        k+=1

    return 0,0

cdef construct_interval(list ret_array,partials,n,primeslist,hmap,hmap2,large_prime_bound,primeslist2,small_primeslist):
   # i=0
   # while i < len(hmap):
     #   k=0
      #  while k < len(hmap[i]):
          #  print("prime: "+str(primeslist[i])+" k: "+str(k+1)+" "+str(hmap[i][k]))
      #      k+=1
     #   i+=1
    
    sievelen=lin_sieve_size
    primelist_f=copy.copy(primeslist)
    primelist_f.insert(0,len(primelist_f)+1)
    primelist_f=array.array('q',primelist_f)
    poly_smooths=[]
    primelist=copy.copy(primeslist)
    primelist.insert(0,2) ##To do: remove when we fix lifting for powers of 2

    pair_list=[]
    co_list=[]
    primelist_solve=copy.copy(primeslist)
    primelist_solve.insert(0,2) ##To do: remove when we fix lifting for powers of 2
    primelist_solve.insert(0,-1)

    pairs_flist=[]
    pair_unfactored=[]
    
    stop=0
    print("[i]Phase one, sieving a-Nk = p")
    mult=1
    while mult < 100_000:
        k_ind=0
        while k_ind < 1:
            a_ind=0
            while a_ind < sievelen: 
                k=1+k_ind
                k_sq=k**2
                a=1+a_ind
                a*=mult
                poly_val=a-n*k_sq
                local_factors, value,seen_primes,seen_primes_indexes = factorise_fast(poly_val,primelist_f)
                if value !=1:
                    a_ind+=1
                    continue        

            
                ##NOTE: This does not have to factorize... 
                local_factors2, value2,seen_primes2,seen_primes_indexes2 = factorise_fast(a,primelist_f)
                if value2 ==a:
                    a_ind+=1
                    continue
                
                pair=[[-1, 1], poly_val, [[1 ,1]], a, [1], [-1,1], 1]
                if pair not in pair_list:
                    pairs_flist.append(seen_primes2)
                    pair_list.append(pair)
                    pair_unfactored.append(value2)
                    if len(pair_list)>30_000:
                        stop=1

                if stop ==1:
                    break

                a_ind+=1  
            if stop ==1:
                break  
            k_ind+=1
        if stop ==1:
            break
        mult+=1


    print("[i]Phase two, checking results")
    orig_pairlist=copy.deepcopy(pair_list)
    d=2
    B=primeslist[-1]

    pairs_used=pair_list
    M=sievelen
    ret_pairs,f_x=find_common_coefficient(pairs_used,pairs_flist,pair_unfactored,n,M,primelist_solve,primelist_f,B,primelist)
    if ret_pairs == 0:
        print("no solutions found")
        sys.exit()
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
    l=0
    while l < 1:
    
        start = default_timer() 

        if g_p !=0 and g_q !=0 and g_enable_custom_factors == 1:
            p=g_p
            q=g_q
            key=p*q
        if key == 0:
            print("\n[*]Generating rsa key with a modulus of +/- size "+str(keysize)+" bits")
            publicKey, privateKey,phi,p,q = generateKey(keysize//2)
            n=p*q
            #key=n
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
        l+=1