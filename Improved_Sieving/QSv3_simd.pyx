##NFS related code is borrowed from: https://github.com/basilegithub/General-number-field-sieve-Python (note: Very impressively written, helped me big time, thanks)


from sympy import symbols, Poly, Matrix
import argparse
import cProfile
from numpy.polynomial.legendre import leggauss
import math, os
from datetime import datetime
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
import gc
import array
import numpy as np
from datetime import datetime
import sys
from datetime import datetime
import sys
import math
from scipy.linalg import det
import math, sys
from datetime import datetime
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

# batch smooth test using a product tree
# more efficient than testing for smoothness every candidate : should be prefered
# Note that in this function, we test for smoothness abs(candidates), see QS.sieve_and_batch_smooth function
# Source: https://cr.yp.to/factorization/smoothparts-20040510.pdf (Bernstein)
    
def batch_smooth_test(candidates, prod_primes, cst_1, cst_2, div):
    L = len(candidates)
    #print("candidates: "+str(candidates))
    # compute the smallest e such that 2**(2**e) >= max(candidates), ie e = int(log2(log2(max(candidates)))
    e = 0
    tmp = 2
    tmp_max = max([abs(candidate[1])//div for candidate in candidates] + [abs(candidate[3]) for candidate in candidates])
    while tmp < tmp_max:
        tmp = tmp*tmp
        e += 1
        
    # build the product tree
    tmp = [abs(candidate[1])//div for candidate in candidates] + [abs(candidate[3]) for candidate in candidates]
    prod_tree = [tmp]
    
    while len(prod_tree[-1]) > 1:
        line = [0]*(len(prod_tree[-1])>>1)

        for i in range(1, len(prod_tree[-1]), 2):
            tmp = prod_tree[-1][i]*prod_tree[-1][i-1]
            if tmp <= prod_primes: line[i>>1] = tmp
            else: line[i>>1] = prod_primes+1

        if len(prod_tree[-1]) & 1: line.append(prod_tree[-1][-1])
        prod_tree.append(line)
    
    # update the product tree by computing the adequate remainders

    prod_tree[-1][0] = prod_primes%prod_tree[-1][0]
    for i in range(len(prod_tree)-1):
        for j in range(len(prod_tree[-i-1])-1):
            prod_tree[-i-2][(j<<1)] = prod_tree[-i-1][j]%prod_tree[-i-2][j<<1]
            prod_tree[-i-2][(j<<1)+1] = prod_tree[-i-1][j]%prod_tree[-i-2][(j<<1)+1]

        tmp = len(prod_tree[-i-1])-1
        prod_tree[-i-2][tmp<<1] = prod_tree[-i-1][tmp]%prod_tree[-i-2][tmp<<1]
        if (tmp<<1)+1 < len(prod_tree[-i-2]):
            prod_tree[-i-2][(tmp<<1)+1] = prod_tree[-i-1][tmp]%prod_tree[-i-2][(tmp<<1)+1]

    # test for smoothness
    smooth = [0]*L
    for i in range(L):
        tmp = 0
        j_1, j_2 = prod_tree[0][i], prod_tree[0][L+i]
        while (j_1 or j_2) and tmp < e:
            j_1 = (j_1*j_1)%(abs(candidates[i][1]//div))
            j_2 = (j_2*j_2)%abs(candidates[i][3])
            tmp += 1

        if not j_1 and not j_2: smooth[i] = [True, 1, 1, 1, 1]

        else:
            tmp_1 = (abs(candidates[i][1])//div)//math.gcd(abs(candidates[i][1]//div), j_1)
            tmp_2 = abs(candidates[i][3])//math.gcd(abs(candidates[i][3]), j_2)
            
            if tmp_1 > 1 and tmp_2 > 1: smooth[i] = [False, 1, 1, 1, 1]

            elif tmp_1 > 1:
                if tmp_1 < cst_1:
                    smooth[i] = [True, 1, tmp_1, 1, 1]
                elif tmp_1 < cst_2 and fermat_primality(tmp_1):
                    tmp_1 = pollard_rho(tmp_1)
                    smooth[i] = [True, min(tmp_1), max(tmp_1), 1, 1]
                else:
                    smooth[i] = [False, 1, 1, 1, 1]

            elif tmp_2 > 1:
                if tmp_2 < cst_1:
                    smooth[i] = [True, 1, 1, 1, tmp_2]
                elif tmp_2 < cst_2 and fermat_primality(tmp_2):
                    tmp_2 = pollard_rho(tmp_2)
                    smooth[i] = [True, 1, 1, min(tmp_2), max(tmp_2)]
                else:
                    smooth[i] = [False, 1, 1, 1, 1]



            else: smooth[i] = [False, 1, 1, 1, 1]
            
    return smooth

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

def dickman(x, table):
    k,res = math.ceil(x),0
    delta = k-x
    if k-1 > len(table): table = get_dickman_table(k)
    tmp = 1
    for i in range(len(table[k-2])):
        res += table[k-2][i]*tmp
        tmp *= delta

    return res, table

def get_dickman_table(k):
    coeffs = [[1-math.log(2)]+[1/(i*(1<<i)) for i in range(1, 30)]]
    for i in range(3, k+1):
        new = [0]*30
        for u in range(1, 30):
            c = 0
            for j in range(u): c += coeffs[-1][j]/(u*pow(i, u-j))
            new[u] = c

        c = 0
        for j in range(1, len(new)): c += new[j-1]/(j+1)
        new[0] = c/(i-1)
        coeffs.append(new)

    return coeffs
    
def my_norm(z):
    if z.real == 0: return math.ceil(abs(z.imag))
    if z.imag == 0: return math.ceil(abs(z.real))

    return math.isqrt(math.ceil(abs(z.imag))**2 + math.ceil(abs(z.real))**2)+1

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

def sieve(length, f_x, rational_base, algebraic_base, m0, m1, b, logs, offset, leading,div):
    offset += div.bit_length()-1
    pairs = []
    tmp_poly = new_coeffs(f_x, b)
  #  print("b: ",b)
    sieve_array = [0]*(length<<1)
    for q, p in enumerate(rational_base):
        tmp_len = length%p
        log = logs[q]

        if m1%p:
            root = (tmp_len+b*m0*invmod(m1, p))%p
            for i in range(root, len(sieve_array), p): sieve_array[i] += log

        for r in algebraic_base[q]:
            root = (tmp_len+b*r)%p
            for i in range(root, len(sieve_array), p): sieve_array[i] += log

    if b&1:
     #   print("hit")
        eval2 = -length*m1-b*m0
        a = -length
        tmp = [0]*len(tmp_poly)
        for j in range(len(tmp_poly)): tmp[j] = evaluate(tmp_poly, -length+j)
        for q in range(1, len(tmp_poly)):
            for k in range(len(tmp_poly)-1, q-1, -1): tmp[k] -= tmp[k-1]

        eval1 = tmp[0]
        for k in range(len(sieve_array)):
            if a and math.gcd(a, b) == 1 and eval2:
                eval = abs(eval1*eval2)
               # print("eval: ",eval)
                if eval != 0 and sieve_array[k] > eval.bit_length()-offset:
                    pairs.append([[-b, a*leading], eval1, [[1 ,1]], eval2, [1], [-b,a], 1])
            a += 1
            eval2 += m1
            eval1 += tmp[1]
            for q in range(1, len(tmp_poly)-1): tmp[q] += tmp[q+1]
    else:
     #   print("hit2")
        init = 0
        if not length&1:
            length -= 1
            init = 1
        eval2 = -length*m1-b*m0
        a = -length
        tmp = [0]*len(tmp_poly)
        for j in range(len(tmp_poly)): tmp[j] = evaluate(tmp_poly, -length+2*j)
        for q in range(1, len(tmp_poly)):
            for k in range(len(tmp_poly)-1, q-1, -1): tmp[k] -= tmp[k-1]

        eval1 = tmp[0]
        for k in range(init, len(sieve_array), 2):
            if math.gcd(a, b) == 1 and eval2:
                eval = abs(eval1*eval2)
                if eval != 0 and sieve_array[k] > eval.bit_length()-offset:
                    pairs.append([[-b, a*leading], eval1, [[1, 1]], eval2, [1], [-b,a], 1])
            a += 2
            eval2 += m1<<1
            eval1 += tmp[1]
            for q in range(1, len(tmp_poly)-1): tmp[q] += tmp[q+1]
 #   print("pairs: ",pairs)
    return pairs
def square_root(f, N, p, m0, m1, leading, bound):
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
    
def create_solution(pairs, null_space, n, len_primes, primes, f_x, m0, m1, inert, f_prime_sq, leading, f_prime_eval, u):
    f_norm = 0
    tmp = 1
    for x in f_x:
        f_norm += x*x*tmp
        tmp *= leading
    f_norm = int(math.sqrt(f_norm))+1
    fd = int(pow(len(f_x)-1, 1.5))+1
    S = 0
    
    x = f_prime_eval
    x = x*create_rational(null_space, n, len_primes, primes, pairs)%n
    
    rational_square = [i for i in f_prime_sq]
    for k in range(len(null_space)):
        if null_space[k]:
            rational_square = div_poly(poly_prod(rational_square, pairs[k][0]), f_x)
            S += pairs[k][6]
            
    x = x*pow(leading, S>>1, n)%n
            
    coeff_bound = [fd*pow(f_norm, len(f_x)-1-i)*pow(2*(leading*u)*f_norm, S>>1) for i in range(len(f_x)-1)]
    
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

    for i in range(len_primes): x = x*pow(primes[i], vec[i]>>1, n)%n

    return x
    


def compute_factors(pairs_used, vec, n, primes, g, g_prime, g_prime_sq, g_prime_eval, m0, m1, leading_coeff,
                    inert_set, M, time_1):

    x, y = create_solution(pairs_used, vec, n, len(primes), primes, g, m0, m1, inert_set[-1], g_prime_sq,
                              leading_coeff, g_prime_eval, M<<1)

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
    for r in pairs:
        index = 1
        
        if r[3] < 0: line = [0]
        else: line = []

        for p in rat:
            tmp = p
            tmp2 = 0
            while not r[3]%tmp:
                tmp *= p
                tmp2 ^= 1
            if tmp2: line.append(index)
            index += 1
            
        if r[1] < 0: line.append(index)
        index += 1

        for p in range(len(alg)):
            tmp = rat[p]
            tmp2 = 0
            while not r[1]%tmp:
                tmp *= rat[p]
                tmp2 ^= 1

            for i in range(len(alg[p])):
                if not eval_mod(r[5], alg[p][i], rat[p]):
                    if tmp2: line.append(index)
                index += 1
                
        for p in range(len(qua)):
            if compute_legendre_character(eval_mod(r[5], qua[p][1], qua[p][0]), qua[p][0]) == -1: line.append(index)
            index += 1
            
        for p in range(len(div_lead)):
            if r[7+p]:
                tmp = div_lead[p]
                tmp2 = 0

                while not r[1]%tmp:
                    tmp *= div_lead[p]
                    tmp2 ^= 1
                if tmp2: line.append(index)

            index += 1
            
        if r[6]&1: line.append(index)
        matrix.append(set(line))

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
                        
# Build the dense binary matrix by computing explicitely every A[i][j]
def build_dense_matrix(pairs, rat, alg, qua, div_lead):
    matrix = []
    for r in pairs:
        line = []
        
        if r[3] < 0: line.append(1)
        else: line.append(0)

        for p in rat:
            tmp = p
            tmp2 = 0
            while not r[3]%tmp:
                tmp *= p
                tmp2 ^= 1
            line.append(tmp2)
            
        if r[1] < 0: line.append(1)
        else: line.append(0)

        for p in range(len(alg)):
            tmp = rat[p]
            tmp2 = 0
            while not r[1]%tmp:
                tmp *= rat[p]
                tmp2 ^= 1
            for i in range(len(alg[p])):
                if not eval_mod(r[5], alg[p][i], rat[p]):
                    line.append(tmp2)
                else:
                    line.append(0)
                
        for p in range(len(qua)):
            if compute_legendre_character(eval_mod(r[5], qua[p][1], qua[p][0]), qua[p][0]) == -1:
                line.append(1)
            else:
                line.append(0)
            
        for p in range(len(div_lead)):
            if r[7+p]:
                tmp = div_lead[p]
                tmp2 = 0
                while not r[1]%tmp:
                    tmp *= div_lead[p]
                    tmp2 ^= 1
                line.append(tmp2)
            else: line.append(0)
            
        line.append(r[6]&1)
        matrix.append(line)
        
    return matrix
    
def siqs_build_matrix_opt(M):
    """Convert the given matrix M of 0s and 1s into a list of numbers m
    that correspond to the columns of the matrix.
    The j-th number encodes the j-th column of matrix M in binary:
    The i-th bit of m[i] is equal to M[i][j].
    """
    m = len(M[0])
    cols_binary = [""] * m
    for mi in M:
        for j, mij in enumerate(mi):
            cols_binary[j] += "1" if mij else "0"
    return [int(cols_bin[::-1], 2) for cols_bin in cols_binary], len(M), m
def find_cycle_fp(graph, init):
    if init[0] == 1: return DFS_fp(graph, init)

    path_p1_to_1 = DFS_fp(graph, [1, init[0]])

    if path_p1_to_1 is not None:
        path_p2_to_1 = DFS_fp(graph, [1, init[1]])

        if path_p2_to_1 is not None:
            index1, index2 = len(path_p1_to_1)-1, len(path_p2_to_1)-1
            while path_p1_to_1[index1] == path_p2_to_1[index2]:
                index1 -= 1
                index2 -= 1
            return path_p2_to_1[:index2+2] + path_p1_to_1[:index1+1][::-1]

    return DFS_fp(graph, init)

# Finds cycles, aka relations, in the graph of relations with one or two large primes in the rational element
def DFS_fp(graph, init):
    path = [init[1]]
    stack = []
    neighbors = [init[1]]
    
    # Initialize the potential paths
    tmp = graph[init[1]]

    if init[0] in tmp:
        return path + [init[0]]
    
    neighbors += tmp
    
    stack.append(neighbors)
    
    # While there is a potential path that has not been explored
    while len(stack):
        while len(stack) and len(stack[-1]) == 1:
            path.pop()
            stack.pop()
        if not len(stack): break

        next_node = stack[-1][-1]
        stack[-1].pop()
        path.append(next_node)
        neighbors = [next_node]

        if init[0] in graph[next_node]:
            return path + [init[0]]
        
        parent = stack[-1][0]

        for i in range(len(graph[next_node])):
            if graph[next_node][i] != parent: neighbors.append(graph[next_node][i])
            else:
                neighbors += graph[next_node][i+1:]
                break

        stack.append(neighbors)
        
# Given a cycle in the fp graph, constructs the full relation
# The graph is always acyclic: when we find a new cycle, the partial relation we just found is not added
# Thus, one new partial relation can only create one cycle at most
def combine_fp(pair, pairs_used, fp, path, divide_leading, g):
    tmp = [i for i in pair]
    for i in range(len(path)-1):
        firstp, secondp = min(path[i],path[i+1]), max(path[i],path[i+1])

        combine = fp[firstp][secondp]
        tmp2 = [div_poly(poly_prod(tmp[0], combine[0]), g),
                tmp[1]*combine[1],
                [[1, 1]],
                tmp[3]*combine[3],
                tmp[4] + [path[i]],
                poly_prod(tmp[5],combine[5]),
                tmp[6]+1]
        for p in range(len(divide_leading)):
            tmp2.append(tmp[7+p] or combine[7+p])

        tmp = [i for i in tmp2]
    pairs_used.append(tmp)
    
# Master function, keeps the partial_relations and possible_smooth lists sorted, find the matching large primes, create the full relations from large primes
def handle_large_fp(pair, large_primes, pairs_used, fp, graph_fp, size_fp, parent_fp, g, divide_leading, cycle_len, full_found, partial_found_fp):
    pair[4] = [large_primes[3]]

    if large_primes[3] == large_primes[4]:
        pairs_used.append(pair)
        full_found += 1

    elif not bool(fp):
        pair[4].append(large_primes[4])

        small_p, big_p = large_primes[3], large_primes[4]

        fp[small_p] = {}
        fp[small_p][big_p] = pair

        graph_fp[small_p] = [big_p]
        graph_fp[big_p] = [small_p]

        parent_fp[small_p] = small_p
        parent_fp[big_p] = small_p

        size_fp += 1

    else:
        small_p, big_p = large_primes[3], large_primes[4]

        flag_small_prime = small_p in graph_fp
        
        # If the smallest prime has not already been seen
        if not flag_small_prime:
            graph_fp[small_p] = [big_p]

            # Find where in partial relations it is needed to insert the new partial relation
            pair[4].append(big_p)

            fp[small_p] = {}
            fp[small_p][big_p] = pair

            size_fp += 1
            
            # Find if largest prime has already been seen
            flag_big_prime = big_p in graph_fp

            if flag_big_prime: # if seen add the smallest prime to the list of links of the largest prime
                graph_fp[big_p].append(small_p)

                parent_fp[small_p] = parent_fp[big_p]

            else: # If it has not been seen, append it
                graph_fp[big_p] = [small_p]

                parent_fp[small_p] = small_p
                parent_fp[big_p] = small_p
            
        # If the smallest prime has already been seen
        else:
            flag_big_prime = big_p in graph_fp
            
            # If big prime has not been seen
            if not flag_big_prime:
                graph_fp[small_p].append(big_p) # Append the largest prime to the list of links of the smallest prime
                graph_fp[big_p] = [small_p] # Create the large prime in the graph

                parent_fp[big_p] = parent_fp[small_p]

                pair[4].append(big_p)
                # Find where to insert the partial relation
                if small_p not in fp:
                    fp[small_p] = {}

                fp[small_p][big_p] = pair

                size_fp += 1

            # If the largest prime has been seen, ie if both primes have already been seen
            else:

                parent_small_p = parent_fp[small_p]
                while parent_fp[parent_small_p] != parent_small_p:
                    parent_fp[parent_small_p], parent_small_p = parent_fp[parent_fp[parent_small_p]], parent_fp[parent_small_p]

                parent_big_p = parent_fp[big_p]
                while parent_fp[parent_big_p] != parent_big_p:
                    parent_fp[parent_big_p], parent_big_p = parent_fp[parent_fp[parent_big_p]], parent_fp[parent_big_p]

                if parent_small_p != parent_big_p:
                    graph_fp[small_p].append(big_p)
                    graph_fp[big_p].append(small_p)

                    parent_fp[parent_big_p] = parent_small_p

                    pair[4].append(big_p)
                    # Find where to insert the partial relation
                    if small_p not in fp:
                        fp[small_p] = {}
                        
                    fp[small_p][big_p] = pair

                    size_fp += 1
                else:
                    path_cycle = find_cycle_fp(graph_fp, [small_p, big_p])
                    if len(path_cycle) < 11: cycle_len[len(path_cycle)-2] += 1
                    else: cycle_len[-1] += 1
                    combine_fp(pair, pairs_used, fp, path_cycle, divide_leading, g)
                    partial_found_fp += 1

    return pairs_used, fp, graph_fp, size_fp, parent_fp, cycle_len, full_found, partial_found_fp

def find_cycle_pf(graph, init):
    if init[0] == (1, 1): return DFS_pf(graph, init)

    path_p1_to_1 = DFS_pf(graph, [(1, 1), init[0]])

    if path_p1_to_1 is not None:
        path_p2_to_1 = DFS_pf(graph, [(1, 1), init[1]])

        if path_p2_to_1 is not None:
            index1, index2 = len(path_p1_to_1)-1, len(path_p2_to_1)-1
            while path_p1_to_1[index1] == path_p2_to_1[index2]:
                index1 -= 1
                index2 -= 1
            return path_p2_to_1[:index2+2] + path_p1_to_1[:index1+1][::-1]

    return DFS_pf(graph, init)


# Finds cycles, aka relations, in the graph of relations with one or two large primes if the norm of the algebraic element
def DFS_pf(graph, init):
    path = [init[1]]
    stack = []
    neighbors = [init[1]]
    
    tmp = graph[init[1]]

    if init[0][0] in tmp and init[0] in tmp[init[0][0]]:
        return path + [init[0]]

    # Initialize the potential paths
    for key in tmp.keys():
        neighbors += tmp[key]
    
    stack.append(neighbors)
    
    # While there is a potential path that has not been explored
    while len(stack):
        while len(stack) and len(stack[-1]) == 1:
            path.pop()
            stack.pop()
        if not len(stack): break
        next_node = stack[-1][-1]
        stack[-1].pop()
        path.append(next_node)
        neighbors = [next_node]

        tmp = graph[next_node]

        if init[0][0] in tmp and init[0] in tmp[init[0][0]]:
            return path + [init[0]]
        
        parent = stack[-1][0]

        for key in tmp.keys():
            if key == parent[0]:
                for i in range(len(tmp[key])):
                    if tmp[key][i] != parent: neighbors.append(tmp[key][i])
                    else:
                        neighbors += tmp[key][i+1:]
                        break

            else:
                neighbors += tmp[key]

        stack.append(neighbors)


# Given a cycle in the pf graph, constructs the full relation
# The graph is always acyclic: when we find a new cycle, the partial relation we just found is not added
# Thus, one new partial relation can only create one cycle at most
def combine_pf(pair, pairs_used, pf, path, divide_leading, g):
    tmp = [i for i in pair]
    for i in range(len(path)-1):
        if path[i][0] <= path[i+1][0]: firstp, secondp = path[i], path[i+1]
        else: firstp, secondp = path[i+1], path[i]

        combine = pf[firstp][secondp]
        tmp2 = [div_poly(poly_prod(tmp[0], combine[0]), g),
                tmp[1]*combine[1],
                tmp[2] + [path[i]],
                tmp[3]*combine[3],
                [1],
                poly_prod(tmp[5], combine[5]),
                tmp[6]+1]
        for p in range(len(divide_leading)):
            tmp2.append(tmp[7+p] or combine[7+p])

        tmp = [i for i in tmp2]

    pairs_used.append(tmp)

def handle_large_pf(pair, large_primes, pairs_used, pf, graph_pf, size_pf, parent_pf, g, divide_leading, cycle_len, full_found, partial_found_pf):
    if large_primes[1] == 1:
        tmp = (1, 1) # [smallest_prime, r_smallest]
        pair[2] = [tmp]

    else:
        tmp = (large_primes[1], pair[5][1]*invmod(-pair[5][0], large_primes[1])%large_primes[1]) # [smallest_prime, r_smallest]
        pair[2] = [tmp]

    tmp2 = (large_primes[2], pair[5][1]*invmod(-pair[5][0], large_primes[2])%large_primes[2]) # [largest_prime, r_largest]

    if large_primes[1] == large_primes[2]:
        pairs_used.append(pair)
        full_found += 1

    elif not len(pf):
        pair[2].append(tmp2)
        pf[tmp] = {}
        pf[tmp][tmp2] = pair

        graph_pf[tmp] = {}
        graph_pf[tmp][tmp2[0]] = [tmp2]
        graph_pf[tmp2] = {}
        graph_pf[tmp2][tmp[0]] = [tmp]

        parent_pf[tmp] = tmp
        parent_pf[tmp2] = tmp

    else:
        flag_small_prime = tmp in graph_pf
        
        # If [smallest_prime, r_smallest] has not already been seen
        if not flag_small_prime:
            graph_pf[tmp] = {}
            graph_pf[tmp][tmp2[0]] = [tmp2] # Create smallest prime in graph

            pair[2].append(tmp2)
            pf[tmp] = {}
            pf[tmp][tmp2] = pair
            
            # Find if largest prime has already been seen
            flag_big_prime = tmp2 in graph_pf

            if flag_big_prime:
                if tmp[0] in graph_pf[tmp2]: # If seen, add the smallest prime to the list of links of the largest prime
                    graph_pf[tmp2][tmp[0]].append(tmp)
                else:
                    graph_pf[tmp2][tmp[0]] = [tmp]

                parent_pf[tmp] = tmp2

            else: # If it has not been seen, create it in the graph
                graph_pf[tmp2] = {}
                graph_pf[tmp2][tmp[0]] = [tmp]

                parent_pf[tmp] = tmp
                parent_pf[tmp2] = tmp

        # If [smallest_prime, r_smallest] has already be seen
        # We look for [largest_prime, r_largest]
        else:
            flag_big_prime = tmp2 in graph_pf

            if not flag_big_prime: # [largest_prime, r_biggest] not in graph
                if tmp2[0] in graph_pf[tmp]:
                    graph_pf[tmp][tmp2[0]].append(tmp2) # Create the link from small to big prime
                else:
                    graph_pf[tmp][tmp2[0]] = [tmp2] # Create the link from small to big prime

                graph_pf[tmp2] = {}
                graph_pf[tmp2][tmp[0]] = [tmp] # Create the large prime in the graph

                parent_pf[tmp2] = tmp

                pair[2].append(tmp2)
                # Find where to insert the partial relation
                if tmp not in pf:
                    pf[tmp] = {}

                pf[tmp][tmp2] = pair

                size_pf += 1

            else:
                parent_tmp = parent_pf[tmp]
                while parent_pf[parent_tmp] != parent_tmp:
                    parent_pf[parent_tmp], parent_tmp = parent_pf[parent_pf[parent_tmp]], parent_pf[parent_tmp]

                parent_tmp2 = parent_pf[tmp2]
                while parent_pf[parent_tmp2] != parent_tmp2:
                    parent_pf[parent_tmp2], parent_tmp2 = parent_pf[parent_pf[parent_tmp2]], parent_pf[parent_tmp2]

                if parent_tmp != parent_tmp2:
                    if tmp2[0] in graph_pf[tmp]:
                        graph_pf[tmp][tmp2[0]].append(tmp2)
                    else:
                        graph_pf[tmp][tmp2[0]] = [tmp2]

                    if tmp[0] in graph_pf[tmp2]:
                        graph_pf[tmp2][tmp[0]].append(tmp)
                    else:
                        graph_pf[tmp2][tmp[0]] = [tmp]

                    parent_pf[parent_tmp2] = parent_tmp

                    pair[2].append(tmp2)

                    if tmp not in pf:
                        pf[tmp] = {}
                        
                    pf[tmp][tmp2] = pair

                    size_pf += 1

                else:
                    path_cycle = find_cycle_pf(graph_pf, [tmp, tmp2])
                    if len(path_cycle) < 11: cycle_len[len(path_cycle)-2] += 1
                    else: cycle_len[-1] += 1
                    combine_pf(pair, pairs_used, pf, path_cycle, divide_leading, g)
                    partial_found_pf += 1

    return pairs_used, pf, graph_pf, size_pf, parent_pf, cycle_len, full_found, partial_found_pf

## Creation of polynomials
def find_relations(f_x, leading_coeff, g, primes, R_p, Q, B_prime, divide_leading, prod_primes, pow_div, pairs_used,
                   const1, const2, logs, m0, m1, M):
    fp, pf = {}, {}
    parent_fp, parent_pf = {}, {}
    full_found = 0
    partial_found_fp, partial_found_pf = 0, 0
    size_fp, size_pf = 0, 0
    graph_fp, graph_pf = {}, {}
    offset = 15+math.log2(const1)
    cycle_len = [0]*10
    
    V = 3 + len(primes) + B_prime + len(Q) + len(divide_leading)
    
    print("sieving...")
    print("need to find at least "+str(V+10)+" relations\n")
    
    time_1 = datetime.now()
    
    b = 1
        
   # print('\r'+"0/("+str(V)+"+10) relations found")
    while len(pairs_used) < V+10:
        
        div = 1
        for q in range(len(divide_leading)):
            p = divide_leading[q]
            if not b%p:
                tmp = p
                while not b%tmp and tmp <= pow_div[q]:
                    div *= p
                    tmp *= p
        
        pairs = sieve(M, f_x, primes, R_p, m0, m1, b, logs, offset, leading_coeff, div)


      #  smooth = batch_smooth_test(pairs, prod_primes, const1, const2, div)

        for i, pair in enumerate(pairs):
            z = trial(pair, primes, const1, const2, div)
         #   z = smooth[i]

            if z[0]:
                tmp = [u for u in pair]
                for p in divide_leading: 
                    tmp.append(not pair[5][0]%p)
                if z[2] == 1 and z[4] == 1:
                    pairs_used.append(tmp)
                    full_found += 1

                elif z[4] > 1 and z[2] == 1:
                    pairs_used, fp, graph_fp, size_fp, parent_fp, cycle_len, full_found, partial_found_fp = handle_large_fp(tmp, z, pairs_used, fp, graph_fp, size_fp, parent_fp, g, divide_leading,
                                                                                                                            cycle_len, full_found, partial_found_fp)

                elif z[4] == 1 and z[2] > 1:
                    pairs_used, pf, graph_pf, size_pf, parent_pf, cycle_len, full_found, partial_found_pf = handle_large_pf(tmp, z, pairs_used, pf, graph_pf, size_pf, parent_pf, g, divide_leading,
                                                                                                                            cycle_len, full_found, partial_found_pf)

        print('\r'+"b = "+str(b)+" "+str(len(pairs_used))+"/("+str(V)+"+10) ; full relations = "+str(full_found)+" | partial found fp = "+str(partial_found_fp)+" ("+str(size_fp)+") | partial found pf = "+str(partial_found_pf)+" ("+str(size_pf)+")")
        b += 1
        
    print("\n")
    
    time_2 = datetime.now()
    
    print("sieving done in "+format_duration(time_2-time_1)+".\n")
    print("Distribution of cycle length:")
    print("1-cycle: "+str(full_found))
  #  for i in range(9):
    #    print(str(i+2)+"-cycle: "+str(cycle_len[i]))
    print("11+-cycle: "+str(cycle_len[-1])+"\n")
    print(str(partial_found_fp)+" partial relations fp found")
    print(str(partial_found_pf)+" partial relations pf found")
                
    return pairs_used, V
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

## Operations on one polynomial
    
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
    
def get_complex_roots(f):
    rho = 0
    for i in range(1, len(f)): rho += abs(f[i])
    rho /= abs(f[0])
    rho = max(1,rho)
    d = len(f)-1
    roots = []
    for i in range(d): roots.append(rho*pow(math.cos(2*math.pi/d)+math.sin(2*math.pi/d)*1j, i))
    next_roots = [None]*d

    for _ in range(1000):
        for i in range(d):
            prod = f[0]
            for k in range(d):
                if k != i: prod *= (roots[i]-roots[k])
            next_roots[i] = roots[i]-evaluate(f, roots[i])/prod
        roots = [i for i in next_roots]
    return roots
    
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
    
def eval_G(x, y, m):
    return x-y*m    
## NFS functions

# This file contains all the functions required to do the polynomial selection step
def minimise_Lnom(f, s, B, m0, m1):
    k, res = 1, get_Lnorm(poly_prod(f,f),s,B)
    while k > 0:
        flag = False
        F = shift(f, -k)
        tmp = get_Lnorm(poly_prod(F, F), s, B)
        if tmp < res:
            res, f, k, flag, m0 = tmp, [i for i in F], k<<1, True, m0+k*m1

        F = shift(f, k)
        tmp = get_Lnorm(poly_prod(F, F), s, B)
        if tmp < res:
            res, f, k, flag, m0 = tmp, [i for i in F], k<<1, True, m0-k*m1
        if not flag: k >>= 1

    return f, m0
    
def evaluate_polynomial_quality(f_x, B, m0, m1, primes):
    _, s = get_sieve_region(f_x, B)

    if len(f_x) > 2:

        if (s*abs(f_x[-3])-abs(f_x[-2]))//m0 > 0:
            f_x = root_sieve2(f_x, [m1,-m0], primes[:150], (s*abs(f_x[-3])-abs(f_x[-2]))//m0, (s*s*abs(f_x[-3])-abs(f_x[-1]))//m0)
            _, s = get_sieve_region(f_x, B)
            f_x, m0 = minimise_Lnom(f_x, s, B, m0, m1)

        elif (s*abs(f_x[-2])-abs(f_x[-1]))//m0 > 0:
            f_x = root_sieve(f_x, [m1,-m0], primes[:150], (s*abs(f_x[-2])-abs(f_x[-1]))//m0)
            _, s = get_sieve_region(f_x, B)
            f_x, m0 = minimise_Lnom(f_x, s, B, m0, m1)

    M, s = get_sieve_region(f_x, B)
    table = get_dickman_table(20)

    alpha = alpha_score(f_x, primes[:350])
    E_score = get_Escore(f_x, [m1, -m0], alpha[0], B, M, round(B/math.sqrt(s)), table)
    E_score_2 = get_Epscore(f_x, [m1, -m0], alpha, B, M, round(B/math.sqrt(s)), table)
    L_norm = get_Lnorm(poly_prod(f_x, f_x), s, B)

    print("f(x) = "+str(f_x)+"    g(x) = "+str([m1, -m0]))
    print("alpha = "+str(alpha[0]))
    print("E-score = "+str(E_score))
    print("E-score 2 = "+str(E_score_2))
    print("L²-norm = "+str(L_norm))
    print("skew = "+str(s)+"\n")
    
    return f_x, m0, M
    
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

def get_Escore(f, g, alpha, B, x_limit, y_limit, table):
    n = 32
    log_B = math.log(B)
    
    x, w = leggauss(n)  # nodes & weights on [-1, 1]
    a, b = 0, math.pi
    xm = 0.5 * (b - a) * x + 0.5 * (b + a)
    wm = 0.5 * (b - a) * w

    res = 0
    for i in range(n):
        angle = xm[i]
        X = x_limit*math.cos(angle)
        Y = (y_limit-1)*math.sin(angle)+1
        res1, table = dickman((math.log(abs(eval_F(X, Y, f, len(f)-1)))-alpha)/log_B, table)
        res2, table = dickman((math.log(abs(eval_F(X, Y, g, len(g)-1))))/log_B, table)
        res += wm[i]*res1*res2

    return res

def get_Epscore(f, g, alpha, B, x_limit, y_limit, table):
    k,l,c = alpha[1],alpha[2],alpha[3]

    n = 16
    log_B = math.log(B)

    # Gauss–Legendre for mu in [c, c + 5c]
    mu_nodes, mu_weights = leggauss(n)
    a_mu, b_mu = c, 6*c
    mu_x = 0.5*(b_mu-a_mu)*mu_nodes + 0.5*(b_mu+a_mu)
    mu_w = 0.5*(b_mu-a_mu)*mu_weights

    # Gauss–Legendre for theta in [0, π]
    th_nodes, th_weights = leggauss(n)
    a_th, b_th = 0, math.pi
    th_x = 0.5*(b_th-a_th)*th_nodes + 0.5*(b_th+a_th)
    th_w = 0.5*(b_th-a_th)*th_weights

    res = 0
    for i in range(n):
        Y = non_central(k, l, mu_x[i])
        for j in range(n):
            angle = th_x[j]
            X_eval = x_limit*math.cos(angle)
            Y_eval = (y_limit-1)*math.sin(angle)+1
            res1, table = dickman((math.log(abs(eval_F(X_eval, Y_eval, f, len(f)-1)))-(mu_x[i]-c))/log_B, table)
            res2, table = dickman((math.log(abs(eval_F(X_eval, Y_eval, g, len(g)-1))))/log_B, table)
            res += mu_w[i]*th_w[j]*res1*res2*Y

    return res

def alpha_score(f, primes):
    E, F = 0, 0
    evals = [0]*primes[-1]
    tmp = [0]*len(f)
    for j in range(len(f)): tmp[j] = evaluate(f, j)
    for q in range(1, len(f)):
        for k in range(len(f)-1, q-1, -1): tmp[k] -= tmp[k-1]

    eval = tmp[0]
    for k in range(primes[-1]):
        evals[k] = eval
        eval += tmp[1]
        for q in range(1, len(f)-1): tmp[q] += tmp[q+1]

    f_prime = get_derivative(f)
    upto = len(f)+10
    baseline_term = 0
    for p in primes:
        log_p = math.log(p)
        baseline_term += log_p/(p-1)
        Ep,Fp = 0,0
        ramified_roots = []
        ramified = False
        for r in range(p):
            if not evals[r]%p:
                if not eval_mod(f_prime, r, p):
                    ramified = True
                    ramified_roots.append(r)
                    Ep += 1/p
                    Fp += 1/p
                else:
                    Ep += 1/(p-1)
                    Fp += (p+1)/((p-1)**2)

        if ramified:
            tmp2 = p
            for i in range(2, upto):
                new = []
                for r in ramified_roots:
                    if not eval_mod(f, r, tmp2*p):
                        for k in range(p): new.append(r+k*tmp2)
                        Ep += 1/tmp2
                        Fp += (2*i-1)/tmp2
                tmp2 *= p
                ramified_roots = new.copy()

            Ep += len(ramified_roots)/(pow(p, upto-2)*(p-1))
            Fp += len(ramified_roots)*(upto*upto+2*upto/(p-1)+(p+1)/((p-1)**2))/tmp2

        if not f[0]%p:
            frev = f[::-1]
            fdrev = get_derivative(frev)
            if eval_mod(fdrev,0,p):
                Ep += 1/(p-1)
                Fp += (p+1)/((p-1)**2)
            else:
                Ep += 1/p
                Fp += 1/p
                ramified_roots = [0]
                tmp2 = p
                for i in range(2, upto):
                    new = []
                    for r in ramified_roots:
                        if not eval_mod(frev, r, tmp2*p):
                            for k in range(p): new.append(r+k*tmp2)
                            Ep += 1/tmp2
                            Fp += (2*i-1)/tmp2
                    tmp2 *= p
                    ramified_roots = new.copy()

                Ep += len(ramified_roots)/(pow(p, upto-2)*(p-1))
                Fp += len(ramified_roots)*(upto*upto+2*upto/(p-1)+(p+1)/((p-1)**2))/tmp2
                
        tmpE, tmpF = p*Ep/(p+1), p*Fp/(p+1)
        E += tmpE*log_p
        F += (tmpF-tmpE*tmpE)*log_p*log_p
        
    k = 2*E-F/2
    lbd = F/2-E
    return E-baseline_term, k, lbd, baseline_term
    
def local_opt(f, g, B):
    m0, m1 = -g[1], g[0]

    d = len(f)-1
    if d <= 5:
        poly = []
        for i in range(3): poly.append(f[i]*binom(2-i, d-i))
        poly = poly_prod(poly, poly)
        poly2 = get_derivative(poly)
        zeros = get_complex_roots(poly2)
        min,mink = f[2]**2,0
        for r in zeros:
            if r.imag == 0:
                if evaluate(poly, math.ceil(r.real)) < min:
                    min, mink = evaluate(poly, math.ceil(r.real)), math.ceil(r.real)
                if evaluate(poly, math.floor(r.real)) < min:
                    min, mink = evaluate(poly, math.floor(r.real)), math.floor(r.real)
    else:
        poly = []
        for i in range(4): poly.append(f[i]*binom(3-i, d-i))
        zero = get_complex_roots(poly)
        min, mink = abs(f[3]), 0
        for r in zero:
            if r.imag == 0:
                if abs(evaluate(poly, math.ceil(r.real))) < min:
                    min, mink = abs(evaluate(poly, math.ceil(r.real))), math.ceil(r.real)
                if abs(evaluate(poly, math.floor(r.real))) < min:
                    min, mink = abs(evaluate(poly, math.floor(r.real))), math.floor(r.real)

    f = shift(f, mink)
    m0 -= mink*m1

    k, u, v, iteration = 1, 1, 1, 0
    _, s = get_sieve_region(f, B)
    min_n = get_Lnorm(poly_prod(f, f), s, B)
    while iteration < 500 and (k > 0 or u > 0 or v > 0):
        F = shift(f, -k)
        tmp = get_Lnorm(poly_prod(F, F), s, B)
        flag = False
        if tmp < min_n: f, min_n, k, flag, m0 = F.copy(), tmp, k<<1, True, m0+k*m1

        F = shift(f, k)
        tmp = get_Lnorm(poly_prod(F, F), s, B)
        if tmp < min_n: f, min_n, k, flag, m0 = F.copy(), tmp, k<<1, True, m0-k*m1
        if flag and not u: u = 1
        if flag and not v: v = 1
        elif not flag: k >>= 1

        if u:
            flag = False
            F = f.copy()
            F[-3] += u*m1
            F[-2] -= u*m0
            tmp = get_Lnorm(poly_prod(F, F), s, B)
            if tmp < min_n: f, min_n, u, flag = F.copy(), tmp, u<<1, True

            F = f.copy()
            F[-3] -= u*m1
            F[-2] += u*m0
            tmp = get_Lnorm(poly_prod(F, F), s, B)
            if tmp < min_n: f, min_n, u, flag = F.copy(), tmp, u<<1, True
            if flag and not k: k = 1
            if flag and not v: v = 1
            elif not flag: u >>= 1

        if v:
            flag = False
            F = f.copy()
            F[-2] += v*m1
            F[-1] -= v*m0
            tmp = get_Lnorm(poly_prod(F, F), s, B)
            if tmp < min_n: f, min_n, v, flag = F.copy(), tmp, v<<1, True

            F = f.copy()
            F[-2] -= v*m1
            F[-1] += v*m0
            tmp = get_Lnorm(poly_prod(F, F), s, B)
            if tmp < min_n: f, min_n, v, flag = F.copy(), tmp, v<<1, True
            if flag and not k: k = 1
            if flag and not u: u = 1
            elif not flag: v >>= 1

        iteration += 1
        _, s = get_sieve_region(f, B)

    return f, m0

def prime_combinations_with_indices(Q, l, B):
    n = len(Q)
    indices = [0] * l  # reuse buffer

    def backtrack(start, depth, product):
        if depth == l:
            yield tuple(indices)
            return
        for i in range(start, n - (l - depth) + 1):
            p = Q[i]
            if product * p >= B:
                break  # sorted Q means all further i will be too large
            indices[depth] = i
            yield from backtrack(i + 1, depth + 1, product * p)

    yield from backtrack(0, 0, 1)

def compute_e(m0, root_used, NB_ROOTS, prod, a_d, n, d):
    e = []
    tmp_m = m_mu(m0, root_used, [0]*NB_ROOTS, NB_ROOTS)
    base = poly(tmp_m, prod, a_d, n, d)[1]%prod
    for i in range(NB_ROOTS):
        line = [0]*d
        for j in range(d):
            if not i:
                tmp_m = m_mu(m0, root_used, [j]+[0]*(NB_ROOTS-1), NB_ROOTS)
                line[j] = poly(tmp_m, prod, a_d, n, d)[1]%prod
            elif j:
                tmp_m = m_mu(m0, root_used, [0]*i+[j]+[0]*(NB_ROOTS-i-1), NB_ROOTS)
                line[j] = (poly(tmp_m, prod, a_d, n, d)[1]-base)%prod
        e.append(line)

    return e

def compute_f(n, a_d, m0, d, prod, root_used, NB_ROOTS, e):
    f0 = (n-a_d*pow(m0, d))/(prod*prod*pow(m0, d-1))
    f = []

    for i in range(NB_ROOTS):
        line = [0]*d
        for j in range(d):
            line[j] = -(a_d*d*root_used[i][j]/pow(prod, 2)+e[i][j]/prod)
        f.append(line)

    return f, f0

def create_first_array(NB_ROOTS, f0, f, d):
    vec = [0]*(NB_ROOTS>>1)
    array1 = []

    while vec[-1] < d:
        U = (f0 + sum([f[j][vec[j]] for j in range(NB_ROOTS>>1)]))%1

        if not len(array1) or U > array1[-1][0]: array1.append([U, [i for i in vec]])
        else:
            tmp_a = 0
            tmp_b = len(array1)-1
            tmp = (tmp_a+tmp_b)>>1
            while tmp_a <= tmp_b:
                if array1[tmp][0] > U: tmp_b = tmp-1
                else: tmp_a = tmp+1
                tmp = (tmp_a+tmp_b)>>1
            array1.insert(tmp_a, [U, [i for i in vec]])
        vec[0] += 1
        for j in range(len(vec)-1):
            if vec[j] == d:
                vec[j] = 0
                vec[j+1] += 1
            else: break

    return array1

def create_second_array(NB_ROOTS, len_vec, d, f):
    vect = [0]*(NB_ROOTS-len_vec)
    array2 = []
    while vect[-1] < d:
        U = -sum([f[len_vec+j][vect[j]] for j in range(len(vect))])%1

        if not len(array2) or U > array2[-1][0]: array2.append([U, [i for i in vect]])
        else:
            tmp_a = 0
            tmp_b = len(array2)-1
            tmp = (tmp_a+tmp_b)>>1
            while tmp_a <= tmp_b:
                if array2[tmp][0] > U: tmp_b = tmp-1
                else: tmp_a = tmp+1
                tmp = (tmp_a+tmp_b)>>1
            array2.insert(tmp_a, [U, [i for i in vect]])
        vect[0] += 1
        for j in range(len(vect)-1):
            if vect[j] == d:
                vect[j] = 0
                vect[j+1] += 1
            else: break

    return array2

def m_mu(m0,roots,vec,l):
    return m0 + sum([roots[i][vec[i]] for i in range(l)])

def poly(m0, m1, a_d, n, d):
    c = [a_d]
    r = [n]
    for i in range(d-1, -1, -1):
        r.append((r[-1]-c[-1]*pow(m0, i+1))//m1)
        delta = -r[-1]*m1*invmod(m1,pow(m0, i))%(m1*pow(m0, i))
        c.append((r[-1]+delta)//pow(m0, i))
    for i in range(1, len(c)):
        if c[i] > m0//2:
            c[i] -= m0
            c[i-1] += m1
        elif c[i] < -m0//2:
            c[i] += m0
            c[i-1] -= m1
    return c

def root_sieve(f, g, primes, U):
    array = [0]*(U<<1)
    for p in primes:
        k = 1
        P = p
        while P < primes[-1]:
            for x in range(P):
                eval1, eval2 = eval_mod(f, x, P), eval_mod(g, x, P)
                for u in range(P):
                    if (eval1+u*eval2)%P == 0:
                        for z in range((u+U)%P, len(array), P): array[z] += math.log(p)/(pow(p, k-1)*(p+1))
            k += 1
            P *= p

    max,maxi = array[0],0
    for i in range(1, len(array)):
        if array[i] > max: max, maxi = array[i], i
    u = maxi-U
    print(u)

    f[-1] += u*g[-1]
    f[-2] += u*g[-2]

    return f

def root_sieve2(f, g, primes, U, V):
    array = [[0 for _ in range(V<<1)] for __ in range(U<<1)]
    for p in primes:
        k = 1
        P = p
        while P < primes[-1]:
            for x in range(P):
                eval1, eval2 = eval_mod(f, x, P), eval_mod(g, x, P)
                for u in range(P):
                    if eval2%p:
                        v = -(u*x+eval1*invmod(eval2, P))%P
                        for i in range((u+U)%P, len(array), P):
                            for j in range((v+V)%P, len(array[0]), P): array[i][j] += math.log(p)/(pow(p, k-1)*(p+1))
            k += 1
            P *= p

    max, maxu, maxv = array[0][0], 0, 0
    for i in range(U<<1):
        for j in range(V<<1):
            if array[i][j] > max: max, maxu, maxv = array[i][j], i, j

    u, v = maxu-U, maxv-V
    print(u, v)

    f[-1] += v*g[-1]
    f[-2] += v*g[-2]+u*g[-1]
    f[-3] += u*g[-2]

    return f


def poly_search(n, primes, nb_roots, prime_bound, c, M, d, NB_POLY_COARSE_EVAL, NB_POLY_PRECISE_EVAL):
    print("Starting polynomial search")


    return Kleinjung_poly_search(n, primes, nb_roots, prime_bound, c, M, d, NB_POLY_COARSE_EVAL,
                                     NB_POLY_PRECISE_EVAL)


# Kleinjung first polynomial search algorithm
def Kleinjung_poly_search(n, primes, NB_ROOTS, PRIME_BOUND, MULTIPLIER, M, d, NB_POLY_COARSE_EVAL, NB_POLY_PRECISE_EVAL):
    t1 = datetime.now()
    P = []
    polys = []
    for p in primes:
        if p > PRIME_BOUND: break
        if p%d == 1: P.append(p)
    a_d = MULTIPLIER
    if d >= 4: admax = round(pow(pow(M, 2*d-2)/n, 1/(d-3)))
    else: admax = M

    cpt = 0
    avg = 0
    while a_d < admax and cpt < NB_POLY_COARSE_EVAL:
        tmp = a_d
        for p in primes:
            while not tmp%p: tmp//= p

        if tmp > 1: # If a_d is not primes[-1] smooth
            a_d += MULTIPLIER
            continue

        mw = math.ceil(pow(n/a_d, 1/d))
        ad1max = round(M*M/mw)
        if d > 2: ad2max = pow(pow(M, 2*d-6)/pow(mw, d-4), 1/(d-2))
        else: ad2max = M

        Q = []
        roots = []
        for p in P:
            if not a_d%p: continue

            f = [a_d]+[0]*d
            f[-1] = (-n)%p # Construct polynomial a_d*x^d - n (mod r)
            root = find_roots_poly(f, p)
            if len(root) > 0:
                Q.append(p)
                roots.append(root)
      #  print("roots: ",roots)
        if len(roots) >= NB_ROOTS:

            combinations = prime_combinations_with_indices(Q, NB_ROOTS, ad1max)

            for set in combinations:
                Q_used = []
                prod = 1
                for i in range(NB_ROOTS):
                    Q_used.append(Q[set[i]])
                    prod *= Q[set[i]]

                root_used = [roots[set[i]] for i in range(NB_ROOTS)]
                for i in range(NB_ROOTS): # Do some CRT
                    x = prod//Q_used[i]
                    tmp2 = x*invmod(x, Q_used[i])
                    for j in range(d): 
                        root_used[i][j] = root_used[i][j]*tmp2%prod

                m0 = mw+(-mw)%prod
                e = compute_e(m0, root_used, NB_ROOTS, prod, a_d, n, d)
                f, f0 = compute_f(n, a_d, m0, d, prod, root_used, NB_ROOTS, e)

                epsilon = ad2max/m0
                array1 = create_first_array(NB_ROOTS, f0, f, d)
                len_vec = NB_ROOTS>>1
                array2 = create_second_array(NB_ROOTS, len_vec, d, f)
                
                min = 0
                for j in range(len(array2)):
                    while min < len(array1) and array2[j][0]-epsilon > array1[min][0]: min += 1
                    if min == len(array1): break
                    z = min
                    while z < len(array1) and abs(array2[j][0]-array1[z][0]) < epsilon:
                       # print("root_used: ",root_used)
                        tmp = [poly(m_mu(m0, root_used, array1[z][1]+array2[j][1], NB_ROOTS), prod, a_d, n, d),
                               m_mu(m0, root_used, array1[z][1]+array2[j][1], NB_ROOTS),
                               prod]
                        cpt += 1
    
                        tmp[0], tmp[1] = local_opt(tmp[0], [prod,-tmp[1]], primes[-1])
                        _, s = get_sieve_region(tmp[0], primes[-1])
                        L_score = get_Lnorm(poly_prod(tmp[0], tmp[0]), s, primes[-1])
                        avg += L_score
                        if not len(polys):
                            tmp.append(L_score)
                            polys.append(tmp)
                        else:
                            if L_score < polys[-1][3]:
                                tmp.append(L_score)
                                a = 0
                                b = len(polys)-1
                                tmpu = (a+b)>>1
                                while a <= b:
                                    if polys[tmpu][3] < L_score: a = tmpu+1
                                    else: b = tmpu-1
                                    tmpu = (a+b)>>1
                                polys.insert(a, tmp)
                            elif len(polys) < NB_POLY_PRECISE_EVAL:
                                tmp.append(L_score)
                                polys.append(tmp)
                            if len(polys) > NB_POLY_PRECISE_EVAL: del polys[-1]

                        if cpt >= NB_POLY_COARSE_EVAL:

                            t2 = datetime.now()
                            print("")
                            print("Polynomial search done in "+format_duration(t2-t1)+".\n")
                            print(str(cpt)+" polynomials created, "+str(len(polys))+" kept for ranking")
                            print("Average L2 score = "+str(avg/cpt))
                            print("Ranking polynomials")

                            return select_best_poly_candidate(polys, primes)

                        z += 1
        a_d += MULTIPLIER

def select_best_poly_candidate(polys, primes):
  #  print("polys: ",polys)
    best_poly = None
    best_E = None
    table = get_dickman_table(10)

    for i in range(len(polys)):
        x_limit,s = get_sieve_region(polys[i][0], primes[-1])
        alpha = alpha_score(polys[i][0], primes[:300])
        polys[i].append(alpha)
        E_score = get_Epscore(polys[i][0], [polys[i][2], -polys[i][1]], alpha, primes[-1], x_limit, round(primes[-1]/math.sqrt(s)), table)
        if best_E == None or E_score > best_E:
            best_poly = polys[i]
            best_E = E_score

    return best_poly

def practical_d(n):
    min,mind = None,None
    for d in range(2, 10):
        f = d*math.log(d)+math.sqrt(pow(d*math.log(d), 2)+4*math.log(n)*math.log(math.log(n)/(d+1))/d+1)
        if min == None or f < min: min, mind = f, d
    return mind
    
def initialize(n):
    
    d = int(pow(3*math.log(n)/math.log(math.log(n)), 1/3))
    first_B = int(math.exp(pow((8/9)*math.log(n), 1/3)*pow(math.log(math.log(n)), 2/3)))
    B = 10*first_B#//10
    
    return d, B
    
def initialize_2(f_x, m1, d, primes, leading_coeff):
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
    
    return pairs_used, R_p, logs, divide_leading, pow_div, B_prime
    
def initialize_3(n, f_x, f_prime, const1, leading_coeff):
    k = 3*n.bit_length()
    Q = []
    q = const1+1
    if not q%2: q += 1
    while len(Q) < k:
        if is_prime(q) and leading_coeff%q:
            if not n%q: return q, n//q
            else:
                tmp = find_roots_poly(f_x, q)
                for r in tmp:
                    if eval_mod(f_prime, r, q): Q.append([q, r])
        q += 2
        
    return Q, k
    
def create_smooth_primes_base(B):
    m = math.isqrt(B)+1
    
    liste = [True, False]*(B>>1)
    base=[2]
    
    for k in (k for k in range (3, m, 2) if liste[k-1]):
        base.append(k)
        for i in range (k*k, B, k<<1):
            liste[i-1]=False

    for k in (k for k in range (m, B) if liste[k-1]):
        base.append(k)

    return base 



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


def NFS(n,z,y,x,polyval,m1_a,primeslist):
    LARGE_PRIME_CONST=10000
    BLOCK_SIZE=8
    NB_POLY_COARSE_EVAL=100
    NB_POLY_PRECISE_EVAL=10
    PRIME_BOUND=300
    NB_ROOTS=2
    MULTIPLIER=1
    
    n = int(n)
    d, B = initialize(n)
    d=2
    primes = create_smooth_primes_base(B)
    for p in primes:
        if not n%p: return p, n//p

    prod_primes = math.prod(primes)
    const1, const2 = LARGE_PRIME_CONST*primes[-1], LARGE_PRIME_CONST*primes[-1]*primes[-1]

  #  f_x,m0,m1,tmp,_ = poly_search(n, primes, NB_ROOTS, PRIME_BOUND, MULTIPLIER,
   #                                                             int(pow(n, 1/(d+1))), d, NB_POLY_COARSE_EVAL,
    #                                                             NB_POLY_PRECISE_EVAL)

    f_x=[z,y,-polyval]
    m0=-m1_a
    m1=z
    f_x, m0, M = evaluate_polynomial_quality(f_x, B, m0, m1, primes)
   # M=2000
    print("poly search completed, parameters : m0 = "+str(m0)+" ; m1 = "+str(m1)+" ; d = "+str(d)+" M: "+str(M))
    leading_coeff = f_x[0]
    f_prime = get_derivative(f_x)
    g = [1, f_x[1]]
    x_test = symbols('x_test')
    if d==3:
        p =Poly(f_x[0]*x_test**3+f_x[1]*x_test**2+f_x[2]*x_test+f_x[3],x_test)
    if d==2:
        p =Poly(f_x[0]*x_test**2+f_x[1]*x_test+f_x[2],x_test)
    else:
        p=None
    q =Poly(m1*x_test-m0,x_test)
   # rx=resultant(p,q)
   # print("Resultant: ",rx)
  #  if rx%n !=0:
     #   print("Catastrophic error")
      #  sys.exit()
    for i in range(2, len(f_x)): 
       # print(str(f_x[i])+" pow(leading_coeff, i-1): "+str(pow(leading_coeff, i-1))+" leading_coeff: "+str(leading_coeff))
        g.append(f_x[i]*pow(leading_coeff, i-1))
    #print("g: "+str(g))
    g_prime = get_derivative(g)
    g_prime_sq,g_prime_eval = div_poly(poly_prod(g_prime, g_prime), g), pow(leading_coeff, d-2, n)*eval_F(m0, m1, f_prime, d-1)%n
    inert_set = []
    for p in primes:
       # print("m1: "+str(m1%p)+" p: "+str(p)+" leading coeff: "+str(leading_coeff%p)+" irreducibility: "+str(irreducibility(g, p))+" g: "+str(g))
        if m1%p and leading_coeff%p and irreducibility(g, p): 
            
            inert_set.append(p)
    if len(inert_set)==0:
        print("fail")
        sys.exit()
    pairs_used, R_p, logs, divide_leading, pow_div, B_prime = initialize_2(f_x, m1, d, primes, leading_coeff)
        
    Q, k = initialize_3(n, f_x, f_prime, B, leading_coeff)
        
    print("components of factor base : ("+str(len(primes))+","+str(B_prime)+","+str(k)+","+str(len(divide_leading))+")")
    print("starting with "+str(len(pairs_used))+" free relations")
    pairs_used, V = find_relations(f_x, leading_coeff, g, primes, R_p, Q, B_prime, divide_leading,prod_primes, pow_div, pairs_used, const1, const2, logs, m0, m1,M)
    print("")
    print("sieving complete, building matrix...")
    matrix = build_sparse_matrix(pairs_used, primes, R_p, Q, divide_leading)
    matrix = transpose_sparse(matrix, V)
    print("matrix built "+str(len(matrix))+"x"+str(len(pairs_used))+" reducing...")
    matrix, pairs_used = reduce_sparse_matrix(matrix, pairs_used)
    print("matrix built "+str(len(matrix))+"x"+str(len(pairs_used))+" finding kernel...")
    time_1 = datetime.now()
    while True:
        null_space = block_lanczos(matrix, len(pairs_used), BLOCK_SIZE)
                
        print(str(len(null_space))+" kernel vectors found")
        for vec in null_space:
            vec = convert_to_binary_lanczos(vec, len(pairs_used))
            flag, res = compute_factors(pairs_used, vec, n, primes, g, g_prime, g_prime_sq, g_prime_eval,
                                                             m0, m1, leading_coeff, inert_set, M,
                                                             time_1)
            if flag: return res








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
        
def QS(n,factor_list,sm,xlist,flist,xy_list,xy_f_list,x_list,x_f_list,primeslist2):#,jsymbols,testl,primeslist2,disc1_squared_list):#,disc_sr_list,pval_list,pflist):
    g_max_smooths=base+2#+qbase
    if len(sm) > g_max_smooths*10000000: 
        del sm[g_max_smooths:]
        del xlist[g_max_smooths:]
        del flist[g_max_smooths:]  
    M2 = build_matrix(factor_list, sm, flist,xy_f_list,x_f_list,primeslist2)#,pflist)
    null_space=solve_bits(M2,factor_list,len(sm))
    f1,f2=extract_factors(n, sm, xlist, null_space,xy_list,x_list)#,disc_sr_list,pval_list,pflist)
    if f1 != 0:
        print("[SUCCESS]Factors are: "+str(f1)+" and "+str(f2))
        return f1,f2   
   # print("[FAILURE]No factors found")
    return 0,0

def extract_factors(N, relations, roots, null_space,xy_list,x_list):#,disc_sr_list,pval_list,pflist):
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
                prod_left *= roots[idx]
                prod_right *= relations[idx]
             #   xy*=xy_list[idx]
                x*=x_list[idx]
                print("polyval:  "+str(relations[idx])+" disc constant "+str(x_list[idx]))
            idx += 1
        prod_right=x
        sqrt_right = math.isqrt(x)
        sqrt_left = math.isqrt(prod_left)#prod_left
        if sqrt_left**2 != prod_left:
            print("horrible error")
        print(" polyval sqrt: "+str(sqrt_left)+" disc constant sqrt: "+str(sqrt_right))#+" zx*zxy: "+str(x))
        ###Debug shit, remove for final version
        sqr1=prod_left%N 
        sqr2=prod_right%N
        if sqrt_right**2 != prod_right:
            print("not a square in the integers")
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

def build_matrix(factor_base, smooth_nums, factors,xy_f_list,x_f_list,primeslist2):#,pflist):
    fb_map = {val: i for i, val in enumerate(factor_base)}

    ind=1

    M2=[0]*((base+2)*2)#+qbase)#+2+qbase)
    for i in range(len(smooth_nums)):
        for fac in factors[i]:
            idx = fb_map[fac]
            M2[idx] |= ind
        ind = ind + ind       

  #  offset=(base+2)-1
  #  ind=1
    
 #   for i in range(len(xy_f_list)):
 #       j=0
 #       while j < len(xy_f_list[i]):
 #           fac=xy_f_list[i][j]
 #           if fac == -1 or fac ==0:
 #               M2[j+offset] |= ind
 #           j+=1
 #       ind = ind + ind

    offset=(base+2)-1
    ind=1
    for i in range(len(x_f_list)):
        for fac in x_f_list[i]:
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
        hmap_p.append(array.array('q',[-1]*prime))
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
                for root in roots:
                    x = root#lift_root(1,y2,n*k,root,prime,1)
                  #  if prime == 3:
                      #  print("roots: "+str(roots)+" y: "+str(y2)+" k: "+str(k)+" x: "+str(x))
                    if (k*x**2+y2*x-n)%prime==0:#**2==0: ##To do: Bug here if we have a single root from solve_quadratic_congruence..
                        hmap_p[k-1][x]=(k*x)+y2
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
def generate_modulus(n,primeslist,seen,tnum,close_range,too_close,LOWER_BOUND_SIQS,UPPER_BOUND_SIQS,tnum_bit,hmap,quad):
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

cdef construct_interval(list ret_array,partials,n,primeslist,hmap,hmap2,large_prime_bound,primeslist2,small_primeslist):
    i=0
  #  while i < len(hmap):
      #  print("prime: "+str(primeslist[i])+" "+str(hmap[i]))
      #  i+=1


    primelist_f=copy.copy(primeslist)
    primelist_f.insert(0,len(primelist_f)+1)
    primelist_f=array.array('q',primelist_f)

    primelist2_f=copy.copy(primeslist2)
    primelist2_f.insert(0,len(primelist2_f)+1)
    primelist2_f=array.array('q',primelist2_f)
    cdef Py_ssize_t size
    primelist=copy.copy(primeslist)
    primelist.insert(0,2) ##To do: remove when we fix lifting for powers of 2
    primelist.insert(0,-1)
    print("[i]Filtering Quadratic Coefficients (quad_size) (to do: can be saved to disk for re-use)")
    valid_quads,valid_quads_factors,qival=filter(primelist_f,n,1,quad_sieve_size+1)
    print("[i]Filtering interval indices (lin_size) (to do: can be saved to disk for re-use)")
    start=round(n**0.25)

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
    sfound=0
    close_range=10
    too_close=5
    LOWER_BOUND_SIQS=1
    UPPER_BOUND_SIQS=4000
   # tnum=int(((n)**0.6) / lin_sieve_size)
    seen=[]
    threshold = int(math.log2((lin_sieve_size)*math.sqrt(abs(n))) - thresvar)
    if threshold < 0:
        threshold = 1
    while 1:
       # new_mod,cfact,indexes=generate_modulus(n,primeslist,seen,tnum,close_range,too_close,LOWER_BOUND_SIQS,UPPER_BOUND_SIQS,bitlen(tnum),hmap,1)
       # if new_mod ==0:
       #     return 0,0
       # col=[]
       # i=0
       # while i < len(indexes):
       #     idx=indexes[i]
       #     if primeslist[idx] != cfact[i]:
       #         print("fatal error")
       #         sys.exit()

       #     prime=primeslist[idx]
       #     col.append(hmap[idx])
       #     i+=1

        i=1
        while i < 1000_0000_000:
            local_factors, value,seen_primes,seen_primes_indexes = factorise_fast(start+i,primelist_f)
            if value != 1:
                i+=1
                continue
            x=start+i
            start_xy=2**(round(keysize*0.60))
            start_y=start_xy-x 
            mult_list=retrieve(hmap,primeslist,x,start_y,n)          
            q=0
            while q < len(mult_list):  
                if q+1 not in valid_quads:
                    q+=1
                    continue

                z=q+1
                xy=start_y+x*z
                local_factors, value,seen_primes,seen_primes_indexes = factorise_fast(xy,primelist_f)
                if value != 1:
                    q+=1
                    continue
                o=1
                while o < lin_sieve_size:
                    k=1#q+1

                    y=(xy*o)-(x*z)
                   # y=(start_xy*o)-(x*z)
                    poly_val=(z*x**2+y*x-n)%(n-1)
                  #  k=((z*x**2+y*x)-poly_val)//n
                    k=1
                    if mult_list[q][o]<bitlen(poly_val)*0.6 or lival[o-1] !=1:
                        o+=1
                        continue
                    

                 #   print("checking")

                    found=[]
                    g=0
                    total=1
                    
                       
  
                    if poly_val==0 or k ==0:
                        o+=1
                        continue
                    local_factors, value,seen_primes,seen_primes_indexes = factorise_fast((poly_val*z)%(n-1),primelist_f)
                    if value == 1:
                       # disc1_squared=(y**2+4*(n*k*z+(poly_val*z)))%(n-1)
                       # disc1=math.isqrt(disc1_squared)
                       # if disc1**2 != disc1_squared:
                         #   print("fatal error ergh")

                        poly_val2=(n*k*z+(poly_val*z))  #note: factorization for this is y+disc1 and y-disc1

                        factors_part1=z*x
                        factors_part2=z*x+y
                        factors_part3=poly_val*z
                        all_parts=factors_part1*factors_part2*factors_part3
                     #   if all_parts != (poly_val2*poly_val*z)%(n-1):
                       #     print("fatal error qq")
                   #     if (4*poly_val2)%((2*z*x))!=0:
                       #     print("error")

                    #    if (4*poly_val2)%(2*(z*x+y))!=0:
                         #   print("error2")
                        local_factors2, value2,seen_primes2,seen_primes_indexes2 = factorise_fast((poly_val*z)%(n-1),primelist_f)
                        test=math.isqrt(value2)
                        if value2 == 1 or test**2==value2:
                            if (poly_val*z)%(n-1) not in coefficients and (poly_val*z)%(n-1) != poly_val2:# and local_factors2 not in factors:


                                local_factors4, value4,seen_primes4,seen_primes_indexes4 = factorise_fast((poly_val2),primelist_f)
                                if value4 !=1:
                                    o+=1
                                    continue
                                x_test = symbols('x_test')
                                p1 =Poly(z*x_test**2+y*x_test-poly_val,x_test)
                                q1 =Poly(z*x_test+(z*x+y),x_test)
                                m1=(z*x+y)
                                rx=resultant(p1,q1)
                                print("rx: ",rx)
                                if rx%n == 0:# and z > 10:
                                    NFS(n,z,y,x,poly_val,m1,primelist)
                                smooths.append((poly_val*z)%(n-1))
                                factors.append(local_factors2)
                                coefficients.append((poly_val*z)%(n-1))
                              #  symbols=[]
                              #  r=0
                               # while r < len(primeslist2):
                                  #  symbols.append(jacobi(z*x+y,primeslist2[r]))

                                  #  r+=1
                                  #  print("super fatal error2")
                                  #  sys.exit()
                                xy_list.append(z*x+y)
                                xy_f_list.append(symbols)
                                x_list.append(poly_val2)
                                x_f_list.append(local_factors4)
                                sfound+=1
                                if g_debug == 1:
                                    print("Smooths #: "+str(len(smooths))+" y: "+str(y)+" zx: "+str(factors_part1)+" zx2 "+str(z*x+y)+" (poly_val*z): "+str(factors_part3)+" final smooth: "+str(all_parts)+" intrvl ind: "+str(o)+" Factors: "+str(local_factors2)+" z: "+str(z))#+" seen_primes: "+str(valid_ind_factors[o1]))#+" seen_primes2: "+str(seen_primes2))#+" test_poly_val: "+str(bitlen(test_poly_val))+" test_zxy: "+str(test_zxy)+" test_zxy_current: "+str(test_zxy_curent))
                                else: 
                                    print("Smooths #: "+str(len(smooths)))
                                if sfound >(base+2)*2:#+qbase:
                                    f1,f2=QS(n,primelist,smooths,coefficients,factors,xy_list,xy_f_list,x_list,x_f_list,primeslist2)
                                    if f1 !=0:
                                        sys.exit()
                                    sfound=0
                        else:
                            print("FATAL ERROR, THIS ONE SHOULD NEVER FAIL BECAUSE z,x,x2 AND poly_val MUST FACTORIZE WHEN ARRIVING HERE:"+str(value2)+" k: "+str(k)+" poly_val: "+str(poly_val)+" o: "+str(o)+" x: "+str(x)+" xy: "+str(xy))
                           # sys.exit()
                   ## else:
                     ##   print("wtf")

                    o+=1
                q+=1
            i+=1

        return
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