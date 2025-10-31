###Author: Essbee Vanhoutte
###WORK IN PROGRESS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

###notes: use with python3 polarbearalg_v18.py -keysize 14

####I am working on my own custom number field sieve implementation using reducible quadratic polynomials.

##See github repo readme for more info



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
p_amount=20 #amount of primes in factor base
p_mod_amount=20 #amount of primes in the modulus
g_enable_custom_factors=0
g_p=107
g_q=41
quad_char_amount=5


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

def launch(n,primeslist1,quad_base):
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

    while z < workers:
        p=multiprocessing.Process(target=find_comb, args=(primeslist1,n,procnum,return_dict,hmap,quad_base,hmap2)) 
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
    modi=inverse(4*n,prime)
    while iN < prime:
        ja= jacobi(iN,prime )
        if ja ==1:
            root=tonelli(iN,prime)
            
            if root > prime // 2:
                root=(prime-root)%prime

            s=(root**2*modi)%prime



            try:
                c=hmap_p[str(s)]
                c.append(root)
                c.append((-root)%prime)
            except Exception as e:
                    c=hmap_p[str(s)]=[root,(-root)%prime]

            try:
                c=hmap_p2[str(root)]
                c.append(s)
            except Exception as e:
                    c=hmap_p2[str(root)]=[s]

        iN+=1  
    return hmap_p,hmap_p2

def create_hashmap(lists,n):
    i=0
    hmap=[]
    hmap2=[]
    while i < len(lists):
        hmap_p,hmap_p2=solve_roots(lists[i],n)
        hmap.append(hmap_p)
        hmap2.append(hmap_p2)
      #  print(str(lists[i])+" "+str(hmap_p))
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


def solve_quadratic_congruence(a, b, c, p):
    if a == 0: #FIX THIS SHIT
      #  print("error, linear congruence, fix me")
        return 0
    else:
        a_inv = inverse(a, p)
        ba = (b * a_inv) % p
        ca = (c * a_inv) % p
        b_div_2 = (ba * inverse(2, p)) % p
        alpha = (pow_mod(b_div_2, 2, p) - ca) % p
     #   print("alpha: ",alpha)
        if jacobi(alpha,p)==-1:
            return -1
        elif jacobi(alpha,p)==0:
            return 0
      #  print("alpha: "+str(alpha))
        y = tonelli (alpha, p)
        if y==None:
            print("should never happen")
            return -1
        x1 = (y - b_div_2) % p
       # x2 = (p - y - b_div_2) % p

        return x1

def gen_comb(collected,mod,z,n,factor_base,ret_array,quad_base,hmap2):
    ret_array[0].append(mod)
    ret_array[0].append([0])
    ret_array[1].append(mod)
    ret_array[1].append([0])
   # ret_array[2].append(mod)
    ret_array[2].append([[0]])
    enum=[]
    i=0
    while i < len(collected):
        enum.append(collected[i+1])
        i+=2

    for idx in enumerated_product(*enum):
        y0_temp=0
        #quad=quad_co
        l=0
        while l < len(idx):
            y0_temp+=enum[l][idx[l]]
            l+=1
        y0_temp%=mod
        i=0
        while i < 1:
            
            y0=y0_temp
            x2=get_root(n,mod,y0,z) 
            x=x2*z
            x%=mod
            n_red=n%mod
            eq=x**2-(n)*z
            eq2=z*x2**2-(n)
            if eq == 0:
                i+=1
                continue
            
            eq_o=eq//mod

            test=(x**2-n*z)%eq_o
            y1=abs(eq2**0.5)
            test2=n_red+mod*i
            lside=z*(x2**2)
            if abs(eq2**0.5)%1 ==0 and lside**0.5%1==0:
                print("y0: "+str(y0)+" x: "+str(x)+" eq: "+str(eq)+" x2: "+str(x2)+" eq2: "+str(abs(eq2)**0.5)+" mod: "+str(mod)+" i: "+str(i)+" z: "+str(z))
           # if n%eq2 != test2%eq2:
             #   i+=1 
             #   continue

            if abs(eq2**0.5)%1 != 0 or lside**0.5%1!=0:# or test != 0:
                i+=1
                continue 

            jac=[]
            for prime in quad_base:
                temp_eq=mod+n*z
                jac.append(jacobi((temp_eq%mod)%prime,prime))
            ret_array[0][-1].append(eq)
            ret_array[1][-1].append(x)
            ret_array[2][-1].append(jac)
           
            new_root=lside**0.5#int((z*(x2**2))**0.5)
            print(new_root)
            factor=gcd(new_root+y1,n)
            if factor != 1 and factor != n:
                print("Factored N: "+str(factor)+" and: "+str(n//factor))
                sys.exit()
            i+=1

def solve_lin_con(a,b,m):
    ##ax=b mod m
    g=math.gcd(a,m)
    #print(g)
    if b%g:
        print("no sols")
    a,b,m = a//g,b//g,m//g
    return pow(a,-1,m)*b%m



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

def extract_factors(N, relations, null_space,relations2,colist,z_list,disc_list,qchar,quad_base,mod_list):
  #  print(" quad_base: "+str(quad_base))
    n = len(relations)
  #  mod1=19
    counter=0
    for vector in null_space:
       # counter+=1
        prod_left = 1
        prod_right = 1
        disc=1
        bits=[]
        y0sqr=1
        ja=0
        check=0
        tmod=1
        for idx in range(len(relations)):
            bit = vector & 1
            vector = vector >> 1
            bits.append(bit)
          #  print("vector: ",bit)
            if bit == 1:
                if check ==0:
                    counter+=1
                check=1
                #prod_left *= roots[idx]
                prod_left *= relations[idx]
                prod_right *= relations2[idx]
                disc*=disc_list[idx]
                z=z_list[idx]
                tmod*=mod_list[idx]
               # g=gcd(z,mod)
                
              #  y0=colist[idx]%mod1
              #  y0sqr*=y0**2


               # if relations[idx]%N != relations2[idx]%N:
                    #print(relations[idx]%N)
                    #print(relations2[idx]%N)
                    
                print("[Debug Info]Adding to left mod N: "+str(relations[idx])+" adding to right mod N: "+str(relations2[idx])+" y0: "+str(colist[idx])+" z: "+str(z_list[idx])+" discriminant: "+str(disc_list[idx])+" qchar: "+str(qchar[idx])+" mod: "+str(mod_list[idx]))#+" ja:"+str(ja))

            idx += 1
        if check == 0:
            continue
        jacobi_symbols=[]
        for prime in quad_base:
            jacobi_symbols.append(jacobi(disc%prime,prime))

        #print("bits: ",bits)

        sqrt_right = math.isqrt(abs(disc))
        sqrt_left = math.isqrt(prod_left)
        print("sqrt_right: "+str(sqrt_right%N)+" sqr_right: "+str(sqrt_right**2%N))
        print("sqrt_left: "+str(sqrt_left%N)+" sqr_left: "+str(sqrt_left**2%N))
        z=((prod_left-prod_right)//N)//4
        y0=sqrt_left#%mod1

        y1=sqrt_right#%N
        testsqr=isqrt(disc)
        print("[i]Debug info: disc%N: "+str(disc%N)+" (disc): "+str(disc)+" jacobi_symbols: "+str(jacobi_symbols)+" prod_left: "+str(prod_left)+" prod_right: "+str(prod_right)+" z: "+str(z%N)+" ja: "+str(ja)+" y0: "+str(y0)+" tmod: "+str(tmod))#+" x: "+str(x)+" eq1: "+str(eq1_full%mod1))

        if testsqr**2 == disc:
            print("[i]Found a square in the integers, the quadratic character base worked")
      
            if sqrt_right**2%N != sqrt_left**2%N:
                print("[i]Our two squares are not congruent mod N")
            sqrt_left = sqrt_left % N
            sqrt_right = sqrt_right % N
            factor_candidate = gcd(N, abs(sqrt_right-sqrt_left))
            if factor_candidate not in (1, N):
                other_factor = N // factor_candidate
                print("[SUCCESS]Factors are: "+str(factor_candidate)+" and "+str(other_factor))#,counter)

                return factor_candidate, other_factor
        else:
            print("[i]Discriminant is not a square in the integers")
        print("\n")
        time.sleep(10000000)
    return 0, 0


def solve_bits(matrix):
    n=p_amount*2+4+quad_char_amount
    lsmap = {lsb: 1 << lsb for lsb in range(n+100)}
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

def build_matrix(factor_base, smooth_nums, factors,smooth_nums2,factors2,qchar):
    print("\n\n")
    fb_len = len(factor_base)
    fb_map = {val: i for i, val in enumerate(factor_base)}
    ind=1
    #factor_base.insert(0, -1)
    M2=[0]*(p_amount*2+4+quad_char_amount)
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


    offset=p_amount*2+4
    ind=1
    i=0
    while i < len(qchar):
       # print("qchar[i]: ",qchar[i])
        j=0
        while j < len(qchar[i]):
            if qchar[i][j]==-1:
                M2[j+offset] |= ind
            
            j+=1
        ind = ind + ind
        i+=1
  #  for vector in M2:
       # bits=[]
        #for idx in smooth_nums:
       # print("row: ",bin(row))

          #  bit = vector & 1
          #  vector = vector >> 1
          #  bits.append(bit)
      #  print(bits)
    return M2

def QS(n,factor_list,sm,colist,flist,sm2,flist2,qchar,z_list,disc_list,quad_base,mod_list):
    global p_amount
    factor_list.insert(0,2)
    factor_list.insert(0,-1)
    p_amount=len(sm)
    g_max_smooths=p_amount*2+4+quad_char_amount
    if len(sm) > g_max_smooths: 
        del sm[g_max_smooths:]
      #  del xlist[g_max_smooths:]
        del flist[g_max_smooths:] 
        del sm2[g_max_smooths:]
        del flist2[g_max_smooths:]   
    M2 = build_matrix(factor_list, sm, flist,sm2,flist2,qchar)
    null_space=solve_bits(M2)
    f1,f2=extract_factors(n, sm, null_space,sm2,colist,z_list,disc_list,qchar,quad_base,mod_list)
    if f1 != 0:
        #print("[SUCCESS]Factors are: "+str(f1)+" and "+str(f2))
        return f1,f2   
    #print("[FAILURE]No factors found")
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

def find_comb(lists,n,procnum,return_dict,hmap,quad_base,hmap2):
   
    z=g_z
    z_max=z+n
    factor_base=copy.deepcopy(lists)
    factor_base.insert(0,len(factor_base)+1)
    factor_base=array.array('Q',factor_base)
    while z < z_max:
        ret_array=[[],[],[],[],[],[],[],[],[]]
        i=0
        while i < p_mod_amount:
            collected=[]   
            mod=1
            try:
                l=hmap[i][str(z%lists[i])]
               # if l[0]==0:
                collected.append(lists[i])
                collected.append(l)
                mod*=(lists[i]**2)
            except Exception as e:
                i+=1
                continue

            
            if gcd(z,mod)!=1: 
                i+=1
                continue
            if mod ==1:
                i+=1
                continue
        

            collected= lift_all(collected,z,n,2)
            collected=get_partials(mod,collected)
            gen_comb(collected,mod,z,n,factor_base,ret_array,quad_base,hmap2)

            i+=1
        z+=1
   # sys.exit()
  #  test=QS(n,lists,ret_array[0],ret_array[1],ret_array[2],ret_array[6],ret_array[7],ret_array[3],ret_array[5],ret_array[8],quad_base, ret_array[4])   
    return 0
     

def get_primes(start,stop):
    return list(sympy.sieve.primerange(start,stop))

def main():
    global key
    global base
    global workers
    k=0
    while k<1:
        start = default_timer() 
        if g_p !=0 and g_q !=0 and g_enable_custom_factors == 1:
            p=g_p
            q=g_q
            key=p*q
        if key == 0:
            print("\n[*]Generating rsa key with a modulus of +/- size "+str(keysize)+" bits")
            publicKey, privateKey,phi,p,q = generateKey(keysize//2)
            n=p*q
          #  key=n
        else:
            print("[*]Attempting to break modulus: "+str(key))
            n=key

        sys.set_int_max_str_digits(1000000)
        sys.setrecursionlimit(1000000)
        bits=bitlen(n)
        primeslist=[]
        plist=[]
        quad_base=[]
        print("[i]Modulus length: ",bitlen(n))
        print("[i]Gathering prime numbers..")
        primeslist.extend(get_primes(3,1000000))

        while len(plist)!=p_amount:
            if n%primeslist[0] !=0:
                plist.append(primeslist[0])
            primeslist.pop(0)

        while len(quad_base)!=quad_char_amount:
            if n%primeslist[0] !=0:
                quad_base.append(primeslist[0])
            primeslist.pop(0)
       # plist=[3,11]#[5,7,11,13,17,19,23]
        launch(n,plist,quad_base)
        duration = default_timer() - start
        print("\nFactorization in total took: "+str(duration))
        k+=1

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


