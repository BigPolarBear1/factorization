I'm broke, to support me: http://ko-fi.com/sandboxescaper

Note: This project is pro-lgbtq, so fck off haters.

#### References (re-used many of the core number theoretical functions from these PoCs to fit my own algorithm): 
https://stackoverflow.com/questions/79330304/optimizing-sieving-code-in-the-self-initializing-quadratic-sieve-for-pypy
https://github.com/basilegithub/General-number-field-sieve-Python 

#### Paper
I'll be doing a rewrite of the paper in the coming weeks as my work is nearing completion now...

#### To run from folder "fu_algo":</br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 35 -base 5000 -debug 0 -lin_size 10 -quad_size 1

Pronounced the fuck you algorithm

Update: A little closer. The thing In noticed is that if k_test is a valid residue modulo the modulus but not in Z, then we can still get the correct root residues.. but when inspecting residues in primes outside the modulus the coefficient will keep growing.. so ideally we shouldn't consider the entire factor base in that case.. maybe just find a bunch of primes that bump up the coefficient inside an interval, and calculate the root using only those primes (and the original root calculated inside the modulus)... anyway, time to go for a run... i'll implement a sieve tomorrow.. until then there isn't really much filtering yet... but I suspect one the interval works properly it will start to take off for real.

Update:

So thing is this, lets say we have a polynomial that solves for 0 in Z:

(41\*31)^2+1832\*(41\*31)-4387\*899 = 0

But that k value (899) can be modulo reduced modulo 729, since the discriminant 1832^2+4\*4387\*899 divides by 729. Which gives us the polynomial: 

(41\*31)^2+14954\*(41\*31)-4387\*170 = 19876185 (note: I use another linear coefficient here, i.e 14954.. as long as its the same mod 729 as 1832... which this one is)

So this root will also have a factor of N. But the polynomial value is only 0 modulo some primes (i.e 5,7).
However if we where to calculate the root modulo 729\*5\*7 then we would succesfully factor N after taking the gcd.
So that gives me a clue as to how this can be solved now. But time for that tomorrow...

Update: Been investigating this deeper. So if the polynomial value is 19876185 (see above example) then this boils down to finding an integer solution to 19876185 = 4387\*729\*x+1271\*729\*y if the integer solution exists, then this polynomial value must also have a factor of N. This is however kind of a chicken and egg problem... however.. I do know from exactly a year ago, how to add an offset to that constant 4387 in my polynomials... AND THEN... the problem becomes much easier to solve. Let me do some digging how exactly I did those things last year.. was around christmas.. should be easy to find those files again. Almost there.. pfft. 



#### To run from folder "psieve":</br>
To build: python3 setup.py build_ext --inplace</br>
To run:  python3 run_qs.py -keysize 150 -base 10_000 -lin_size 10_000 -quad_size 1</br>

note: Up to 120-bit is near instant now.

Some experimentation to optimize psieve_old with claude. Furthermore this is by far the fastest algorithm in the SQUFOF complexity class. But I suspect there is an even faster algorithm hiding and that is what the newer fu_algo is working toward.

#### To run from folder "Coefficient_Sieve" (For use with the paper):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run:  python3 run_qs.py -keysize 40 -base 50 -debug 1 -lin_size 10_000 -quad_size 100</br>

Just demonstrates the math from the paper using quadratics. For educational purposes. And rather then taking a square root over a large prime we can also just calculate the discriminant. But this demonstrates the interesting relation between these quadratics and the factors of N.

#### To run debug.py" (Prints the linear and quadratic coefficients to solve for 0 in the integers, for use with my paper):</br></br>

To run: python3 debug.py -keysize 12

This basically creates a system of quadratics. Solving them mod p is easy. But there is only one root solution (the factor of N) which solves the system for 0 for any mod p (aka solves it in the integers). Figuring out how to exactly do this quickly is still an ongoing area of research for me. And if a polynomial time algorithm for factorization exists, it is likely done by solving this system of quadratics. 

(PERFORMING HEATHEN CYBER RITUALS IN THE BELGIAN MOORLANDS TO DESTROY MY ENEMIES)

