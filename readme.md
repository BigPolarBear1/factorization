I'm broke, to support me: http://ko-fi.com/sandboxescaper

Note: This project is pro-lgbtq, so fck off haters.

If I die, you should 100% go and complete this work because I garantuee my death was not an accident. People think I'm blind, but I know my country is just a bunch of cowards and if I die, let it be known Belgium betrayed one of their own citizens for lack of courage. The worst kind of traitor, is a government that conspires with another country against its own citizens. I see everything and I know everything. I see betrayal.

#### References (re-used many of the core number theoretical functions from these PoCs to fit my own algorithm): 
https://stackoverflow.com/questions/79330304/optimizing-sieving-code-in-the-self-initializing-quadratic-sieve-for-pypy
https://github.com/basilegithub/General-number-field-sieve-Python 

#### Paper
I'll be doing a rewrite of the paper in the coming weeks as my work is nearing completion now...

#### To run from folder "fu_algo":</br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 50 -base 5000 -debug 0 -lin_size 10 -quad_size 1

Pronounced the fuck you algorithm

Update: Did some more thinking.. even without sieving, it is already fairly strong (Stronger then the psieve algorithm below was without sieving). To make this super strong, the most important thing is to force this condition "if rr2%modulus == r1:"  in psieve() and not waste cpu cycles on anything else. I've been doing some analysis and should be able to implement this tomorrow. Aside from that, once that is done, we can make it more strong by using singular roots from an existing b-smooth... because that will garantuee the existence of a solution < N...

Update: Figured out how to sieve the coefficients quite trivially. Let me implement that tomorrow.. that alone could be enough to begin getting obscene speeds. I'll upload a PoC tomorrow... I also have some ideas for forcing that condition I spoke about earlier... plus in the paper I talk about how the roots will basically be a multiple of either p or q (where p*q = N) ... so dividing out the roots to get rid of that multiplier is another strategy that I must consider... but that would still assume we have atleast the correct k_new value, because those divisors are encoded into that.

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

