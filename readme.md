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

Update: After recent furthering of my understanding, I'm rewriting the algorithm. Its late so just uploading my work for today. Actually, we can have a perfect filter within a finite field and can build up a solution exponent by exponent... AND the ONLY condition is that the coefficient is correct in Z. Everything else, the root, k_test (see psieve()) it can all be modulo reduced inside the finite field.. and we can build up both of them exponent by exponent..... and it works.... it fucking works... what the hell.

The current version only expands one exponent up so its not working on big keys yet.. I'll continue tomorrow.... wtf man. Coming for RSA now. Go to hell shit world. FUck you.

Update: Almost there. Let me keep digging. 

Update: Its 1am right now. I was doing some further analysis, looking at different exponents. As a matter of fact, the higher we can lift a root, without the root changing.. the stronger this algorithm will become.... because it won't bump up that k_new value. Headed toward a major calamity and nobody seems to be all the wiser. I knew what I had. I knew it nearly 4 years ago and I still know it today. Just took me too damn long. Let me add this tomorrow.... good thing I'm all powerful and nobody can harm me, good to be born a polar bear.

Update: Added some more updates to the code. Now lets you control minimum bounds for the finite field and maximum bounds for the root. Not entirely convinced yet that it helps, but it gives some more parameters to test with. Sieving the "while u < 10:" (in psieve()) loop is quite trivially done, but that loop should be kept small. The more important loop that needs to get sieved is the coefficient loop..."while i < 100_000:" (in psieve())... either sieved or straight up calculated modulo factor base primes.. I'll explore both.

The most important condition that must be met, and the largest indicator of a root that will solve in Z is this line: "if rr2%modulus == r1:" in psieve() .... so any sieving/calculations to find the right coefficient should leveraege this condition.

Update: Did some more thinking.. even without sieving, it is already fairly strong. To make this super strong, the most important thing is to force this condition "if rr2%modulus == r1:"  in psieve() and not waste cpu cycles on anything else. I've been doing some analysis and should be able to implement this tomorrow. Aside from that, once that is done, we can make it more strong by using singular roots from an existing b-smooth... because that will garantuee the existence of a solution < N...

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

