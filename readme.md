I'm broke, to support me: http://ko-fi.com/sandboxescaper

Note: This project is pro-lgbtq, so fck off haters.

#### References (re-used many of the core number theoretical functions from these PoCs to fit my own algorithm): 
https://stackoverflow.com/questions/79330304/optimizing-sieving-code-in-the-self-initializing-quadratic-sieve-for-pypy
https://github.com/basilegithub/General-number-field-sieve-Python 

#### Paper
I'll be doing a rewrite of the paper in the coming weeks as my work is nearing completion now...

#### To run from folder "fu_algo":</br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 20 -base 5000 -debug 0 -lin_size 10 -quad_size 1

Pronounced the fuck you algorithm

Update: After recent furthering of my understanding, I'm rewriting the algorithm. Its late so just uploading my work for today. Actually, we can have a perfect filter within a finite field and can build up a solution exponent by exponent... AND the ONLY condition is that the coefficient is correct in Z. Everything else, the root, k_test (see psieve()) it can all be modulo reduced inside the finite field.. and we can build up both of them exponent by exponent..... and it works.... it fucking works... what the hell.

The current version only expands one exponent up so its not working on big keys yet.. I'll continue tomorrow.... wtf man. Coming for RSA now. Go to hell shit world. FUck you.

Update: Digging deeper into what I found.. so basically we inspect polynomial values and if atleast the coefficient is correct in Z, then figuring out how much we need to add the modulus to the root and k_test is trivial. But I need to generalize all of this a bit more... its still very buggy and inefficient. More in the right direction then what I was doing a week ago though.. really zero-ing in on the polynomial values of the quadratics representing the discriminants.. and using the fact that polynomial values can also be represented as a linear diophantine equation. Which is basically what I'm exploiting in this new code.. but much of it can likely be simplified enormously... which is what I'll focus on next.... in addition, one can use an iterative approach and keep growing the modulus like this I suspect. Let me see... very close now... 

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

