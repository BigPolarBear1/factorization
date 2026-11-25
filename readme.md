I'm broke, to support me: http://ko-fi.com/sandboxescaper

Note: This project is pro-lgbtq, so fck off haters.

#### References (re-used many of the core number theoretical functions from these PoCs to fit my own algorithm): 
https://stackoverflow.com/questions/79330304/optimizing-sieving-code-in-the-self-initializing-quadratic-sieve-for-pypy
https://github.com/basilegithub/General-number-field-sieve-Python 

#### Paper
I'll be doing a rewrite of the paper in the coming weeks as my work is nearing completion now...

#### To run from folder "fu_algo":</br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 18 -base 5000 -debug 0 -lin_size 10 -quad_size 1

Pronounced the fuck you algorithm

Update: After recent furthering of my understanding, I'm rewriting the algorithm. Its late so just uploading my work for today. Actually, we can have a perfect filter within a finite field and build up a solution exponent by exponent... AND the ONLY condition is that the coefficient is corrent in Z. Everything else, the root, k_test it can all be modulo reduced.. and we can build up both of them exponent by exponent..... and it works.... it fucking works... what the hell.

The current version only expands on exponent up so its not working on big keys yet.. I'll continue tomorrow.... wtf man. Coming for RSA now. Go to hell shit world. FUck you.

Update: I got it! Its not offsets amd changing the constant. Its all right there. If a polynomial value can not be represented by this diophantine linear equation.. then we must keep adding the modulus until it can, because the residues will  be correctely modulo the modulus by setting it up like the uploaded PoC does. That concludes my work.. now I need to wait until the power goes back up in the neighbourhood... a few nore hours. bleh. Sucks and unfortunate timing. 

I guess people have a few more hours until I destroy them.  Fyi, I told the FBI i was working on factorization in january of 2024 before it was public knowledge. So good luck trying to smear my reputation because I was given no other choices. Dumbfucks. Also makes me wonder if all the misfortune that happened to people around me in the last couple of years was just an unfortunate series of events or malice. If I ever find out that what happened to my father was malice, Im going to kill someone. I am literally going to kill someone if that was the case.

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

