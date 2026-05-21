Disclaimer: At no point did AI contribute anything to this research project. Copilot can't even calculate basic congruences without making mistakes. People who think AI can do novel math research are delusional. I did not use AI for the paper, and I did not use AI for any of the code. I attacked this problem, from the very basics 3 years ago, slowly burrowing deeper and deeper. I did experiment from time to time with AI, but it always failed to be useful at math. It is useful for writing simple and short python functions for simple things.. but that is the extend of its usefulness. Saying this here because I am closing in fast now and I know people will claim I just used AI or something bc of my non-math background and being transgender. And I'm honestly tired of the trolls feeding my work into AI and then attacking my work using the garbarge analysis AI makes of it.</br>

Disclaimer2: I never finished highschool. I started working on this 3 years ago, not even knowing basic algebra (I can explain exactly how I can to all the conclusions I came to if someone ever cares to hear about it). This is 3 years. I'm only going to get better. I'm still far away from plateauing skill wise. 

The day I break factorization will be day 0 of the gay future. A better world will be born. A world of human creativity. A world of art, science and spirituality. A world without suffering. I will see that world become reality, no matter the cost. 

Been unemployed for years without income, still looking for work: big_polar_bear1@proton.me ... able to relocate to anywhere in the world except the US.

#### References (re-used many of the core number theoretical functions from these PoCs to fit my own algorithm): 
https://stackoverflow.com/questions/79330304/optimizing-sieving-code-in-the-self-initializing-quadratic-sieve-for-pypy
https://github.com/basilegithub/General-number-field-sieve-Python 

#### About the paper
Math paper is a work in progress. The final chapters are a bit rushed and building an algorithm around p-adic lifting isnt as straight forward as I had assumed. I do think there is an angle there I can exploit, but I'll do some further experimentation first and get a working PoC before I make edits to the paper again.

#### To run from folder "polysieve" WORK IN PROGRES...extremely early version:</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 60 -base 500 -debug 0 -lin_size 100_000 -quad_size 1_00

Did a quick refactor of find_same() ... the goal there wont be fast sieving but to find near identical b-smooths and terminate factorization much sooner.
Next we need to use the residue map for lside (poly without -nk) as a hashmap to quickly test the factorization.. the idea being that we can precalculate residues and p-adically lift them to squares for as many primes as we want. once that's done.. we also need a quick heuristic to keep the polynomial value for the polynomial with -nk as small as possible once divided by the modulus, because we dont want to introduce new large factors there. Once that is done.. just optimize and experiment with higher degree polynomials.

Update: Quickly added a check which compares roots to an optimal root value and only proceed if the bitlength of the difference is smaller then 25 bits. Such that when the polynomial value is divided by the modulus, we end up with a small value that is more likely to only contain small-ish factors. Now the only thing left to solve is the left side... and for that we're going to construct a hashmap with large squares (these large squares dont need to be part of the factor base.. so we can calculate as much as we have disk space.. that's the critical part here.. plus these residues can be re-used for any N.. that's the very critical part here :) ), indexed by roots. But we'll add all the possible roots within a certain range into a single hashmap for all primes.. so that we need to query it only once and immediatly see if a lside value contains any large squares. We want to avoind doing any kind of trial factorization there... just a quick hashmap query. The speed at which it finds smooths in find_same will be very slow.. but once this is implemented, it will yield smooths that have a high chance of succeeding during the linear algebra step. That will be good as an initial proof of concept... then we just come up with ways to gain speed for the algorithm.

#### To run from folder "Coefficient_Sieve" (For use with the paper):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run:  python3 run_qs.py -keysize 40 -base 50 -debug 1 -lin_size 10_000 -quad_size 100</br>

Just demonstrates the math from the paper using quadratics. For educational purposes. And rather then taking a square root over a large prime we can also just calculate the discriminant. But this demonstrates the interesting relation between these quadratics and the factors of N.

#### To run from folder "CUDA_QS_variant" (Failed Experiment):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run: To run:  python3 run_qs.py -keysize 240 -base 100_000 -debug 1 -lin_size 100_000_000 -quad_size 100</br></br>
 
Prerequisites: </br>
-Python (tested on 3.13)</br>
-Numpy (tested on 1.26.2)</br>
-Sympy</br>
-cupy-cuda13x</br>
-cython</br>
-setuptools</br>
-h5py</br>
(please open an issues here if something doesn't work)</br></br>

Additionally cuda support must be enabled. I did this on wsl2 (easy to setup), since it gets a lot harder to access the GPU on a virtual machine.

This was an attempt at finding smooths with similar factorization using an SIQS variant. By using quadratic coefficients. But it didnt end up working as I had hoped so I abondoned this approach, but perhaps someone will get some use out of it.

#### To run debug.py" (Prints the linear and quadratic coefficients to solve for 0 in the integers, for use with my paper):</br></br>

To run: python3 debug.py -keysize 12

This basically creates a system of quadratics. Solving them mod p is easy. But there is only one root solution (the factor of N) which solves the system for 0 for any mod p (aka solves it in the integers). Figuring out how to exactly do this quickly is still an ongoing area of research for me. And if a polynomial time algorithm for factorization exists, it is likely done by solving this system of quadratics. Finding a polyomial time algorithm is my ultimate goal, as this would make progress toward solving p = np as well. 


