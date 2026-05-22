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
To run: python3 run_qs.py -keysize 100 -base 10_000 -debug 0 -lin_size 100 -quad_size 1_00

OMG. The actual solution was super simple and basically what I had been doing already a year ago. ASSYMETRY BETWEEN FACTOR BASE AND REQUIRED NUMBER OF B-SMOOTHS, ACHIEVED! Will optimize now and then "formally" announce my breakthrough (well, as formally as an amateur can).

Update: Slept very poorly last night. Might take a break today and just go running soon. Something we can also do in the PoC, is p-adically lift those roots to any other odd exponent.. but by adding bits to the root, we need to offset that using that k multiplier to N so we keep in that sweet spot where the bits after dividing by the modulus are as few as possible. Plus I'm also seeing some potential other big improvements based on the stuff I was doing in recent weeks... let me do some thinking.

What this PoC does:

1. Uses residue sieving and sieves for b-smooths with a garantueed square modulus as factor.
2. Once a b-smooth is found, bc of the square modulus we can +/- ignore half the factor.. since even exponents can be ignored during guassian elimination over gf(2).
3. Now we jump into find_same(), and we set as modulus, what remains after dividing the b-smooth we found by the square modulus and now we go looking for b-smooths with this modulus.
We use the leading coefficient and multiples of N to keep control over the size of generated b-smooth candidates.
4. Because this results in b-smooth pairings where half the factors can be ignored, we reduce the amount of B-smooths required by half. THIS IS DEMONSTRATED BY THE POC SO I'M NOT GOING TO ARGUE ABOUT THIS WITH IDIOT NERDS USING AI. READ THE F*CKING POC. Dont have patience for idiots.

Now we can reduce the required amount of B-smooths even further (although reducing it by half is already a big milestone/breakthrough).. by making sure in step 3, any polynomial values generated are as small as possible after dividing by the modulus.. meaning fewer new factors are introduced that need to be matched. Which we can do by using polynomials... but I'll demonstrate this later on when everything is optimized.

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


