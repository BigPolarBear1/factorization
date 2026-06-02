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
To run: python3 run_qs.py -keysize 50 -base 5_000 -debug 0 -lin_size 10_000 -quad_size 1

UPDATE: Ignore this version for a few days. Yes, using square moduli for sieving is better then using non-square moduli. Because this always garantuees that a large part of the B-smooth can be ignored as far as linear algebra over Z/2 is concerned. So that gives us a small improvement over SIQS, where SIQS can only achieve something similar if it finds a lot of B-smooth using the same modulus. This in itself could be worth publishing a paper over... I have however also linked all these findings back to NFS in the last few months and if you can introduce square moduli to NFS's approach, reduce the required amount of B-smooths there.. then that jump in performance could give us a real shot at threatening RSA-1024.

As CUDA_QS_variant already shows how to do a SIQS style variant using squaring moduli and all the math behind it, I will transform this one into an NFS variant instead. I should be ready to upload it soon.

Update: Quickly replaced the PoC with a old NFS variant I had modified a few months ago. This needs a lot of fixing and more features. First thing first will be to add support for arbitrary degree.. but in a way that it links back to my own work and not just does what NFS does. Then after that I add lifting for the modulus... optimize everytime and that should be good as a demonstration...

Update: Also got it to work with 4th degree polynomials. I wonder if I just sieve using polynomials generated from the same binomial terms if I cant just combine all of them during linear algebra... its too warm today... going to take a break until later tonight. This heat is giving me brain damage.

Update: Been doing a little bit more work toward implementing higher degree polynomials.. will finish it in the coming days.

Update: Uploaded PoC now uses 4th degree polynomials. I do still need to fix some shit, mainly doing it so that it sieves a proper parabola everytime it calls into sieve3()... especially important when using 4th degree polynomials and larger.. since we can end up generating very big values otherwise.. and ofcourse once that is fixed I should experiment with square moduli. When all of that is done I should see if I can combine sieve results from polynomials of different degrees, but generated with the same binomial term... and see what else can be done there to improve sieving.

Update: So going up in degree.. it keeps that rational side small. The algebraic side is a big problem.. thats exactly why NFS does base-m expansion rather than what I'm doing. However, with the way everything is organized now.. I do have an idea... based on what I was doing prior in the last few days.  

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


