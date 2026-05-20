Disclaimer: At no point did AI contribute anything to this research project. Copilot can't even calculate basic congruences without making mistakes. People who think AI can do novel math research are delusional. I did not use AI for the paper, and I did not use AI for any of the code. I attacked this problem, from the very basics 3 years ago, slowly burrowing deeper and deeper. </br>
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
To run:  python3 run_qs.py -keysize 80 -base 500 -debug 0 -lin_size 100_000 -quad_size 1_00

Update: Did a quick refactor of find_same(). This illustrates how to sieve the polynomial value of the polynomial without -n*k. These residues can actually be saved on disk and re-used since they are not influenced by semiprime N. 
Still to do now... implement support for higher degrees.. and selecting optimal coefficient ranges to sieve... such that the polynomial value of the polynomial with -n\*k is as small as possible.. because there we only want as factors, small factor and the big factors that were found in the original B-smooth (the B-smooth that ends up calling into find_same)
Anyway.. going for a run I think. I'll do some more work soon.

Update: Added some code that will calculate the optimal coefficients. But right now we're only working with degree=2 and calculating the optimal coefficient for the linear term. Next step is to use higher degree and more coefficients.. I'll need to do some thinking how to best sieve that though... because the idea is that by using higher degree polynomials.. we'll have smaller coefficients.. so I'm hoping to use that somehow.

Update: Yeap definitely need to move up in degree.. let me do some testing tomorrow. If it will work like it think it might.. then this will be over shortly.

Update: Did some more 3am pondering. So definitely need to use p-adic lifting and higher degrees here. And just precompute a bunch of stuff and have it sitting on disk. We dont need to calculate all possible coefficients mod p ... just up to a certain range. so its not too bad. And the thing is, with higher degree polynomials... the coefficients we want wont be obscured by the modulus so we can directly query a database. 

Update: Too many AI assisted "math expert" trolls on twitter shitting on my work. Anyway.. time to grind. Fck those losers. The one challenge left now.. is use p-adic lifting and building residue maps for polynomials of arbitrary degree so that we can actually quickly look up if for a given root a solution is found near some polynomial without the constant -Nk which factors completely over the factor base (well only small factors from the factor base, since we dont want to introduce large factors there). That part I'm stll working on.. because right now its just simply shifting the burden from one side of the square relation x^2-Nk=y^2 to the other, but this setup should give use a lot more fine grain control to get specific B-smooths. Just need to make it practical now. Morrons dont have the brains to understand even half of what I'm doing. Literally just harassing me with what is clearly AI generated shit that is completely missing the point.

Update: I'll be AFK from socials, email, everything.. going to focus for a week now on finishing this, at the end of it now. I do see how this allows me to achieve what I've been trying to do for a long time... rather then just shifting the burden of having to do trial factorization for x^2-Nk = y^2 from x to y. These morrons dont have the imagination to see what I see and neither does their AI.

Update: Pushed another updates. Now what I need to do is have a small root. And coefficients with factors in common... bc that way, no matter if we sieve the root, one side of the congruence is garantueed to have those shared factors in the coefficients.. Thats going to be the easiest way to appraoch this. Cant really see another way that this would be practical otherwise. So next step is, just increasing the size of the coefficients to be closer to this optimal coefficients that are going to produce small polynomial values. Then before deciding wether or not to build a sieve interval... just calculate the GCD between coefficient and only proceed is its bigger then some bit length.. then later we can come up with something more elaborate.

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


