Disclaimer: At no point did AI contribute anything to this research project. Copilot can't even calculate basic congruences without making mistakes. People who think AI can do novel math research are delusional. I did not use AI for the paper, and I did not use AI for any of the code. I attacked this problem, from the very basics 3 years ago, slowly burrowing deeper and deeper. </br>
Disclaimer2: I never finished highschool. I started working on this 3 years ago, not even knowing basic algebra (I can explain exactly how I can to all the conclusions I came to if someone ever cares to hear about it). This is 3 years. I'm only going to get better. I'm still far away from plateauing skill wise.

The day I break factorization will be day 0 of the gay future. A better world will be born. A world of human creativity. A world of art, science and spirituality. A world without suffering. I will see that world become reality, no matter the cost. 

#### References (re-used many of the core number theoretical functions from these PoCs to fit my own algorithm): 
https://stackoverflow.com/questions/79330304/optimizing-sieving-code-in-the-self-initializing-quadratic-sieve-for-pypy
https://github.com/basilegithub/General-number-field-sieve-Python 

#### About the paper
Math paper is a work in progress. The final chapters are a bit rushed and building an algorithm around p-adic lifting isnt as straight forward as I had assumed. I do think there is an angle there I can exploit, but I'll do some further experimentation first and get a working PoC before I make edits to the paper again.

#### To run from folder "binomial_sieve":</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 75 -base 500 -debug 1 -lin_size 100_000 -quad_size 10 -t 60

This is how this is supposed to be sieved! Finally got it. Let me begin improving this now :).
This differs from quadratic sieve in that we are able to sieve with polynomials with non-zero derivatives, hence giving us more precision (adds a linear offset).

To do:

-Ive already demonstrated how polynomials of degree bigger then 2 can also be used. But I need to re-implement this into the PoC. This is quite an important feature thats still missing.</br>
-Implement p-adic lifting.</br>
-Combine the above two items to hunt for similar factorizations such that success at the linear algebra step becomes independent from the factor base size, which happens to be the largest bottleneck in QS variants.</br>
-Bunch of optimizations.. including size of the binomials we sieve with.. and we can also sieve at different multiples of N using the k variable or the leading coefficient (i'll need to think if there is ever any point in using that leading coefficient).

Update: Going to take a break for today. Technically this code will also work with degrees > 2 now. But tomorrow I will write a function, if we find a B-smooth at the second degree, we then call into a new function thats going to calcualte all the roots modulo the B-smooth for degrees > 2 and different multiples of N (the constant in our polynomials) and try to find B-smooths with similar factorization. It's going to work pretty smoothly, I've got it figured out now.

Update: Ok ok, work out the details for next steps. Tomorrow I'll implement a proper find_same() function. If it is going to work like I hope, then we should be able to succeed with just a handful of B-smooths. Big day tomorrow.... lets see. Will definitely need p-adic lifting for it.

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


