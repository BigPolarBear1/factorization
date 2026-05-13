Disclaimer: At no point did AI contribute anything to this research project. Copilot can't even calculate basic congruences without making mistakes. People who think AI can do novel math research are delusional. I did not use AI for the paper, and I did not use AI for any of the code. I attacked this problem, from the very basics 3 years ago, slowly burrowing deeper and deeper. </br>
Disclaimer2: I never finished highschool. I started working on this 3 years ago, not even knowing basic algebra (I can explain exactly how I can to all the conclusions I came to if someone ever cares to hear about it). This is 3 years. I'm only going to get better. I'm still far away from plateauing skill wise.

The day I break factorization will be day 0 of the gay future. A better world will be born. A world of human creativity. A world of art, science and spirituality. A world without suffering. I will see that world become reality, no matter the cost. 

Been unemployed for years without income, still looking for work: big_polar_bear1@proton.me ... able to relocate to anywhere in the world except the US.

#### References (re-used many of the core number theoretical functions from these PoCs to fit my own algorithm): 
https://stackoverflow.com/questions/79330304/optimizing-sieving-code-in-the-self-initializing-quadratic-sieve-for-pypy
https://github.com/basilegithub/General-number-field-sieve-Python 

#### About the paper
Math paper is a work in progress. The final chapters are a bit rushed and building an algorithm around p-adic lifting isnt as straight forward as I had assumed. I do think there is an angle there I can exploit, but I'll do some further experimentation first and get a working PoC before I make edits to the paper again.

#### To run from folder "binomial_sieve":</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run:  python3 run_qs.py -keysize 75 -base 500 -debug 1 -lin_size 1_000 -quad_size 1_000_000 -t 70

Update: Removed some things from the paper. I did some testing with B-smooths where a small factor is added to be the polynomial value and N multiplier, k. Thinking these might count as unique B-smooths, but they dont and contribute nothing extra to the linear algebra step. Hence the only real option is to try and use our setup to get more factor overlap then what SIQS achieves while ignoring these trivial cases. I'll focus only this now. I.e Find a B-smooth...then use all these different dimensions we can sieve in (binomial expansions,multiples of N, different binomial terms, roots generating non-zero deratives.... etc) to just keep going until we get a hit on nearly identical b-smooths. It's the only way that it can be done. I've exhausted everything else at this point. So this is the final thing I'll try before continueing with my research in private.

Update: Been doing additional thinking. The best strategy will be starting like we do in CUDA_QS_Variant, using square moduli, so that any B-smooth we find will have as few odd exponent factors as possible (while still being able to sieve efficiently). Then we can sieve for those odd exponent factors in all these different dimensions using binomial expansions and everything.  This will be quite key, making sure our b-smooth with all the squares factored out, is as small as possible compared to N. 

Update: Alright, so first thing I'm going to do is optimize find_same() to find the initial b-smooth as fast as possible and with as few odd exponent factors as possible. Once that is done we will shift our focus to find_same2().

Update: Eureka. That was it. Will upload PoC tomorrow :)

Update: Need one more day. My intuition was good.. but the implementations I attempted these last few months were way to overcomplicated. The solutions ended up being fairly trivial in retrospect. I feel like standing on a tall mountain now after 3 years and finally have clarity. Will upload tomorrow. 

You know... today is one of those days, where I feel an important milestone has been reached in my math career. I guess 3 years of hard work, starting from nothing is what it takes to master a niche. Tomorrow I am going to upload something beautiful. Perhaps after tomorrow life will somehow turn out alright... perhaps people will finally see what I see. 

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


