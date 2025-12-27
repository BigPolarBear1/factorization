Disclaimer: At no point did AI contribute anything to this research project. Copilot can't even calculate basic congruences without making mistakes. People who think AI can do novel math research are delusional.

NOTE: Starting 2026, none of my research will be published. Only people who treat me with respect will be allowed access to my work. And NATO countries/big tech are very low on that list after harassing me for years and treating me like shit. And I garantuee you, I will succeed at finding a polynomial time algorithm. There is no one else alive in this fucking world more determined then me to succeed at this. Fucking losers.


#### To run debug.py" (Prints the linear and quadratic coefficients to solve for 0 in the integers, for use with my paper):</br></br>

To run: python3 debug.py -keysize 12

This basically creates a system of quadratics. Solving them mod p is easy. But there is only one root solution (the factor of N) which solves the system for 0 for any mod p (aka solves it in the integers). Figuring out how to exactly do this quickly is still an ongoing area of research for me. And if a polynomial time algorithm for factorization exists, it is likely done by solving this system of quadratics. Finding a polyomial time algorithm is my ultimate goal, as this would make progress toward solving p = np as well. 

#### To run from folder "CUDA_QS_variant" (WIP):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run:  python3 run_qs.py -keysize 140 -base 1000 -debug 1 -lin_size 10_000_000 -quad_size 1_000</br></br>

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

UPDATE: Just uploading my work in progress for the day. Right now all the number theory except p-adic lifting is added to the PoC. I also realized, rather then having a 2d sieve interval of multiple quadratic/linear coefficients in memory, we should keep it on disk. That way can build sieve intervals for all quadratic coefficients at the same time! Because there will be a lot of duplication with residues mod p... so we can really speed things up that way. I'll rewrite that 2d interval to have it all on disk tomorrow...

So goal for tomorrow in build_database2interval() we are going to build sieve intervals for ALL quadratic coefficients that are valid for a modulus in one go and save them on disk. That way we can save on function calls for quadratic coefficients that have the same tiling mod p and really just start using the GPU and vector addition to our advantage. The uploaded code is really bad, once I'm happy with the overall strategy I'm going to need to do a lot of optimization work. I can already see its potential though..

