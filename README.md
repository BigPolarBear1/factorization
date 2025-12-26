Disclaimer: At no point did AI contribute anything to this research project. Copilot can't even calculate basic congruences without making mistakes. People who think AI can do novel math research are delusional.

NOTE: Starting 2026, none of my research will be published. Only people who treat me with respect will be allowed access to my work. And NATO countries/big tech are very low on that list after harassing me for years and treating me like shit. And I garantuee you, I will succeed at finding a polynomial time algorithm. There is no one else alive in this fucking world more determined then me to succeed at this. Fucking losers.


#### To run debug.py" (Prints the linear and quadratic coefficients to solve for 0 in the integers, for use with my paper):</br></br>

To run: python3 debug.py -keysize 12

This basically creates a system of quadratics. Solving them mod p is easy. But there is only one root solution (the factor of N) which solves the system for 0 for any mod p (aka solves it in the integers). Figuring out how to exactly do this quickly is still an ongoing area of research for me. And if a polynomial time algorithm for factorization exists, it is likely done by solving this system of quadratics. Finding a polyomial time algorithm is my ultimate goal, as this would make progress toward solving p = np as well. 

#### To run from folder "CUDA_QS_variant" (WIP):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run:  python3 run_qs.py -keysize 80 -base 500 -debug 1 -lin_size 1_000_000 -quad_size 1_000</br></br>

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

To do: Use larger moduli, and save them to disk. So we don't need to do as much tiling. Additionally we can move workload to the GPU. We should also add support for p-adic lifting to the PoC</br></br>
To do: Paralell reduction in the gpu on the completed sieve interval to find large values. Basically like performing a type of binary search for large values....</br></br>
To do: Chunk the sieve interval .. like first do 0 to 10_000_000 then do 10_000_000 to 20_000_000 and so on.. because if we use a very small step size for the interval then we can keep going for much longer without increasing the smooth value size by much. Additionally I should also reimplement support to use moduli for the step size</br></br>

Ok. Headed into the right direction now. I'll go for a run now. When I come back, I'm going to rewrite chapter 8 in the paper one final time.. then tomorrow I'll start going over that to-do list.. which shouldn't take long. Tomorrow I will add composite moduli support for the database.. this is mostly important for the small primes.. because we don't want to tile mod 3 or some small prime over a huge interval... better to just use a larger composite modulus.
