Disclaimer: At no point did AI contribute anything to this research project. Copilot can't even calculate basic congruences without making mistakes. People who think AI can do novel math research are delusional.

NOTE: Starting 2026, none of my research will be published. Only people who treat me with respect will be allowed access to my work. And NATO countries/big tech are very low on that list after harassing me for years and treating me like shit. And I garantuee you, I will succeed at finding a polynomial time algorithm. There is no one else alive in this fucking world more determined then me to succeed at this. Fucking losers.


#### To run debug.py" (Prints the linear and quadratic coefficients to solve for 0 in the integers, for use with my paper):</br></br>

To run: python3 debug.py -keysize 12

This basically creates a system of quadratics. Solving them mod p is easy. But there is only one root solution (the factor of N) which solves the system for 0 for any mod p (aka solves it in the integers). Figuring out how to exactly do this quickly is still an ongoing area of research for me. And if a polynomial time algorithm for factorization exists, it is likely done by solving this system of quadratics. Finding a polyomial time algorithm is my ultimate goal, as this would make progress toward solving p = np as well. 

#### To run from folder "CUDA_QS_variant" (WIP):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run:  python3 run_qs.py -keysize 140 -base 2000 -debug 1 -lin_size 10_000_000 -quad_size 1_00</br></br>

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

Update: I changed a few things and added p-adic lifting... then wasted 90% of the day tracking down a bug. Anyway... next I will speed up build_database2interval by making use of the fact that sieve intervals for multiple quadratic coefficients would have the same tiling mod p, especially for smaller primes. So we can speed that up. Goal is to really quickly produce many intervals at the same time and then process them in bulk in the gpu using putmask() and non_zero() to quickly find values above a threshold.

So what need to be changed in build_database2interval:

We should only iterate primes_list a single time. Iterate the possible residues mod p and grab all sieve intervals (across all quadratic coefficients) to which this residue applies (this residue is the root mod p minus the root for the modulus.. or something like that). This has to be really fast... I'm still thinking what the most optimal way is..  then just load them into a 2d interval and increase those values in the gpu in one go. Do this for all residues mod p and move to the next prime until we hit the end of primes_list ... at which point we should have succesfully build all sieve intervals for all quadratic coefficients. Then just move to processing those intervals.. and that function we will optimize later. Anyway.. I'll do some brain storming now.. what the best way is to pull sieve intervals that will have the same tiling mod p without building complicated residue hashmap each time we change the modulus... there's probably some fancy vector math I can use.
