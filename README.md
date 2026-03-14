Disclaimer: At no point did AI contribute anything to this research project. Copilot can't even calculate basic congruences without making mistakes. People who think AI can do novel math research are delusional. I did not use AI for the paper, and I did not use AI for any of the code. I attacked this problem, from the very basics 3 years ago, slowly burrowing deeper and deeper. 

#### To run from folder "NFS_WIP", "NFS_WIP2" and "NFS_Variant_Simple" (Will be adding more NFS related code in the coming days/weeks):</br></br>

##NFS related code is borrowed from: https://github.com/basilegithub/General-number-field-sieve-Python (note: Very impressively written, helped me big time, thanks)

Note: Experimental WORK IN PROGRESS.</br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 50 -base 500 -debug 1 -lin_size 100_000 -quad_size 1</br></br>

You can run the above command from either NFS_Variant_Simple or NFS_WIP.

NFS_Variant_Simple represents chapter Chapter VII in the paper. This is an intermediate step between QS and NFS.
NFS_WIP is my first attempt at porting these findings to a proper NFS algorithm. The downside being that we are restricted to sieving a single pair of coefficients. 

NFS_WIP2 is a first version trying to leverage the final chapter. You can run this using (should find a solution most of the time): 
python3 run_qs.py -keysize 20 -base 15 -debug 1 -lin_size 200 -quad_size 1

UPDATE: Ignore NFS_WIP2 and the final chapter in the paper for now. I just realized something. It is very easy to improve NFS_WIP (not NFS_WIP2) massively.. Its not a problem to keep both the polynomial value and zx+y small... its actually quite easy. I dont know why I overcomplicate things so much lol. Let me first thing tomorrow just improve that previous PoC a bit... I understand all the moving parts a little better now... 

Update: Deleted NFS_WIP2. Going to improve NFS_WIP instead... I just realized it is actually trivial to keep both zx+y and the polynomial value small. Hence that solves all my issues and there also wouldn't be a reason to use higher degree polynomial.

#### To run from folder "CUDA_QS_variant":</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run:  python3 run_qs.py -keysize 240 -base 100_000 -debug 1 -lin_size 100_000_000 -quad_size 100</br></br>
 
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

