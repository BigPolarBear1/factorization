Disclaimer: At no point did AI contribute anything to this research project. Copilot can't even calculate basic congruences without making mistakes. People who think AI can do novel math research are delusional. I did not use AI for the paper, and I did not use AI for any of the code. I attacked this problem, from the very basics 3 years ago, slowly burrowing deeper and deeper. 

#### To run from folder "NFS_WIP" and "NFS_Variant_Simple" (Will be adding more NFS related code in the coming days/weeks):</br></br>

Note: Experimental WORK IN PROGRESS.</br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 50 -base 500 -debug 1 -lin_size 100_000 -quad_size 1</br></br>

You can run the above command from either NFS_Variant_Simple or NFS_WIP.
To run from NFS_WIP2 use (super super early draft): python3 run_qs.py -keysize 20 -base 50 -debug 1 -lin_size 200 -quad_size 1

NFS_Variant_Simple represents chapter Chapter VII in the paper. This is an intermediate step between QS and NFS.
NFS_WIP is my first attempt at porting these findings to a proper NFS algorithm. The downside being that we are restricted to sieving a single pair of coefficients. 
NFS_WIP2 is a very early draft of trying to implement Chapter VIII in the paper, where we sieve many different polynomials.

To do: NFS_WIP2 is very ugly and rushed for now, I struggled with depression a lot the last week and did some blogging about bug hunting instead. I just hacked this together really quickly.

1. The main ticket item is during phase two, after sieving... we need to multiply our results with a quadratic coefficient, such that a majority of sieving results will share the same linear coefficient. The uploaded PoC does a bruteforce attempt on this... but its written very poorly. Ideally I'll need to solve systems of linear congruences to figure out a good quadratic coefficient. Thats by far the biggest item to solve.
2. Once 1 is finished, we need to reimplment sieving by moduli, like in NFS_WIP
3. Its discarding all the negative linear coefficient right now.. should add a negative siev interval as well.
4. Just the entire sieving logic needs to be fixed. We need 'a' values such that we know we're going to get linear coeficients within a certain range. I know how that math works, so I need to be a bit smarter about it.

I'll continue uploading as I progress.. and will delete the older versions once this isn't in such a sad state anymore. 


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

