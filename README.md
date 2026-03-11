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

TO DO: In NFS_WIP2 the main bottleneck is finding a good "k" value in find_common_coefficient(). Right now it is doing it in a really bad bruteforce way (plus it is also not calculating all possible coefficients, since it only takes into account individual primes and not combinations of primes). I have found a formula to calculate this "k" value for 2 smooths, something similar should work for 3 and more smooths as well.. but it is quite complicated.. as it goes back to those first chapters in my paper, the relation between p+q and p*q (where p and q are factors of semiprime N).. there's this link between addition and multiplication.. but it is very complicated.. but I know it will be possible based on what I'm seeing, the problem is simply that this k value will become a big number very quickly, even though the factors of k will remain predictable... so that doesn't end up being an issue. In addition once that formula is completed we will need a lot less smooths to be found in "phase one" (sieving a-Nk = p).. since we should be able to then find a common coefficient for most, if not all, smooths found.
Then once that is implemented.. we just add all the same bells and whistles we have in CUDA_QS_variant.

Update: I got it... so finding a k value that works for many smooths.. it comes from the factorization of the roots (which in itself comes from the factorization of a, as in a-Nk=p). I should be good now to clear that bottleneck and finally bring my work to completion.. just now, I was mentally going some really bad places... it is always in my moments of greatest anguish that I finally see it. 

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

