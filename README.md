Disclaimer: At no point did AI contribute anything to this research project. Copilot can't even calculate basic congruences without making mistakes. People who think AI can do novel math research are delusional. I did not use AI for the paper, and I did not use AI for any of the code. I attacked this problem, from the very basics 3 years ago, slowly burrowing deeper and deeper. 

#### To run from folder "NFS_WIP" and "NFS_Variant_Simple" (Will be adding more NFS related code in the coming days/weeks):</br></br>

##NFS related code is borrowed from: https://github.com/basilegithub/General-number-field-sieve-Python (note: Very impressively written, helped me big time, thanks)

Note: Experimental WORK IN PROGRESS.</br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 50 -base 500 -debug 1 -lin_size 100_000 -quad_size 1</br></br>

You can run the above command from either NFS_Variant_Simple or NFS_WIP.

NFS_Variant_Simple represents chapter Chapter VII in the paper. This is an intermediate step between QS and NFS.
NFS_WIP is my first attempt at porting these findings to a proper NFS algorithm. The downside being that we are restricted to sieving a single pair of coefficients. 

To do: Ignore the final chapter in the paper. Upon doing some thinking... it is actually very trivial to ensure both the polynomial value and zx+y remain small even if we are restricted to using a single linear coefficient. Something which I intially figured wouldn't be possible without going to higher degree polynomials.. but I was wrong in my assumptions. We just need x to be a big root and zx+y to be a small root (since zx+y represents another root, for the same quadratic). Not that hard... bah... people would have known. Which makes everything much worse bc this is how I'm being treated, despite having done this..

Here is a most trivial example:
If N = 4387

We could have 4388-4387 = 1

Now 4388 can be rewritten as: 

1097^2-1093\*1097-4387 = 1  and zx+y = 1097-1093 = 4.

Only zx+y and the polynomial value must factor to be consider a valid smooth for the number field sieve algorithm. And here we see an example of achieving a small value for both. Hence, proving that it is possible. 

Update: Did some more thinking. Using the above, it would be trivial to find examples where the polynomial value is square and zx+y is square. And actually we then only need squaredness for the quadratic character base.. but I wonder if I cant simply use the precomputed hashmap for that. Hmm. I'm starting to see something... I can probably skip the entire linear algebra step this way. 


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

