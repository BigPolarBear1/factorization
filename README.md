Disclaimer: At no point did AI contribute anything to this research project. Copilot can't even calculate basic congruences without making mistakes. People who think AI can do novel math research are delusional. I did not use AI for the paper, and I did not use AI for any of the code. I attacked this problem, from the very basics 3 years ago, slowly burrowing deeper and deeper. 

The day I break factorization will be day 0 of the gay future. A better world will be born. A world of human creativity. A world of art, science and spirituality. A world without suffering. I will see that world become reality, no matter the cost. 

#### To run from folder "NFS_WIP2" (Experimental WORK IN PROGRESS):</br></br>

##NFS related code is borrowed from: https://github.com/basilegithub/General-number-field-sieve-Python (note: Very impressively written, helped me big time, thanks)

To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 50 -base 70 -debug 1 -lin_size 100_000 -quad_size 50</br></br>

Update: Just made some more improvements to use a larger zx+y value. Its kind of trash. Lots of brute force shit which should be done with number theory instead (i.e see create_sieve_interval).. but it doesn't matter.. I'll use this PoC now as a research tool to figure out how I can transform these ideas into using linear algebra to solve for this solution instead.. without having to restrict myself to a single set of coefficients as normal NFS does. That's the main thing now. I have an idea how to do it.. but I'll need to grind out the details now.. final step in the algorithm though.. once that works.. that's it.. it will be done. The uploaded PoC is basically doing fermat's factorization method with legendre symbols.. its kind of shit until I add some linear algebra.

To use linear algebra you would use the quadratic coefficient similar to what I do in CUDA_QS_Variant.. but I have to run the numbers and write the code.. I already know it can be done though.

UPDATE: 
 zx^2+yx-Nk = pv</br>
If the polynomial value, pv, is square then the legendre symbols generated will be the same if for example z = 1 and k = 4 or z = 4 and k = 1 ... as long as the polynomial value is square, these Legendre symbols will be the same.

Hence this proves that we can use these quadratic coefficients (z) to sieve... improving the existing number field sieve algorithm.

I will begin implementing this now... as soon as I can get a functional PoC then this project will be finished.

An intermediate step between QS and NFS (representing chapter VII in the paper) can be found here: https://github.com/BigPolarBear1/factorization/tree/7deba681fc78c349ea70e514a36ab723399f8e96/NFS_Variant_simple

#### To run from folder "CUDA_QS_variant":</br></br>
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

