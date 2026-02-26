Disclaimer: At no point did AI contribute anything to this research project. Copilot can't even calculate basic congruences without making mistakes. People who think AI can do novel math research are delusional. I did not use AI for the paper, and I did not use AI for any of the code. I attacked this problem, from the very basics 3 years ago, slowly burrowing deeper and deeper. 

#### To run from folder "NFS_WIP" and "NFS_Variant_Simple" (Will be adding more NFS related code in the coming days/weeks):</br></br>

Note: Experimental WORK IN PROGRESS.</br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 50 -base 500 -debug 1 -lin_size 100_000 -quad_size 1</br></br>

You can run the above command from either two. NFS_Variant_Simple is an intermediate step between QS and NFS and NFS_WIP is merging all my findings into a proper NFS algorithm...

!!!NFS related code is borrowed from: https://github.com/basilegithub/General-number-field-sieve-Python (note: Very impressively written, helped me big time, thanks)

To do for NFS_WIP:

Update: I spent all day trying to figure out how I'm going to sieve arbitrary quadratic in NFS_WIP, similar to how I'm doing it in NFS_Variant_Simple. I think I've mostly figured it out now.

So the first phase we will just sieve quadratics and find cases where the polynomial value factors over the factor base.
Once we have enough of them... its just a matter of multiplying and dividing roots and coefficients until we find common coefficients that are shared among enough of these smooths we found. We should be able to use that precalculated hashmap for that in a way that makes algorithmically sense.

Let me begin writing some code tomorrow that attemps to do this.

Yea that should work... i.e if we have zx^2+yx-4387 = 9, then it doens't matter what coefficients and roots we are using for zx^2+yx, we can always change them using the hashmap. So as long as we find cases like that where the polynomial value factors.. then finding among that similar coefficients should be trivial. Ofcourse we also need to take into account the factorization of zx+y but that should be manage-able. I like this approach. Feels right somehow. My head is going bad places lately, those days I lived in Vancouver and all the friends I knew there.. those days are never coming back. I wish I could go back in time and just relive those years over and over again. It all seems like a distant memory now. 

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

