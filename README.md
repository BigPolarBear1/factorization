Disclaimer: At no point did AI contribute anything to this research project. Copilot can't even calculate basic congruences without making mistakes. People who think AI can do novel math research are delusional. I did not use AI for the paper, and I did not use AI for any of the code. I attacked this problem, from the very basics 3 years ago, slowly burrowing deeper and deeper. 

The day I break factorization will be day 0 of the gay future. A better world will be born. A world of human creativity. A world of art, science and spirituality. A world without suffering. I will see that world become reality, no matter the cost. 

#### To run from folder "NFS_WIP2" and "NFS_Variant_Simple" (Will be adding more NFS related code in the coming days/weeks):</br></br>

##NFS related code is borrowed from: https://github.com/basilegithub/General-number-field-sieve-Python (note: Very impressively written, helped me big time, thanks)

Note: Experimental WORK IN PROGRESS.</br>
To build: python3 setup.py build_ext --inplace</br>
To run NFS_Variant_Simple: python3 run_qs.py -keysize 50 -base 500 -debug 1 -lin_size 100_000 -quad_size 1</br>
To run NFS_WIP2: python3 run_qs.py -keysize 30 -base 500 -debug 1 -lin_size 100_000 -quad_size 1</br></br>

NFS_Variant_Simple represents chapter Chapter VII in the paper. This is an intermediate step between QS and NFS.</br>
NFS_WIP2 in the meanwhile is making progress toward a better NFS variant. Right now it has stripped out the linear algebra step, and simply sieves legendre symbols while keeping what is in NFS the rational and algebraic side square. This needs a lot more work... especially zx+y should be bigger values.. because this relates to the step size (you'll see what I mean). And also we need to use different k values as in zx^2+yx-Nk by multiplying both zx+y, k and the polynomial value with squares. Plus we can also use something similar like generate_modulus in CUDA_QS_variant.. but for Legendre symbols instead.

Now the real endgame will be to implement some type of linear algebra step again to solve these Legendre symbols. I can somewhat see a way to potentially do it, but its fairly complicated.. but I'll get there eventually.

Update: AHA! I got it. So in zx^2+yx-Nk ... actually moving factors of k over to the quadratic coefficient z, changes the Legendre symbols!!!! Which makes sense now that I think of it. Plus I wonder if I cant figure out which factors should be in k and which should be in z using linear algebra. 

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

