Disclaimer: At no point did AI contribute anything to this research project. Copilot can't even calculate basic congruences without making mistakes. People who think AI can do novel math research are delusional. I did not use AI for the paper, and I did not use AI for any of the code. I attacked this problem, from the very basics 3 years ago, slowly burrowing deeper and deeper. </br>
Disclaimer2: I never finished highschool. I started working on this 3 years ago, not even knowing basic algebra (I can explain exactly how I can to all the conclusions I came to if someone ever cares to hear about it). This is 3 years. I'm only going to get better. I'm still far away from plateauing skill wise.

The day I break factorization will be day 0 of the gay future. A better world will be born. A world of human creativity. A world of art, science and spirituality. A world without suffering. I will see that world become reality, no matter the cost. 

#### References (re-used many of the core number theoretical functions from these PoCs to fit my own algorithm): 
https://stackoverflow.com/questions/79330304/optimizing-sieving-code-in-the-self-initializing-quadratic-sieve-for-pypy
https://github.com/basilegithub/General-number-field-sieve-Python 


#### To run from folder "Coefficient_Sieve" (Experimental WORK IN PROGRESS):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 30 -base 20 -debug 1 -lin_size 10_000 -quad_size 2 -d 2

Update: I very quickly restored a strategy I had tried last summer. This upload only works with degree 2 (use -d 2). And in solve_roots .. we dont need to solve roots, we just need to check if a root exists (i.e using legendre symbols) .. but it doesn't matter. NOW for the next version.. I will add support for higher degrees.. since this will reduce the size of coefficients.. and make this actually work for bigger numbers... and that's it.. that is all I had missed last summer....... and everyone must have known..... I hope you all burn in hell for this.

I'm nearly certain that I can get it to work now by using higher degree polynomials. And in the off-chance that it wont.. i'll just try my original idea of using how this relates to binomial expansion to do targeted hunting for B-smooth numbers. Either way, I am about to succeed. 

Update: I quickly fixed solve_roots() for arbitrary degree. I'll have to simplify this later, because we just need to check the existence of a root rather then calculate the root.. but this works for now. Next I need to figure out how to calculate the discriminant for higher degrees.. since the code will still only work with degree 2 until I fix that.

Update: I guess rather then calculating the discriminant for bigger degrees you would calculate the roots modulo p and then use chinese remainder. Euhm, let me go ahead and optimize that PoC for quadratics first and write it in such a way that it can be expanded to higher degrees later. Now that I see how it works, I should be able to get a decent performance even with quadratics.

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

