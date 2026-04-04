Disclaimer: At no point did AI contribute anything to this research project. Copilot can't even calculate basic congruences without making mistakes. People who think AI can do novel math research are delusional. I did not use AI for the paper, and I did not use AI for any of the code. I attacked this problem, from the very basics 3 years ago, slowly burrowing deeper and deeper. </br>
Disclaimer2: I never finished highschool. I started working on this 3 years ago, not even knowing basic algebra (I can explain exactly how I can to all the conclusions I came to if someone ever cares to hear about it). This is 3 years. I'm only going to get better. I'm still far away from plateauing skill wise.

The day I break factorization will be day 0 of the gay future. A better world will be born. A world of human creativity. A world of art, science and spirituality. A world without suffering. I will see that world become reality, no matter the cost. 

#### To run from folder "Higher_Degree_QS" (Experimental WORK IN PROGRESS):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 30 -base 500 -debug 1 -lin_size 100_000 -quad_size 1 -d 4

Update: I just uploaded a very barebones version. You can use -d to change the degree used to sieve with. Now we improve this. 
If we sieve a-N = pv, where a is substituted by a polynomial of arbitrary degree.. then based on the value of the coefficients, we know if a will be square or not. And the coefficients can also be used to control the factorization of a. 

Its super simple actually. I don't quite understand why I'm only seeing this now. I guess I just needed to advance enough with my math skills. Either way...  now the fun begins, because understanding all of this, how this works.. we can actually fix CUDA_QS_variant and properly implement the strategy of hunting for smooths with similar factorization. I know how to do it now.. finally seeing the big picture. Hold on.. I'll start improving this quickly now.

#### (abandoning this approach in favor for the one above) To run from folder "NFS_WIP2" (Experimental WORK IN PROGRESS):</br></br>

##NFS related code is borrowed from: https://github.com/basilegithub/General-number-field-sieve-Python (note: Very impressively written, helped me big time, thanks)

To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 50 -base 70 -debug 1 -lin_size 100_000 -quad_size 50</br></br>

An incomplete attempt at merging my findings with NFS. I will delete this once the version above matures more.

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

