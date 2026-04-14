Disclaimer: At no point did AI contribute anything to this research project. Copilot can't even calculate basic congruences without making mistakes. People who think AI can do novel math research are delusional. I did not use AI for the paper, and I did not use AI for any of the code. I attacked this problem, from the very basics 3 years ago, slowly burrowing deeper and deeper. </br>
Disclaimer2: I never finished highschool. I started working on this 3 years ago, not even knowing basic algebra (I can explain exactly how I can to all the conclusions I came to if someone ever cares to hear about it). This is 3 years. I'm only going to get better. I'm still far away from plateauing skill wise.

The day I break factorization will be day 0 of the gay future. A better world will be born. A world of human creativity. A world of art, science and spirituality. A world without suffering. I will see that world become reality, no matter the cost. 

#### References (re-used many of the core number theoretical functions from these PoCs to fit my own algorithm): 
https://stackoverflow.com/questions/79330304/optimizing-sieving-code-in-the-self-initializing-quadratic-sieve-for-pypy
https://github.com/basilegithub/General-number-field-sieve-Python 


#### To run from folder "Coefficient_Sieve" (Experimental WORK IN PROGRESS):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 30 -base 40 -debug 1 -lin_size 5000 -quad_size 100 -t 240

This is an algorithm I had already discovered a year ago (and had uploaded to github), but being a math novice, I dismissed it as it didn't seem practical and deleted it again. 

Update: I have made some changed. The code now uses a proper threshold value to sieve with, this will be needed later one when we add the linear algebra step. I have also added some p-adic lifting code.. but I need to improve it.. there is a bunch of ways to speed up that code. I'll make some improvements later. It will need to be reasonably fast before I implement linear algebra since it will be a core element to get the linear algebra step to work.

There's also a bug I need to resolve.. because that log value should always be the same for squares in Z.. its probably skipping some solutions where one of the coefficients is divisible by Zp.

And also lifting for p = 2 needs to be implemented.. because that one has a lot of singular roots that lift nicely. 

UPDATE: OH YEA.. I should probably only be lifting singular roots... find one that lifts very high with small coefficients. 

UPDATE2: I AM AN IDIOT! I should just lift singular roots to hunt for large squares to reduce the size of smooth candidates.. like what I tried with CUDA_QS_Variant.. but setting it up like this will probably work better. I'll work through the night to produce a PoC.

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


