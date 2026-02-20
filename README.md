Disclaimer: At no point did AI contribute anything to this research project. Copilot can't even calculate basic congruences without making mistakes. People who think AI can do novel math research are delusional. I did not use AI for the paper, and I did not use AI for any of the code. I attacked this problem, from the very basics 3 years ago, slowly burrowing deeper and deeper. 

#### To run from folder "NFS_Variant_Simple" (Will be adding more NFS related code in the coming days/weeks):</br></br>

Note: Experimental WORK IN PROGRESS.</br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 50 -base 300 -debug 1 -lin_size 100_000 -quad_size 10</br></br>

NFS related code is borrowed from (Main sieve(), just to demonstrate we can sieve the same way): https://github.com/basilegithub/General-number-field-sieve-Python (note: Very impressively written, helped me big time, thanks)

I was trying to port some of my insights to the number field sieve algorithm. NFS_Variant_Simple is an intermediate step between QS and NFS. I suspect it may be possible to use the number field sieve approach while sieving many different polynomials rather then just a single polynomial and its multiples.

But it needs more work. I'm very close to figuring it all out. But I cant do it anymore, living like this. Everyday I am tormented by memories of better days. Rollerskating with my former teamlead around the seawall in Vancouver, best days of my life. Places and people I can now not see anymore. Distant memories of better times.. I am stuck. I cant find employment. And its looking like there is only one way out of this now. 

I dont know what people would say. Late 2023, it was the worst time of my life. I was harassed and threatened with a gun (unprovoked). Then microsoft fired me also. I had a utter and complete breakdown as I desperately tried to find a way to stay near my friends. I dont know what microsoft will say to smear my character, but its easy to destroy someone's life and then use the fall-out retroactively to justify destroying someone's life. My life has never been the same since. its easy to point fingers when you dont have context. Should talk to people who actually knew me, like my former managers. Bye.

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

