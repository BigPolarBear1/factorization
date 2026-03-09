Disclaimer: At no point did AI contribute anything to this research project. Copilot can't even calculate basic congruences without making mistakes. People who think AI can do novel math research are delusional. I did not use AI for the paper, and I did not use AI for any of the code. I attacked this problem, from the very basics 3 years ago, slowly burrowing deeper and deeper. 

#### To run from folder "NFS_WIP" and "NFS_Variant_Simple" (Will be adding more NFS related code in the coming days/weeks):</br></br>

Note: Experimental WORK IN PROGRESS.</br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 50 -base 500 -debug 1 -lin_size 100_000 -quad_size 1</br></br>

You can run the above command from either NFS_Variant_Simple or NFS_WIP.

NFS_Variant_Simple represents chapter Chapter VII in the paper. This is an intermediate step between QS and NFS.
NFS_WIP is my first attempt at porting these findings to a proper NFS algorithm. The downside being that we are restricted to sieving a single pair of coefficients. 

Update: I spent a few hours this evening exploring the best approach. That approach described in the last chapter of the paper... using that "k" variable so that we end up with smooths having the same coefficients , is feasible. In fact I just figured out that there exists a very quick calculation to find what this "k" value should be. Hmm. I'll move on this quicly tomorrow.

Update: Good day. My mood is better today. I implemented some code that will calculate new roots and a k value (as in x^2+yx-Nk) such that a pair of smooths with distinct coefficient can be made to share the same coefficient (prerequisite for NFS). This code works for two smooths all the time.. now I need to generalize that same math so it works for an arbitrary number of smooths. Which should be straightforward enough. Once that is done.. I'll upload a PoC and the dominoes should start falling quickly. I don't understand this.. getting ignored like this. It doesnt matter bc once I succeed I will make sure everyone knows, and I'm super close now.

Update: Hmpf, trying to combine more smooths then 2 seems to break my formula. I wonder if I simply need to increase the degree of my polynomials. That would make the most sense actually... I'll figure it out tomorrow. My mood just massively crashed again. I feel absolute anguish about my situation. It is like people are being told to ignore me on purpose and make it impossible for me to build a future or find employment. It literally feels like the pre-war years before WWII.. banning undesirables from public life. 

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

