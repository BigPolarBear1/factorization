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

Update: I struggled a bit trying to go from 2 smooths and using that k-variable to get similar coefficients to 3 smooths.. but almost there now, its definitely possible just using quadratics, its just that these numbers get large really fast which is why it was giving me such a headache.. then from 3 smooths to an arbitrary number of smooths should follow from those results. Bit of a shit day today. Just the absolute agony and hopelessness of my situation. Either people know I am right and closing in fast, and this is their strategy they decided to deal with that (in which case, I'm going to Asia as soon as I can)... or people are oblivious to what is happening.. I dont know which one is worse. Just driving me nuts. I just want to live with some dignity, not like this, I want my own life, my independence.. just wasting the years sitting in a small room at my parents home. If people see me being angry at the west.. this is exactly why I am so angry. And its just the same shit from before 2020. Only my former manager had the courage to do anything about it.. but then microsoft retaliated against him I guess.

Update: I know it is possible what I'm trying to do in a way that would yield a superior version of NFS. So close. This last part is a bit complicated though. 

Update: It is very complicated. There is something going on that ties back to the first chapters of my paper.. with p+q and p*q (where p an q are factors of semiprime N) .. this relation between addition and multiplication... and goes into some really deep number theoretical stuff. I guess what I want to do tomorrow.. let me just write a PoC.. which will represent the final version.. with the only thing missing a function to calculate that proper k value. But lets just bruteforce that k value tomorrow and upload it. Then that is going to get me as close to a finshed PoC as I've ever been with one final challenge left to solve. I know it is possible.. but it is just boggling my mind right now and I dont know how long this final step will take for me to solve.

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

