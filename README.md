Disclaimer: At no point did AI contribute anything to this research project. Copilot can't even calculate basic congruences without making mistakes. People who think AI can do novel math research are delusional. I did not use AI for the paper, and I did not use AI for any of the code. I attacked this problem, from the very basics 3 years ago, slowly burrowing deeper and deeper. </br>
Disclaimer2: I never finished highschool. I started working on this 3 years ago, not even knowing basic algebra (I can explain exactly how I can to all the conclusions I came to if someone ever cares to hear about it). This is 3 years. I'm only going to get better. I'm still far away from plateauing skill wise.

The day I break factorization will be day 0 of the gay future. A better world will be born. A world of human creativity. A world of art, science and spirituality. A world without suffering. I will see that world become reality, no matter the cost. 

#### References (re-used many of the core number theoretical functions from these PoCs to fit my own algorithm): 
https://stackoverflow.com/questions/79330304/optimizing-sieving-code-in-the-self-initializing-quadratic-sieve-for-pypy
https://github.com/basilegithub/General-number-field-sieve-Python 


#### To run from folder "Coefficient_Sieve" (Experimental WORK IN PROGRESS):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 40 -base 40 -debug 1 -lin_size 5000 -quad_size 1000 -d 2

This is an algorithm I had already discovered a year ago (and had uploaded to github), but being a math novice, I dismissed it as it didn't seem practical and deleted it again. 

I have added information about how this works to the newest version of my math paper, found in this github repo. My mistake last spring/summer was that I didn't quite understand the relation to the binomial theorem.. so when I tried to find the correct multiple of N in the discriminant.. I couldn't pull it off, but now to I do understand it.. I'm pretty sure I can do it now. I need a few days to refactor and include this in my PoC. If that works like I believe it should.. then that's it.. my work will be done.

If it still doesn't work, then it is time to start investigating using higher degree polynomials and how we can sieve more effectively with those.

Update: Some pretty bad headache right now. I'll figure it out tomorrow. I'm fairly sure I can use these binomial expansions to eliminate incorrect solutions mod p and build up a solution in the integers that way. I see a way to potentially do it. If at all possible... I'll know tomorrow fur sure. Today was just a wasted day, just couldn't focus.

Update: Making rapid progress now. I'll share more tomorrow. But singular roots in a primefield they give information about the divisors of the binomial terms.. I just had an ephiphany.

Update: Paper is somewhat done. Ergh, its kind of rushed. I dont want to waste to much time on it for now until I fix the code.. either that code is going to work or it wont.. in which case I'll have to change my strategy and end up having to rewrite much of that paper.. so not going to waste more time on it. Tomorrow.. coding begins :).

Update: I spent the day brainstorming this. I'm not happy with the paper. It describes correctly how it relates to the binomial term and all that.. but I think the way to create a powerful algorithm is infact by moving up to higher degree polynomials. But now that I understand those binomial terms I'm not stuck anymore having to calculate the discriminant for higher degrees.. so atleast it wasnt a waste of time. I'll try and update the uploaded version so it works with degree 4, for example.

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


