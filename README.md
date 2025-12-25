Disclaimer: At no point did AI contribute anything to this research project. Copilot can't even calculate basic congruences without making mistakes. People who think AI can do novel math research are delusional.

NOTE: Starting 2026, none of my research will be published. Only people who treat me with respect will be allowed access to my work. And NATO countries/big tech are very low on that list after harassing me for years and treating me like shit. And I garantuee you, I will succeed at finding a polynomial time algorithm. There is no one else alive in this fucking world more determined then me to succeed at this. Fucking losers.


#### To run debug.py" (Prints the linear and quadratic coefficients to solve for 0 in the integers, for use with my paper):</br></br>

To run: python3 debug.py -keysize 12

This basically creates a system of quadratics. Solving them mod p is easy. But there is only one root solution (the factor of N) which solves the system for 0 for any mod p (aka solves it in the integers). Figuring out how to exactly do this quickly is still an ongoing area of research for me. And if a polynomial time algorithm for factorization exists, it is likely done by solving this system of quadratics. Finding a polyomial time algorithm is my ultimate goal, as this would make progress toward solving p = np as well. 

#### To run from folder "CUDA_QS_variant" (WIP):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 120 -base 500 -debug 1 -lin_size 100_000 -quad_size 100</br></br>

Prerequisites: </br>
-Python (tested on 3.13)</br>
-Numpy (tested on 1.26.2)</br>
-Sympy</br>
-cupy-cuda13x</br>
-cython</br>
-setuptools</br>
(please open an issues here if something doesn't work)</br></br>

Additionally cuda support must be enabled. I did this on wsl2 (easy to setup), since it gets a lot harder to access the GPU on a virtual machine.

To do: I will finish this PoC, and after that any work I do will be kept private as my career is currently going nowhere.
A lot of minor optimizations still need to be done. A big ticket item which I hope to implement tomorrow is building a sieve interval with many quadratic coefficients at once. So we only have to call sum() once for many quadratic coefficient. One easy way would be to have 64bit numbers and have 16bit represent a quadratic coefficient, so we can cram 4 different quadratic coefficients in there. In addition... we can also just add more elements to accomadate more space for more quadratic coefficients.. so it just becomes a matter of how much memory the GPU can handle. After that, the main bottleneck will just be trial factorization. We should probably have two different processes. One creating sieve intervals in the GPU and then a bunch of worker threads on the cpu doing trial factorization.
Oh and ofcourse also in roll_interval2d() ... we shouldnt be making copies of an empty 2d interval... because then we may aswell use a single interval array... we need some fast way to create a 2d interval using slicing or something... and prevent arrays from being created. I'll figure something out.. I guess I may just start by using shift_right tomorrow... see how that works. From reading (i have 0 experience with CUDA) the main bottleneck seems to be transferring to and from the GPU.. so I dont want to create any new arrays.. I just need to reuse the same one by shifting.

Update: I was doing some thinking. And actually, if we use a fairy small lin_size but we have many quadratic coefficient in our 2d sieve interval at the same time... we can use much larger moduli that way. Lets get a proof of concept working today. Goal for today is, have multiple quadratic coefficients represented in the 2d interval. 
