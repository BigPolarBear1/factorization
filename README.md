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

UPDATE: OMG. I was out running and I suddenly was hit with a realization. You know, my brain has this weird thing, where I'll wildly overcomplicate things in my head. Its great for doing creative research because sometimes it leads to unexpected insights. Anyway.. I was overcomplication stuff wildly again. LOL. Its actually very simple. We just need to create small sieve intervals mod m. Write them to disk. Then once we calculate many of them.. we just sum() them together onto one very large interval in the GPU. Easy! So easy! Plus we can use p-adic lifting easily aswell this way. Its so easy and elegant and somehow I didn't see it until I was out running in the cold. I might be mentally challenged for real. It was infront of me for so long.... prepare to get a proper PoC before new year... going to start 2026 the right way! 

Lol. Its so simple. Calculate sieve intervals mod m, write them to disk. And then just sum them all together in the gpu onto one really big big interval. I may be writing "lol" but internally I'm screaming.. because this was my one chance to get justice for what microsoft did to the one person who supported and believed in me... and I've wasted to much time. Sound the war drums, we march on my enemies now. Hahahahahaha. Prepare yourselves. I won. You lost. Game over.
