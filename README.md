One final time before I do what I must, looking for work, big_polar_bear1@proton.me. Cant work for big tech bc of what happened in 2023. Cant travel to the US anymore bc I would get arrested (hence cant work for US companies either). Cant say I didnt try when all the infosec shitheads start pointing fingers again.

Disclaimer: At no point did AI contribute anything to this research project. Copilot can't even calculate basic congruences without making mistakes. People who think AI can do novel math research are delusional.

NOTE: Starting 2026, none of my research will be published. Only people who treat me with respect will be allowed access to my work. And NATO countries/big tech are very low on that list after harassing me for years and treating me like shit. And I garantuee you, I will succeed at finding a polynomial time algorithm. There is no one else alive in this fucking world more determined then me to succeed at this. Fucking losers.


#### To run debug.py" (Prints the linear and quadratic coefficients to solve for 0 in the integers, for use with my paper):</br></br>

To run: python3 debug.py -keysize 12

This basically creates a system of quadratics. Solving them mod p is easy. But there is only one root solution (the factor of N) which solves the system for 0 for any mod p (aka solves it in the integers). Figuring out how to exactly do this quickly is still an ongoing area of research for me. And if a polynomial time algorithm for factorization exists, it is likely done by solving this system of quadratics. Finding a polyomial time algorithm is my ultimate goal, as this would make progress toward solving p = np as well. 

#### To run from folder "CUDA_QS_variant" (WIP):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run:  python3 run_qs.py -keysize 160 -base 2000 -debug 1 -lin_size 10_000_000 -quad_size 100</br></br>

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

UPDATE: Changed some more stuff. So right now, we just precalculate a 2d root map that we keep in memory. Its really fast now with all the number theory I've figured out. Plus it can be further optimized since we only really need to calculate roots where quad < prime (and then just repeat that value across the 2d interval in the height)... but I'll change that later as its not really bottlenecking anything yet.
P-adic lifting is broken again. But Ill fix that after I've done everything else. 
So all linear coefficients for a modulus, with the same quadratic coefficient, they actually share the same tiling mod p.. just at different offsets... which means we can speed that up with vector addition in the GPU (still need to do this). And we can even do multiple quadratic coefficients in bulk.. but it would depend on how many of those we can keep in memory at one time.. since we must avoid reading and writing to disk to many times... since that is hugely time consuming. Anyway... I will begin optimizing building the sieving intevals now .. so it builds them in bulk utilizing the GPU. It took me a few days of organizing and re-organizing everything to settle on how exactly I'm going to appraoch this.. but I think I got it down now.

UPDATE: I do wonder. So I know for sure I can bulk build all linear coefficient for a modulus for the same quadatic coefficient. Since they have the same tiling mod p (meaning the distance between root_dist1 and root_dist2 stays the same) ... however, I want to bulk build much more then that to really have this appraoch take off. So we could either bulk build multiple quadratic coefficients..... or actually... multiple moduli..which may be the best since using different moduli doesnt change the tiling mod p. Let me go for a run today and think about how Im going to code this exactly.
