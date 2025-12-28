Disclaimer: At no point did AI contribute anything to this research project. Copilot can't even calculate basic congruences without making mistakes. People who think AI can do novel math research are delusional.

NOTE: Starting 2026, none of my research will be published. Only people who treat me with respect will be allowed access to my work. And NATO countries/big tech are very low on that list after harassing me for years and treating me like shit. And I garantuee you, I will succeed at finding a polynomial time algorithm. There is no one else alive in this fucking world more determined then me to succeed at this. Fucking losers.


#### To run debug.py" (Prints the linear and quadratic coefficients to solve for 0 in the integers, for use with my paper):</br></br>

To run: python3 debug.py -keysize 12

This basically creates a system of quadratics. Solving them mod p is easy. But there is only one root solution (the factor of N) which solves the system for 0 for any mod p (aka solves it in the integers). Figuring out how to exactly do this quickly is still an ongoing area of research for me. And if a polynomial time algorithm for factorization exists, it is likely done by solving this system of quadratics. Finding a polyomial time algorithm is my ultimate goal, as this would make progress toward solving p = np as well. 

#### To run from folder "CUDA_QS_variant" (WIP):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run:  python3 run_qs.py -keysize 140 -base 1000 -debug 1 -lin_size 10_000_000 -quad_size 1_000</br></br>

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

Update: I did some thinking on how to best approach this. I should use p-adic lifting for each prime, to prepare a sieve interval on disk.. and lift up to a certain bound, determined by how much disk space we have.. Once that is written to disk.. we only need to read it.. which is reasonably fast. Its writing to disk which is what really sucks. After that we can pull them from disk.. arrange them in a 2d matrix (use tiling and slicing to get the correct shape) and call sum() in one go in the GPU. I think that is probably the best route. And once the sieve interval is fully done.. we could opt to temporarily write it to disk and when we move to processing the sieve interval.. just stack a bunch of them ontop when calling putmask() in the gpu.. because the more we do in one go there the better it seems. Like all the stuff we do in the GPU... we need to really do as few calls as possible with as much memory as possible. Anyway.. thats my goal for tomorrow. Then I'll check what still needs improvement. 


I will implement what I described above. I'm not 100% sure yet if I'll just save root info to disk and construct sieve intervals in memory or save entire sieve intervals to disk.. I'm going to do some numbers today... like if we have a bound for quadratic coefficients of say 1 to 100_000, and a factor base of 10_000, then what is the worst case scenario for disk space requirements if we use entire sieve intervals... that I need to figure out the exact math for first. Then re-implement p-adic lifting and do the rest like I said above... just have three major calls in the gpu sum(), putmask() and non_zero() for creating and processing sieve intervals... where we want to process as much memory in bulk as we can because thats where GPU shines. Once that is all good I'll finish the paper. Make a post on reddit or some shit about my work.. and stop publishing future research. There is plenty of people interested in my work.. and strangely enough they are all from non-nato countries. It is hilarious, the moment I began working on factorization I became a complete outcast in the west, fucking tech bro losers jealous bc they couldnt do this kind of shit in a 100 years. I really detest western infosec. I was never able to have a real career in the west, the only person to support me and believe in me got fired for supporting and believing in me... all these losers at msrc and the microsoft red teams, they always hated me because I was finding the bugs they didnt. Western infosec is in a pitiful state. They will lose. I'm going to make the people who respect me, the most powerful cyber force in the world, and I dont care if that means I will need to leave the west and go live elsewhere. There is nothing left for me here.

And yea, if one day the western police tries to arrest me, I fully intend to go out in a blaze of glory. I'm not ever seeing the inside of a jail cell alive. Either way, its just a matter of time until I permanently move to a country where they cant get me anyway.
