One final time before I do what I must, looking for work, big_polar_bear1@proton.me. Cant work for big tech bc of what happened in 2023. Cant travel to the US anymore bc I would get arrested (hence cant work for US companies either). Cant say I didnt try when all the infosec shitheads start pointing fingers again.

Disclaimer: At no point did AI contribute anything to this research project. Copilot can't even calculate basic congruences without making mistakes. People who think AI can do novel math research are delusional.

NOTE: Starting 2026, none of my research will be published. Only people who treat me with respect will be allowed access to my work. And NATO countries/big tech are very low on that list after harassing me for years and treating me like shit. And I garantuee you, I will succeed at finding a polynomial time algorithm. There is no one else alive in this fucking world more determined then me to succeed at this. Fucking losers.

#### To run debug.py" (Prints the linear and quadratic coefficients to solve for 0 in the integers, for use with my paper):</br></br>

To run: python3 debug.py -keysize 12

This basically creates a system of quadratics. Solving them mod p is easy. But there is only one root solution (the factor of N) which solves the system for 0 for any mod p (aka solves it in the integers). Figuring out how to exactly do this quickly is still an ongoing area of research for me. And if a polynomial time algorithm for factorization exists, it is likely done by solving this system of quadratics. Finding a polyomial time algorithm is my ultimate goal, as this would make progress toward solving p = np as well. 

#### To run from folder "CUDA_QS_variant" (WIP):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run:  python3 run_qs.py -keysize 180 -base 4000 -debug 1 -lin_size 10_000_000 -quad_size 50</br></br>

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

Update: Bah wasted day. It doesn't make sense to sieve both sides of the congruences with zx^2+N and zx^2-N. There may be some weird number theoretical trick there.. but thats going to require further research and I'll investigate that some more in the future (wont publish, pay me) .. let me finish what I was doing yesterday.

Update: Fixed a small bug in the PoC. It wasn't using the largest factors to construct a modulus... but i fixed that now so it does use the largest smooth factors to construct a modulus...
I'm going to try again to add p-adic lifting for when we try to find smooths with large factors of an earlier smooth.. because for those smooth candidates.. we dont want to introduce any new large factors.. we only want the large factors from the original smooth there.. and the only way I know how to solve that is to use p-adic lifting...we absolutely cant introduce any new large factors there... because that defeats the entire purpose of what I'm doing... and its the large factors which are the most problematic since the bigger they get they rarer their occurance. I'll do it tomorrow... feeling optimistic... maybe tonight I can sleep.

Like.. it doesnt matter if sieving for smooths with large factors from an earlier smooth is slow due to sieving with p-adic lifting and skipping very large primes that didn't occure in the original smooth.. because anytime we find a smooth with that.. it will dramatically increase the odds of eventually succeeding early at the linear algebra step. And thats the only way to punch past the 110 digit ceiling for quadratic sieve based algorithms.. hmm. I know I tried this before.. but I don't think I did it the proper way last time.. let me try again.
