Disclaimer: At no point did AI contribute anything to this research project. Copilot can't even calculate basic congruences without making mistakes. People who think AI can do novel math research are delusional.

NOTE: Starting 2026, none of my research will be published. Only people who treat me with respect will be allowed access to my work. And NATO countries/big tech are very low on that list after harassing me for years and treating me like shit. And I garantuee you, I will succeed at finding a polynomial time algorithm. There is no one else alive in this fucking world more determined then me to succeed at this. Fucking losers.


#### To run debug.py" (Prints the linear and quadratic coefficients to solve for 0 in the integers, for use with my paper):</br></br>

To run: python3 debug.py -keysize 12

This basically creates a system of quadratics. Solving them mod p is easy. But there is only one root solution (the factor of N) which solves the system for 0 for any mod p (aka solves it in the integers). Figuring out how to exactly do this quickly is still an ongoing area of research for me. And if a polynomial time algorithm for factorization exists, it is likely done by solving this system of quadratics. Finding a polyomial time algorithm is my ultimate goal, as this would make progress toward solving p = np as well. 

#### To run from folder "QS_variant" (Standard SIQS with our number theory as backend):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 220 -base 20_000 -debug 1 -lin_size 10_000_000 -quad_size 1</br></br>

See below for an improved way of performing what this PoC does.. I'll delete this Proof of Concept once the PoC for Improved_QS_Variant matures a little more<br><br>

Note: With a large enough -base and lin_size this PoC will find smooths for 110 digits. Albeit very slowly, but this is a highly unoptimized cython PoC. However, to push beyond that into novel terroritory for Quadratic Sieve-based algorithms we need to use quadratic coefficients and p-adic lifting, and that is what the PoC below (Improved_QS_Variant) will be for. 

#### To run from folder "Improved_QS_Variant" (Building smooths from the ground up... see paper on clues to finish the code):</br></br>

To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 14 -base 20 -debug 1 -lin_size 100_000 -quad_size 1</br></br>

Note: Broken for now... this is some work in progress research

#### To run from folder "CUDA_QS_Variant"  (QS variant using mainly vector operations and CUDA support):</br></br>

To run: python3 QScu.py -keysize 50 -base 100 -lin_size 1_000_000 -quad_size 10_000</br></br>

Note: You need to run from the host. I'm using wsl2 and ubuntu with CUDA support. I don't believe this will work from hyper-v.

Update: Bah. I solved some big bottlenecking... I should only really run sum() inside of the GPU lol. Let me add some moduli for the step size of the interval. We also shouldnt run cp.copy() in the gpu lol.. because then we are restricted by GPU memory. I also wonder if we cant calculate multiple quadratic coefficients in one go hmm... because we could just like use the first 16bit for quadratic coefficient = 1, the next 16bit  for quadratic coefficient=2 ... and so on... this actually may work well with a GPU since I assume it doesn't struggle with register size limitations like the CPU. And since we use log base 2 values.. it wouldn't overflow 16bit.










