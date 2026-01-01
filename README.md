One final time before I do what I must, looking for work, big_polar_bear1@proton.me. Cant work for big tech bc of what happened in 2023. Cant travel to the US anymore bc I would get arrested (hence cant work for US companies either). Cant say I didnt try when all the infosec shitheads start pointing fingers again.

Disclaimer: At no point did AI contribute anything to this research project. Copilot can't even calculate basic congruences without making mistakes. People who think AI can do novel math research are delusional.

NOTE: Starting 2026, none of my research will be published. Only people who treat me with respect will be allowed access to my work. And NATO countries/big tech are very low on that list after harassing me for years and treating me like shit. And I garantuee you, I will succeed at finding a polynomial time algorithm. There is no one else alive in this fucking world more determined then me to succeed at this. Fucking losers.


#### To run debug.py" (Prints the linear and quadratic coefficients to solve for 0 in the integers, for use with my paper):</br></br>

To run: python3 debug.py -keysize 12

This basically creates a system of quadratics. Solving them mod p is easy. But there is only one root solution (the factor of N) which solves the system for 0 for any mod p (aka solves it in the integers). Figuring out how to exactly do this quickly is still an ongoing area of research for me. And if a polynomial time algorithm for factorization exists, it is likely done by solving this system of quadratics. Finding a polyomial time algorithm is my ultimate goal, as this would make progress toward solving p = np as well. 

#### To run from folder "CUDA_QS_variant" (WIP):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run:  python3 run_qs.py -keysize 160 -base 2000 -debug 1 -lin_size 10_000_000 -quad_size 1</br></br>

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

Update: Alright... the code is still very rough now... but it works like this:

We find smooths using x^2+N and square moduli (to reduce to amount of odd exponent factors)
Then we use a factor with odd exponent from the smooth we found as quadratic coefficient and construct a modulus with the rest of the factors and find more smooths with zx^2+N .. this means that we end up with a lot of smooths with a similar factorization. I still need to work a bit on more on which factors to use as quadratic coefficient and which to use as modulus.. so we get a good parabola when generating smooth candidates. But the base idea is implemented.

Anyway... my hope is that that cuts down the required amount of smooths for really really large numbers when working with really really large factor bases. As a proof of concept this works.. but it still requires a bit more work to get it all setup correctly. Fuck it.. all I wanted was some dignity, respect.. an income.. justice for my former manager... instead the west treats me like complete trash and then gaslights me into believing im insane or some shit. These last two years, I have been stuck in a small room, without income, unable to see my friends anymore, socially isolated. And everyone stays silent. Even though I know people know I'm right about my math. Yet I am treated like this. If a visa can be arranged for any country, be it russia, china.. fuck, iran.. I dont care anymore... i'll go. Because every day, I struggle with suicidal thoughts due to being stuck in life here in the west, and if I don't do something, find a place where I'm respected for my work, I know how this will end.
