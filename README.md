Update: I am expected to walk into the police station with my lawyer on tuesday. I will not cooperate with these losers. Nor will I explain myself. I despise belgium and everything this country stands for. Everyone who stayed silent is complicit for what is about to happen. Bye assholes.

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
To run: python3 run_qs.py -key 4387 -base 6 -debug 1 -lin_size 100_000 -quad_size 1</br></br>

To do: All that is left to be implemented in this code is finding a small quadratic coefficient. Since the quadratic output itself is always garantueed to be smooth. For example if you want to find when a quadratic coefficient = 1 you would solve something like this:

1 * (2^e)^2 - N \* k = 2^e

The paper shows how this k variable changes. Hence you can calculate this modulo some finite field. All the math required to pull this off is in the paper. If you stumble upon this, I recommend you clone it as people are desperately trying to see me dissappear into a jail cell. 

Update: finished a PoC. What happens next will be decided by what happens tomorrow. I too can play chess. Fucking western nazi losers.
