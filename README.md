Disclaimer: At no point did AI contribute anything to this research project. Copilot can't even calculate basic congruences without making mistakes. People who think AI can do novel math research are delusional.

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

The paper shows how this k variable changes. Hence you can calculate this modulo some finite field. All the math required to pull this off is in the paper.

Update: Still have my freedom. Which means that the road to finish my work is now open. One final mad dash across the finish line now..  and if I succeed, I don't care anymore about the consequences. These have been the hardest years of my life.  Lost it all in the pursuit of this. I have nothing left to lose now. All that is left now is factorization.

Update: I started taking this apart today. I actually noticed something. You can actually derive p+q from a single smooth number. I know I explored this earlier in september. But I actually see how it works now. Just one smooth number, and boom, you have factorization due to the fact we can use that to immediatly find p+q. Atleast I have gotten these calculations to work already with something of the form: 1 * (prime^e)^2 - N \* k = prime^e   ... let me begin generalizing these findings and I'll upload them soon....
