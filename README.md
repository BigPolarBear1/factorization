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
To run: python3 run_qs.py -keysize 14 -base 20 -debug 1 -lin_size 100_000 -quad_size 1</br></br>

UPDATE Saturday 20 DEC: I'll upload a better PoC shortly as I work out the details.

After spending the afternoon thinking about this, I think the following must be done now:

We must combine quadratic polynomials, such that z\*x\*2 = y.
I need to work out the exact math for this, so that when we combine them, no extra factors that don't already occur in the quadratic outputs that we are combining get added. 
Ideally rather than using quadratic polynomials to perform sieving such that the discriminant is smooth, we should combine them such that we generate a square and can take the GCD.. so somehow that combining of quadratic polynomials must be fed into guassian elimination over GF(2). So that's what I'll begin exploring now. I'm fairly sure it's do-able though. This does feel like the right approach. The main issue will be approaching it in a way that makes sense.. rather then desperately trying to mimic steps in the number field sieve algorithm.. because that's just going to confuse me.

