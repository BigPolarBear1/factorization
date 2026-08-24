Disclaimer: No AI was used for any of this, except for reviewing my paper these last 2 weeks, but AI has not written a single sentence.
None of the code is written by AI either, except for one or two functions like lift_root2(),
which is just hensel, something I had already implemented before but with a coefficient list as input.
I don't believe AI is quite there yet to do math research. It's very rigid and can't think outside the bounds of existing literature
and often just makes very dumb conceptual mistakes and it has a total lack of creating abstractions.
It is a tool, well suited for basic tasks, nothing more.
This research is also still ongoing, and especially some of the things stated in the last chapter might be missing the mark. 
I'll also be properly learning about number fields now, and see how all of that can be fit into my work.

#### References (re-used many of the core number theoretical functions from these PoCs to fit my own algorithm): 
https://stackoverflow.com/questions/79330304/optimizing-sieving-code-in-the-self-initializing-quadratic-sieve-for-pypy
https://github.com/basilegithub/General-number-field-sieve-Python 

#### About the paper
Math paper is a work in progress. Ignore the final chapter for now.. that one I'll rewrite if and when I can get "psieve" below working correctly.

#### To run from folder "psieve" WORK IN PROGRES...extremely early version:</br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 45 -base 10_000 -debug 0 -lin_size 10_000 -quad_size 1</br></br>

Uploaded a very first commit. This basically merges cuda_qs_variant and coefficient_sieve. Calling into coefficient sieve (see psieve()) whenever a b-smooth is found. 
Whats left to do is setting the proper sieve region to check.. because now its just doing 0 to 5000 in a 2d plane.. which is not going to work if we increase key size. Let me do some analysis how to find the best region. In addition we could also tune "div" that's being passed into psieve as another parameter. So this seems like a very promising research direction so far. Lets see.

Update: Let me focus next on either optimizing the sieve region or divisor value in psieve(). If I can find a divisor value that's going to have more "valid" solutions on average per prime.. then thats going to create a better sieve interval. So it can definitely be optimized.. question is if it can be optimized in a way that actually matters I guess. May also need to include p-adic lifting for this.

Update: So the whole thing about lifting solutions to quadratics.. and how these quadratics link to the squaredness of the discriminant.. we can lift solutions for small primes.. and the density of those solutions sets will actually let us optimize this divisor I believe. Its really those very small primes that matter most.. should be able to get it done this weekend hopefully. Just add a function to generate some type of scoring for each divisor.. then either run the psieve logic when we find a divisor with a good score or just multiply the divisor with squares until it sits in a good score.. 

Update: Did some late night testing.. seems the number of solutions, it depends whether or not the divisor is a square residue mod p... and I also know that if I multiply primes with eachother to construct a divisor.. their legendre symbols carry over into the product.. hmm.. this hence shouldn't be very complicated.
What this means, is that in psieve, where it either checks the hashmap or calculates a legendre symbol based on whether or not the divisor is a square residue.. those two code paths actually yield very different amounts of solutions.. with the non-QR case (Where we just calculate legendre on the discriminant) being more favorable.
Hence if we simply maximize that.. find a divisor that is a non-QR for as many small primes as possible (and we can even integrate lifting) .. then that will yield a much better interval and higher probability of hitting a b-smooth with a large square in it. 
OH SHIT. I'm starting to see it now. 
Tomorrow could be the day. Tomorrow could be the day THIS ENDS. This 3+ years nightmare. Lets go. 

#### To run from folder "Coefficient_Sieve" (For use with the paper):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run:  python3 run_qs.py -keysize 40 -base 50 -debug 1 -lin_size 10_000 -quad_size 100</br>

Just demonstrates the math from the paper using quadratics. For educational purposes. And rather then taking a square root over a large prime we can also just calculate the discriminant. But this demonstrates the interesting relation between these quadratics and the factors of N.

edit: You can modify this to use the leading coefficient to divide the discriminant and find B-smooths with a large square in it... I will experiment with it soon.. but first I want to spent a few more days just trying direct computation of a solution without sieving, because that could potentially straight up break factorization and I must know for sure if it is possible or not first.

#### To run from folder "CUDA_QS_variant" (Failed Experiment):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run: To run:  python3 run_qs.py -keysize 240 -base 100_000 -debug 1 -lin_size 100_000_000 -quad_size 100</br></br>
 
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

This was an attempt at finding smooths with similar factorization using an SIQS variant. By using quadratic coefficients. But it didnt end up working as I had hoped so I abondoned this approach, but perhaps someone will get some use out of it.

#### To run debug.py" (Prints the linear and quadratic coefficients to solve for 0 in the integers, for use with my paper):</br></br>

To run: python3 debug.py -keysize 12

This basically creates a system of quadratics. Solving them mod p is easy. But there is only one root solution (the factor of N) which solves the system for 0 for any mod p (aka solves it in the integers). Figuring out how to exactly do this quickly is still an ongoing area of research for me. And if a polynomial time algorithm for factorization exists, it is likely done by solving this system of quadratics. 


