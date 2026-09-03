The purpose of this research is to destroy spewers of anti-lgbtq hate (Russia, MAGA) and bring forth the gay future. Also, pete hegseth is a little man and a coward.

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
To run: python3 run_qs.py -keysize 40 -base 10_000 -debug 0 -lin_size 1_00 -quad_size 1</br></br>

Update: Bit slow on progress with my father being in the hospital. Changed a few things though. The modulus/interval stepsize for psieve (refered to as div) ...is now square.. and the factors we want to find are now the quadratic coefficient.

Now I can actually change things further and use the quadratic representations of the discriminant and replace the whole interval setup with that... since we just need to find quadratics that evaluate to 0... as this maps to a discriminant that is a square residue modulo any prime.

Update: fixed the code that shows the quadratics and their roots at line 2184. Next I'll get rid of the legendre based interval and we're going to build an algorithm using these quadratics.. 

UPDATE: OH YEA. WAIT. I should really just be able to reveal the factor of N once I have a b-smooth with a large enough square in it using that logic at line 2184 (2214 in new version). That would reduce factorization to finding a single B-smooth and put RSA-1024 well within scope. Ok.. this won't be very difficult, give me a day. I went to visit my father at the hospital and when I came back the battery of my research laptop stopped working. Most likely some morrons tried to insert an implant or some shit.

Update: Works with squares.. now to get it to work with non-squares.. this can't be so hard... I'll have it before the end of the week. I just know it can be done. 

Update: I might have it. Instead of using some random large prime.. we should take whats non-square in the b-smooth, ie if we're looking for a discriminant where a\*x^2+4\*N\*k is square and a is not square, then we lift using a (and use CRT if a has multiple primes ofcourse). That could work hmm.. because then this a becomes a 0 solution..basically canceling it out and leaving x^2+4\*N\*k... I like this idea. I'll try it tomorrow. If that works I probably need to rework my PoC and just use the psieve approach.. because all we will need to do is find a single b-smooth with it with a large enough square in it... or could probably also use the SIQS method to sieve for these types of b-smooth... which might be faster.. and then just finish with the approach I just mentioned.

ps: Hello Russia, Hello loser pete hegseth. How does it feel to get destroyed by a single polar bear? Retrocausality b*tches. You lost before the game even began.

Update: One step closer to completion yet again. Replaced the square root over a large prime with p-adic lifting now. However, rather then lifting some small prime where "a" is a quadratic residue (a as in a\*x^2+4\*N\*k .. or see how it is used in psieve_process_interval()), we need to solve the case where the prime divides a... I'm guessing using the theory about singular and non-singular roots that I described in the paper.. once that works.. hopefully I can solve cases where "a" is non-square... and IF that works.. I have an immediate breakthrough of historical proportions and RSA-1024 would fall for sure.. and who knows what else.

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


