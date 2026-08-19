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
Math paper is a work in progress. The final chapters are a bit rushed and building an algorithm around p-adic lifting isnt as straight forward as I had assumed. I do think there is an angle there I can exploit, but I'll do some further experimentation first and get a working PoC before I make edits to the paper again.

#### To run from folder "psieve" WORK IN PROGRES...extremely early version:</br>
To run:  python3 psieve.py -base 30 -keysize 40 -lin_size 10_000</br>

Update: Doing one final hail mary to complete this one. So this extends on what I had before in Coefficient Sieve by dividing the discriminant by a divisor d. So we can find b-smooths which are products of d and a large square by sieving with legendre symbols or checking for the existence of a root solution to a quadratic (which is the same thing, since it means the discriminant is a QR).
Now the key thing is, when divisor d is composite, there is a direct relation between the legendre symbols of each factor. So in theory it should be possible to sieve with a divisor set to primes and combine them so that the legendre symbols all indicate square residues.. while keeping coefficients small (which was my main issue in previous attempts).. I checked the math today and I actually see a way it can be done. It's complicated as hell though. But just going all-in this week..
I want my father to see I succeeded incase anything happens with his surgery next week... all that matters in my life right now. Everything else can wait now.

I am fairly sure I got it figured out and can finish it this week.. if someone wants to stop me from uploading more PoCs... because I know people must know I'm right and what I'm about to upload next.. just contact me and let my father know what I succeeded at before he has surgery next week.
If anything happens to my father.. and his last memories of me are this bullshit, me being unemployed for years... unable to even secure job interviews.. that will be something that can absolutely never be forgiven and that will definitely change me deeply as a person.

Update: Very hard to focus today with my father about to have a really invasive major surgery. Plus if it does end up being malign cancer.. he has about a 50% survival rate, even with a succesful surgery (because the cancer cells may have spread already to lungs, liver, etc and he will need chemo if it ends up being malign to reduce that possibility..). I feel numb inside my head. I did notice, if two prime divisors (div in the PoC) yield b-smooths within a coefficient region.. then its product will also yield a b-smooth within the coefficient region.. or atleast this seems to work extremely reliably beyond statistical chance... one thing left to verify is unique of b-smooths found by just taking the product of divisors... and what I can do once I verify this.. is just plug this logic into cuda_qs_variant to generate more smooths from just a couple of them. It's actually quite straight forward in that case... 

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


