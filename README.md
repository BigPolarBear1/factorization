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

Update: Decided to just take the path of least resistance and build on something that I already know will work, and that is the coefficient Sieve PoC. However, this time around, I will use the leading coefficient to divide the discriminant.. this way we can sieve for small multiples of a large square.. so basically a different way to sieve for B-smooths... I'll add support to divide the discriminant with the leading coefficient soon... 

Update: Ok, that wasn't too hard.. I'll upload how to do it after I come back from running. Kinda struggled for a few days figuring out how to divide the discriminant while using a setup like coefficient_sieve does..  but you just multiply your linear coefficient by the square root of the divisor mod p. So I'll upload that soon.. then I need to figure out how I can use everything else to even quicker find B-smooths with a large square in it by using residues and perhaps p-adic lifting. I'm sure a superior QS variant can be constructed this way.

Update: Just really quickly uploaded a version using that math... it's not very clean, but I'll improve it now. The real thing I need to figure out is if there exists an algorithm to more quickly find these cases.. so an algorithm centered around quadratics and the legendre symbol of their discriminant mod p.. rather factorization over a factor base as traditional QS does..

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


