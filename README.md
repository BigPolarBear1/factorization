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

Only one thing left to do. In psieve_process_interval() "a" is square. And then we can simply calculate a root in a finite field for some prime such that "a" is a quadratic residue for that prime. However, if we can make these calculations for some prime that divides "a" .. then the discriminant a\*x^2+4\*N\*k then "a" becomes a 0 solution and we might be able to reveal the factor of N even when "a" is not square. That would reduce factorization to simply finding some discriminant a\*x^2+4\*N\*k that is square and where the factorization of a is known, but not necessarily square. If I can achieve this one thing, it would be a complete break in the factorization problem. Might be able to do it using the singular vs non-singular theory in my paper. 

Update: Did some initial analysis.. there definitely is a good chance that it can be done. I'll continue shortly although the next two days are a bit busy.

Update: EUREKA!!!!!!!!!!!! I figured it out!!!! If for the discriminant a\*x^2+4\*N\*k, "a" is non-square.. we need to calculate everything in the finite field of some divisor of k :). I got it! ITS DEFINITELY POSSIBLE!

Update: Yeap. I got it. This was the correct way to approach it. FINALLY. I'm about to break reality. This is going to be fcking surreal. And either everyone underestimates me and this is going to come out of nowhere.. or people know.. and I guess.. then they got what they deserve. Have to do some shit tomorrow.. but will start uploading a more finished version soon. 
 
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


