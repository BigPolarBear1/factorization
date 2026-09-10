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
To run:  python3 run_qs.py -keysize 40 -base 1000 -debug 0 -lin_size 1_000 -quad_size 1</br></br>

Update: Alright, rather then taking a square root of a non-square leading coefficient I'm just going to try and find some "k" such that we have a full square on both sides. The premise being that as long as we know the factorization of both linear coefficients, and both linear coefficients are large enough... then we can calculate the correct "k" value (k as in x^2-4Nk).

Right now its only calculating possible k values for one linear coefficient. BUT... if we know the factorization of the other linear coefficient, then we know that whatever k we use, must also be a valid k for every prime that divides that other linear coefficient. This feels right now. Should be it. Get fcked losers.

AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAH. Started refactoring. Fuck I'm an idiot. If you know the factorization of both linear coefficients, and both are large enough.. then ofcourse you can find the solution. It literally becomes finding a correct k plus you can add the modulus to the linear coefficient of the near square b-smooth. Its not difficult. Expect a solution to come online soon.

Update: Going to run 20k. Then few more hours on my math. As long as I know the factorization of the root that produced a B-smooth with a large enough square in it.. then I can work in the reverse direction and compute residues to make it completely square. I know it can be done. Watch me. #nevergiveup

Update: Alright. Added the "main idea" very roughly to psieve_factor() ... lots more work needs to be done. I'll do that tomorrow. It's really about trying to find the root that generated the original b-smooth from the other side around. And.. we know how that works... how quadratics with a solution are a QR when taking the discriminant. Code is very bad right now.. but going to expand on this idea quickly now and streamline it. THIS IS IT! I KNOW IT! I GOT IT! LETS GO! TIME TO FINISH THIS! MY FATHER WILL SEE ME SUCCEED AT THIS.

#### To run from folder "Coefficient_Sieve" (For use with the paper):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run:  python3 run_qs.py -keysize 40 -base 50 -debug 1 -lin_size 10_000 -quad_size 100</br>

Just demonstrates the math from the paper using quadratics. For educational purposes. And rather then taking a square root over a large prime we can also just calculate the discriminant. But this demonstrates the interesting relation between these quadratics and the factors of N.

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

