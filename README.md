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

Update: I made a bit of a mistake... it's not the legendre symbol using the divisor that determines chance of success. Let me actually add some better metrics and see how often the sieve interval actually gets marked per prime and compare. Thinking in the right direction though.. need to find cases when the interval is least saturated by invalid solutions from small primes (i.e flipped to 0).

Update: Doing some more analysis. So we should actually combine both methods... and we should also be able to roughly calculate where these solutions in the interval ends up being.. because if div is the odd exponent factors of an original b-smooth.. then this actually does create some algebraic structure that is going to garantuee solutions... it's just at higher keysizes predicting where these will be becomes a little harder. So I need to figure out some kind formula for this.

Update: Yeap, I'm certain now. Don't touch div.. leave it as it was from the original b-smooth.. don't multiply or divide it... it's the sieve region that we must adjust.. but where the next solutions will be should be calculate-able... what I'm seeing isn't random at all.. it's a clear pattern.. so close to victory now... tomorrow could be the day.

My father is going to the hospital tomorrow... I must succeed. I have to succeed. My father has to see me succeed at this. I won't accept any other reality.

Update: Absolute dreadful days with my father about to have surgery. I did some more brainstorming.. so the best way to approach this I believe is to actually multiply div with squares until there is a small solution that shows up for every prime (well, every prime where div is a quadratic residue). That's how you would do that... so let me get rid of the interval approach of favor of a residue based approach.. then optimize the hell out of it.. and then I'll be done..

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


