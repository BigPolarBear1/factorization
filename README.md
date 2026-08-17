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

Update: Pushed some updates. It's a bit shitty. At line 1168 ( if len(sqr)==0 or math.gcd(div,prime)!=1:) .. these cases are the "lazy" calculation. We could just calculate legendre symbols for all the primes. But what I really want to do is look at the cases where we can actually construct a quadratic.. where if it has a root solution then the discriminant will also be a quadratic residue mod p or divide by p (depending if its singular or non-singular).. so I want to lift these quadratics to p^e next.. create a residue map of these... then use some fast math to figure out the existence of B-smooths with large squares in it.. I think I know how.. it should be possible.. but I should really implement p-adic lifting...

But anyway.. the code shows how to get quadratics now.. even when we divide the discriminant in advance by "div". Because only looking at full squares doesn't give us a large enough set.. but this way we can also consider small-ish multiples of a square. Which still gives us a rank reduction during the linear algebra step (I believe rank reduction is the correct term.. with that I just mean, require less B-smooths to succeed).

The main thing that needs to be figured out is if there exists a way to find these "almost" large squares without using a sieve interval and using residues instead. I'll delegate this research to my weekends for a little while.. I've also been working on some other shit... not going to waste my life slowly perishing in a tiny attic room.

Update: Alright, after spending a few days on another project (preparing the destruction of Russia, haha, suka blyat), its weekend.. time to do some math.

Update: Quickly changed a few things so it just uses a residue map. Going to implement p-adic lifting. Will probably just implement some placeholder bruteforce logic for a while.. doesn't really matter for testing purposes plus I can use it later when implementing proper p-adic lifting to make sure it's working and not skipping solutions. I hope with p-adic lifting I can skip having to calculate legendre symbols for the primes that I cant convert into a quadratic because div isn't a quadratic residue... and then all that's left really is finding a "good div" that maps to a valid solution in my residue map for each prime. Easy enough. I can do this actually.. lol. Fuck them all. Fuck Russia, Fuck MAGA. You went after the gays and queers. Get destroyed by one now. Fucking losers. I piss on all of you. I fucking despise you. Fucking traitors. Fucking russian lapdogs. R*tards.

Update: Oh shit.. I was doing analysis.. there are "divisors" that will yield an enormous amount of b-smooths... oh shit... and I think I can calculate when this happens. oh shit.... I swear if this works like I think it will.. go to hell for forcing me to burn my work when it could have been used to fucking pown russia instead. 

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


