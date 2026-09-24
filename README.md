The purpose of this research is to destroy spewers of anti-lgbtq hate (Russia, MAGA) and bring forth the gay future. Also, pete hegseth is a little man and a coward. And so is that pig erdogan and all those other "traditional family values right wing nut jobs" rounding up LGBTQ people simply for existing. Don't travel to Europe. In Europe, nobody gets persecuted for who they are, and people who do shouldn't come here. And maybe people in European politics have forgotten what Europe stands for, but fck them, grow a spine already, because Europe is all that is left of the free world now. People used to take pride in defending freedom and protecting the innocent. What has this world become? Where did courage go? Sick and tired of spineless cowards.

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
To run:  python3 run_qs.py -keysize 50 -base 1000 -debug 0 -lin_size 100 -quad_size 1</br></br>

Minimized a lot of the code now. Gutting all the SIQS style code. 
I'm trying to figure out how to get some type of linear algebra implemented. 
I don't think finding an analogue to "b-smooths" is the correct approach... but rather finding a correct a and k such that b shows up modulo enough primes to ensure squaredness in the integers, for discriminant ab^2+4Nk. This has to be something that can be solved for with linear algebra... a and k basically just multiply/divide b mod p.

UPDATE: DAMNIT!!!!! Let me revert and add my SIQS variant again. I just noticed something. If we find one solution with a large square (via SIQS style sieving for example)... we can use that to find other solutions.... its not just some gimmick or imagined thing. There's an actual pattern here.

UPDATE: Re-uploaded the SIQS style code.... I think I know how to do it now.... let me investigate. 

Update: Just run with the above command. There seems to be this interesting pattern where we are garantueed to find other solutions at the same "k" value when the leading coefficient "a" is a single prime. I'm not sure if the size of the prime matters.. plus I need to do some more testing vs composites. But if this assumption does end up holding true... then its simply a matter of using linear algebra to construct a b-smooth with a single prime... which is massively better then needing a full square. That would be quite interesting... let me do more testing.. might be on the precipice of a monumental breakthrough... at last..

Update: Probably just yields more solutions when "a" is small.. since that generates smaller discriminants. But still seems to be a pattern here that we are garantueed to find more solutions if we find atleast one solution already (which we can find with SIQS style sieving). That I need to start zero-ing in hard on now... because I just know that's something I can leverage somehow.

I'll go for a run.. let me do some analysis of the residues of the coefficients on both sides (so for +Nk and -Nk) when it finds these multiple (but non-trivial) solutions with the psieve() logic. There might be something there that I can leverage. Its definitely the case that when "a" is small.. these multiple solutions are easier to find... but that is also as expected. 

I should also add to my siqs variant support for different k values (which is the multiplier to N in the discriminant)... then when we jump into the psieve logic... there may be a trick to find these other solutions... since k is shared on both sides (for the side with -Nk and with +Nk) .. anyway... I got this... so close now. I feel it. I'm a polar bear.

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

(PERFORMING HEATHEN CYBER RITUALS IN THE BELGIAN MOORLANDS TO DESTROY MY ENEMIES)

