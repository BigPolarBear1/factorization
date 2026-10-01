The purpose of this research is to destroy spewers of anti-lgbtq hate (Russia, MAGA) and bring forth the gay future. Also, pete hegseth is a little man and a coward. And so is that pig erdogan and all those other "traditional family values right wing nut jobs" rounding up LGBTQ people simply for existing. Don't travel to Europe. In Europe, nobody gets persecuted for who they are, and people who do shouldn't come here. And maybe people in European politics have forgotten what Europe stands for, but fck them, grow a spine already, because Europe is all that is left of the free world now. People used to take pride in defending freedom and protecting the innocent. What has this world become? Where did courage go? Sick and tired of spineless cowards.

Disclaimer: No AI was used for any of this. I tried using it at times, but the only real application I got out of it was reviewing my paper for minor mistakes (it has however not written a single word in my paper). A context aware search engine.. good for learning things or producing simple code for well known mathematical functions, but it quickly falls apart when it gets into research territory.

#### References (re-used many of the core number theoretical functions from these PoCs to fit my own algorithm): 
https://stackoverflow.com/questions/79330304/optimizing-sieving-code-in-the-self-initializing-quadratic-sieve-for-pypy
https://github.com/basilegithub/General-number-field-sieve-Python 

#### About the paper
Math paper is a work in progress. Ignore the final chapter for now.. that one I'll rewrite if and when I can get "psieve" below working correctly.

#### To run from folder "psieve" WORK IN PROGRES...extremely early version:</br>
To build: python3 setup.py build_ext --inplace</br>
To run:  python3 run_qs.py -keysize 70 -base 10_000 -debug 0 -lin_size 100 -quad_size 1</br></br>

This starts with an SIQS variant where the stripped away factors are square to minimize the odd exponent factors in the bsmooth and after that we use the research from the paper in psieve(). This then utilizes a two sided approach and looks for similar b-smooths by calculating quadratic residues and applying hensel's lifting (I need to fix some of that hensel code still, its just bruteforcing roots mod p^e for now). 

To do: Performance now is decent enough. There are two options now to finish this project:

Option 1: Find a way to select a better "a" and "k" parameter. We can just add squares to "a", this won't impact the algorithm but will change how the solutions calculated in "resmaps" are shifted in psieve_build_interval(). Optimizing these parameters should be possible.. but I need to dig a little deeper into how to best approach this. And doing that "shifting" of solutions from "resmaps" until a shared solution is found can probably be done with some type of linear algebra.. but aside from guassian elimination over GF(2) I havnt had too much exposure to linear algebra yet.

Option 2: Find a way, to for example use hensel's lifting, to just straight up calculate solutions. 

There is some small things in psieve() that can increase performance a little more.. like the roots for the primes that divide the discriminant are trivially precalculated and we can also use graycodes there. Thats minor speed increases though.

Update: Quickly added an "a_mul" loop in psieve() which adds squares to "a". This doesn't change the primes where a is a quadratic residue, but does generate a unique interval. Hence this is a perfect variable to build an optimizer function for... let me try to brainstorm how I can do this.

Update: Doing some more thinking.. if we're not at the correct "k" then it doesn't matter, we will never find an "a_mul" that will generate solutions in our interval. This line:  if kronecker_symbol(a,n*k) == -1: used to skip bad "k" values does seem to filter out some of the bad "k" values.. if I can narrow it down further.. thats going to be a massive speedboost and only then I should start thinking about optimizing a_mul. So let me explore that first...

#### To run from folder "Coefficient_Sieve" (For use with the paper):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run:  python3 run_qs.py -keysize 40 -base 50 -debug 1 -lin_size 10_000 -quad_size 100</br>

Just demonstrates the math from the paper using quadratics. For educational purposes. And rather then taking a square root over a large prime we can also just calculate the discriminant. But this demonstrates the interesting relation between these quadratics and the factors of N.

#### To run debug.py" (Prints the linear and quadratic coefficients to solve for 0 in the integers, for use with my paper):</br></br>

To run: python3 debug.py -keysize 12

This basically creates a system of quadratics. Solving them mod p is easy. But there is only one root solution (the factor of N) which solves the system for 0 for any mod p (aka solves it in the integers). Figuring out how to exactly do this quickly is still an ongoing area of research for me. And if a polynomial time algorithm for factorization exists, it is likely done by solving this system of quadratics. 

(PERFORMING HEATHEN CYBER RITUALS IN THE BELGIAN MOORLANDS TO DESTROY MY ENEMIES)

