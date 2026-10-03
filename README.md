The purpose of this research is to destroy spewers of anti-lgbtq hate (Russia, MAGA) and bring forth the gay future. Also, pete hegseth is a little man and a coward. And so is that pig erdogan and all those other "traditional family values right wing nut jobs" rounding up LGBTQ people simply for existing. Don't travel to Europe. In Europe, nobody gets persecuted for who they are, and people who do shouldn't come here. And maybe people in European politics have forgotten what Europe stands for, but fck them, grow a spine already, because Europe is all that is left of the free world now. People used to take pride in defending freedom and protecting the innocent. What has this world become? Where did courage go? Sick and tired of spineless cowards.

Disclaimer: No AI was used for any of this. I tried using it at times, but the only real application I got out of it was reviewing my paper for minor mistakes (it has however not written a single word in my paper). A context aware search engine.. good for learning things or producing simple code for well known mathematical functions, but it quickly falls apart when it gets into research territory.

#### References (re-used many of the core number theoretical functions from these PoCs to fit my own algorithm): 
https://stackoverflow.com/questions/79330304/optimizing-sieving-code-in-the-self-initializing-quadratic-sieve-for-pypy
https://github.com/basilegithub/General-number-field-sieve-Python 

#### About the paper
Math paper is a work in progress. Ignore the final chapter for now.. that one I'll rewrite if and when I can get "psieve" below working correctly.

#### To run from folder "psieve" WORK IN PROGRES...extremely early version:</br>
To build: python3 setup.py build_ext --inplace</br>
To run:  python3 run_qs.py -keysize 60 -base 10_000 -debug 0 -lin_size 100 -quad_size 1</br></br>

This starts with an SIQS variant where the stripped away factors are square to minimize the odd exponent factors in the bsmooth and after that we use the research from the paper in psieve(). This then utilizes a two sided approach and looks for similar b-smooths by calculating quadratic residues and applying hensel's lifting (I need to fix some of that hensel code still, its just bruteforcing roots mod p^e for now). 

Update: I added a new datastructure called "resmaps2" this contains all the residues for the linear coefficient of the quadratic on the other side. I've added this on purpose.. tomorrow I'll make sure the indexes of that coefficient list also match the one from resmaps and visa versa. We know that a_mul is going to divide the residues from resmaps2... so we'll need to use those residues to calculate a good "a_mul".

Update: PoC is a bit convoluted now and a little slower then yesterday's version, but this demonstrates how to switch between sides. I do see now how I can optimize that a_mul variable... the primes that are used to mark wether a discriminant is a QR after dividing out a... we can construct multiple small moduli there. For a correct solution there will be one overlapping solution. However with a_mul we can divide out coefficients on the ab^2+4Nk side. I see it now. This is about to end very soon... people saw how determined I was, yet I was left to rot in an attic room for year. I have no regrets. 

#### To run from folder "Coefficient_Sieve" (For use with the paper):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run:  python3 run_qs.py -keysize 40 -base 50 -debug 1 -lin_size 10_000 -quad_size 100</br>

Just demonstrates the math from the paper using quadratics. For educational purposes. And rather then taking a square root over a large prime we can also just calculate the discriminant. But this demonstrates the interesting relation between these quadratics and the factors of N.

#### To run debug.py" (Prints the linear and quadratic coefficients to solve for 0 in the integers, for use with my paper):</br></br>

To run: python3 debug.py -keysize 12

This basically creates a system of quadratics. Solving them mod p is easy. But there is only one root solution (the factor of N) which solves the system for 0 for any mod p (aka solves it in the integers). Figuring out how to exactly do this quickly is still an ongoing area of research for me. And if a polynomial time algorithm for factorization exists, it is likely done by solving this system of quadratics. 

(PERFORMING HEATHEN CYBER RITUALS IN THE BELGIAN MOORLANDS TO DESTROY MY ENEMIES)

