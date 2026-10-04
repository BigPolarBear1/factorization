The purpose of this research is to destroy spewers of anti-lgbtq hate (Russia, MAGA) and bring forth the gay future. Also, pete hegseth is a little man and a coward. And so is that pig erdogan and all those other "traditional family values right wing nut jobs" rounding up LGBTQ people simply for existing. Don't travel to Europe. In Europe, nobody gets persecuted for who they are, and people who do shouldn't come here. And maybe people in European politics have forgotten what Europe stands for, but fck them, grow a spine already, because Europe is all that is left of the free world now. People used to take pride in defending freedom and protecting the innocent. What has this world become? Where did courage go? Sick and tired of spineless cowards.

Disclaimer: No AI was used for any of this. I tried using it at times, but the only real application I got out of it was reviewing my paper for minor mistakes (it has however not written a single word in my paper). A context aware search engine.. good for learning things or producing simple code for well known mathematical functions, but it quickly falls apart when it gets into research territory.

#### References (re-used many of the core number theoretical functions from these PoCs to fit my own algorithm): 
https://stackoverflow.com/questions/79330304/optimizing-sieving-code-in-the-self-initializing-quadratic-sieve-for-pypy
https://github.com/basilegithub/General-number-field-sieve-Python 

#### About the paper
Math paper is a work in progress. Ignore the final chapter for now.. that one I'll rewrite if and when I can get "psieve" below working correctly.

#### To run from folder "psieve" WORK IN PROGRES...extremely early version:</br>
To build: python3 setup.py build_ext --inplace</br>
To run:  python3 run_qs.py -keysize 50 -base 10_000 -debug 0 -lin_size 100 -quad_size 1</br></br>

This starts with an SIQS variant where the stripped away factors are square to minimize the odd exponent factors in the bsmooth and after that we use the research from the paper in psieve(). This then utilizes a two sided approach and looks for similar b-smooths by calculating quadratic residues and applying hensel's lifting (I need to fix some of that hensel code still, its just bruteforcing roots mod p^e for now). 

Update: Re-uploaded yesterday's version after fcking around today. Just realized something.

1. Should use bit-packing for the interval. Big improvement.
2. The smaller "a" is ... the more we can add padding with that "a_mul" variable.. and as long as jacobi(-nk,sqrt(a)) is a quadratic residue, we can add more solutions to the interval trivially like this. Seeing something nice here now... time to finish this. Also we can generate a small "a" using linear algebra on multiple b-smooths rather then directly using sieving results from the SIQS variant.

Anyway, hope the guys at the NSA have a good day, dont blow your brains out. Go to the arctic like a normal person instead.

Update: Added some a_mul optimizer function. Does seem to mostly increase the solutions in the interval. Although it needs a little more research. Good enough for starters. Next I will write an intermediate linear algebra step that will produce nearly square b-smooths. and use that instead of direct sieving results from the SIQS variant. Because we need room to build up that a_mul...

Update: Quickly re-added yesterdays CRT code at the bottom of psieve. This right now is just slowing down the algorithm... instead this should be a single prime that we lift with hensel..... then I should experiment if this can reveal solutions outside of the interval. etc. Just adding it to show how thats done... for now. 

#### To run from folder "Coefficient_Sieve" (For use with the paper):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run:  python3 run_qs.py -keysize 40 -base 50 -debug 1 -lin_size 10_000 -quad_size 100</br>

Just demonstrates the math from the paper using quadratics. For educational purposes. And rather then taking a square root over a large prime we can also just calculate the discriminant. But this demonstrates the interesting relation between these quadratics and the factors of N.

#### To run debug.py" (Prints the linear and quadratic coefficients to solve for 0 in the integers, for use with my paper):</br></br>

To run: python3 debug.py -keysize 12

This basically creates a system of quadratics. Solving them mod p is easy. But there is only one root solution (the factor of N) which solves the system for 0 for any mod p (aka solves it in the integers). Figuring out how to exactly do this quickly is still an ongoing area of research for me. And if a polynomial time algorithm for factorization exists, it is likely done by solving this system of quadratics. 

(PERFORMING HEATHEN CYBER RITUALS IN THE BELGIAN MOORLANDS TO DESTROY MY ENEMIES)

