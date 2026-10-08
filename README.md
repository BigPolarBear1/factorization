The purpose of this research is to destroy spewers of anti-lgbtq hate (Russia, MAGA) and bring forth the gay future. Also, pete hegseth is a little man and a coward. And so is that pig erdogan and all those other "traditional family values right wing nut jobs" rounding up LGBTQ people simply for existing. Don't travel to Europe. In Europe, nobody gets persecuted for who they are, and people who do shouldn't come here. And maybe people in European politics have forgotten what Europe stands for, but fck them, grow a spine already, because Europe is all that is left of the free world now. People used to take pride in defending freedom and protecting the innocent. What has this world become? Where did courage go? Sick and tired of spineless cowards.

Disclaimer: No AI was used for any of this except for generating some of the code for well-known mathematical functions (it is indicated in my code when it's AI generated), and others are copy pasted from the references listed below, but all the novel logic is written by me by hand. I've also not used it for any of the research itself. Most of the foundation I already came up with 3 years ago (if people still remember that bad LLL implementation I wrote way back.. that one abused exactly the same principles I'm still exploiting here.. it just took time for me to mature mathematically and make the correct connections with the framework I had come up with). 

#### References (re-used many of the core number theoretical functions from these PoCs to fit my own algorithm): 
https://stackoverflow.com/questions/79330304/optimizing-sieving-code-in-the-self-initializing-quadratic-sieve-for-pypy
https://github.com/basilegithub/General-number-field-sieve-Python 

#### About the paper
Math paper is a work in progress. Ignore the final chapter for now.. that one I'll rewrite if and when I can get "psieve" below working correctly.

#### To run from folder "psieve" WORK IN PROGRES...extremely early version:</br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 100 -base 1000 -lin_size 10_000 -quad_size 1 -mode nfs</br></br>

Note: For a real world implementation you should have the SIQS variant and Psieve() running in two different threads. SIQS collects b-smooths, Psieve then uses a quadratic residue based approach to try and complete the linear algebra step much sooner. Still a work in progress and a lot still need to be done. Performance of Psieve is still subpar with what I suspect it should be able to achieve. An alternative approach would be to just have Psieve grind B-smooths for arbitrary small 'a' values to reduce the rank of the matrix use for the linear algera step.. since this lets best_a_mul() optimize the interval more.. but there is some quadratic reciprocity related stuff happening that dictates when this will work or not, which I'll need to study a bit deeper.

This starts with an SIQS variant where the stripped away factors are square to minimize the odd exponent factors in the bsmooth and after that we use the research from the paper in psieve(). This then utilizes a two sided approach and looks for similar b-smooths by calculating quadratic residues.

Update: It's now splitting the odd exponent factor part of a b-smooth between a and u. Where ay^2+4Nk = ub^2 (and visa versa). Because we can do this: (ay)^2+4Nka = aub^2... hence if the existing of au is proven via the SIQS main logic, then it can be split into a and u in the psieve logic. Now how exactly to split the factors in au and distribute them over a and u is still something I'm looking at. But since (ub)^2-4Nku = uay^2 ..the this means that ku must be a valid multiplier to N. Things whose residues we can compute... etc. I'm still trying to think how to approach all of this in an algorithmic way that actually works rather then just sieving.

#### To run from folder "Coefficient_Sieve" (For use with the paper):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run:  python3 run_qs.py -keysize 40 -base 50 -debug 1 -lin_size 10_000 -quad_size 100</br>

Just demonstrates the math from the paper using quadratics. For educational purposes. And rather then taking a square root over a large prime we can also just calculate the discriminant. But this demonstrates the interesting relation between these quadratics and the factors of N.

#### To run debug.py" (Prints the linear and quadratic coefficients to solve for 0 in the integers, for use with my paper):</br></br>

To run: python3 debug.py -keysize 12

This basically creates a system of quadratics. Solving them mod p is easy. But there is only one root solution (the factor of N) which solves the system for 0 for any mod p (aka solves it in the integers). Figuring out how to exactly do this quickly is still an ongoing area of research for me. And if a polynomial time algorithm for factorization exists, it is likely done by solving this system of quadratics. 

(PERFORMING HEATHEN CYBER RITUALS IN THE BELGIAN MOORLANDS TO DESTROY MY ENEMIES)

