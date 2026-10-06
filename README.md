The purpose of this research is to destroy spewers of anti-lgbtq hate (Russia, MAGA) and bring forth the gay future. Also, pete hegseth is a little man and a coward. And so is that pig erdogan and all those other "traditional family values right wing nut jobs" rounding up LGBTQ people simply for existing. Don't travel to Europe. In Europe, nobody gets persecuted for who they are, and people who do shouldn't come here. And maybe people in European politics have forgotten what Europe stands for, but fck them, grow a spine already, because Europe is all that is left of the free world now. People used to take pride in defending freedom and protecting the innocent. What has this world become? Where did courage go? Sick and tired of spineless cowards.

Disclaimer: No AI was used for any of this except for generating some of the code for well-known mathematical functions (it is indicated in my code when it's AI generated), and others are copy pasted from the references listed below, but all the novel logic is written by me by hand. I've also not used it for any of the research itself. Most of the foundation I already came up with 3 years ago (if people still remember that bad LLL implementation I wrote way back.. that one abused exactly the same principles I'm still exploiting here.. it just took time for me to mature mathematically and make the correct connections with the framework I had come up with). 

#### References (re-used many of the core number theoretical functions from these PoCs to fit my own algorithm): 
https://stackoverflow.com/questions/79330304/optimizing-sieving-code-in-the-self-initializing-quadratic-sieve-for-pypy
https://github.com/basilegithub/General-number-field-sieve-Python 

#### About the paper
Math paper is a work in progress. Ignore the final chapter for now.. that one I'll rewrite if and when I can get "psieve" below working correctly.

#### To run from folder "psieve" WORK IN PROGRES...extremely early version:</br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 100 -base 10_000 -debug 0 -lin_size 100 -quad_size 1</br></br>

Note: For a real world implementation you should have the SIQS variant and Psieve() running in two different threads. SIQS collects b-smooths, Psieve then uses a quadratic residue based approach to try and complete the linear algebra step much sooner. Still a work in progress and a lot still need to be done. Performance of Psieve is still subpar with what I suspect it should be able to achieve. An alternative approach would be to just have Psieve grind B-smooths for arbitrary small 'a' values to reduce the rank of the matrix use for the linear algera step.. since this lets best_a_mul() optimize the interval more.. but there is some quadratic reciprocity related stuff happening that dictates when this will work or not, which I'll need to study a bit deeper.

This starts with an SIQS variant where the stripped away factors are square to minimize the odd exponent factors in the bsmooth and after that we use the research from the paper in psieve(). This then utilizes a two sided approach and looks for similar b-smooths by calculating quadratic residues.

Update: Getting very close now. There's some additional research to be done still... its almost finished though. One thing I want to research is the commented out code at the bottom of psieve() where it checks roots by CRT. I need to rewrite that into a smart hensel implementation and see if it can find valid solutions that lie outside the bounds of the interval easily somehow. There is also a very deep pattern based on quadratic reciprocity here that dictates when psieve succeeds that needs further exploration. This is why you cannot set "a" as an aribtrary small prime in psieve.. it has to come from existing b-smooths.. I know some of the conditions that must be met but others are still hidden and unless I find them using existing b-smooths is the only way. But very happy were its at now.

Note: Aaaalot of dead code that needs to be pruned. Will do it tomorrow.

Finally hit the 100 bit threshold. Will spent the day optimizing the code.. and after that I'll continue research and improvements (lots of small things to fix still like the code still doesnt work on negative polynomial values in psieve() etc... , I need to fix that so I can also properly optimize the sieve region there)

Update: So managed to filter out bad "a_mul" and "k" values with jacobi/kronecker symbols as well as possible. Any further filtering likely has to be done now with hensel. I do know that primes that divide the discriminant will have singular roots and they have unique behavior when lifting, and I know how to trivially lift them as non-singular by switching side.. so I think I can crack this problem wide open now. Hensel is the way forward. 

Going for a run... when I come back I'm going to further eliminate "k" values that cant possibly have a solution by lifting singular roots and seeing how they behave at higher exponents.... easy enough. Just got to keep hammering this now.... almost there.

Update: Quickly added some lifting logic when a solution is found. You can see how lifting singular roots work... and that it can be used as a tool to solve this.. 

I have multiple ideas, most straightforward is to use this to see if a solution exist outside the interval within some larger bound for solutions that survived marking but are not square after dividing out a.

Update: Let me try a more straightforward approach tonight. So with singular root lifting (for primes that divide a) we can generate a much sparser solution set within a bound... so rather then an interval based approach we use these solutions.. and we can filter them with those other primes if needed. Plus if we keep track of the coefficients residues for the other side around (ay^2+4Nk instead of ab^2-4Nk).. we know that a_mul divides out these coefficients.. so I'm sure I can come up with something there and replace that best_a_mul() function with it. I should also study some abstract algebra, rings and modules specifically.. I havnt yet found much time.. but I'm running out of money so I make a miracle happen now or its over for me.

#### To run from folder "Coefficient_Sieve" (For use with the paper):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run:  python3 run_qs.py -keysize 40 -base 50 -debug 1 -lin_size 10_000 -quad_size 100</br>

Just demonstrates the math from the paper using quadratics. For educational purposes. And rather then taking a square root over a large prime we can also just calculate the discriminant. But this demonstrates the interesting relation between these quadratics and the factors of N.

#### To run debug.py" (Prints the linear and quadratic coefficients to solve for 0 in the integers, for use with my paper):</br></br>

To run: python3 debug.py -keysize 12

This basically creates a system of quadratics. Solving them mod p is easy. But there is only one root solution (the factor of N) which solves the system for 0 for any mod p (aka solves it in the integers). Figuring out how to exactly do this quickly is still an ongoing area of research for me. And if a polynomial time algorithm for factorization exists, it is likely done by solving this system of quadratics. 

(PERFORMING HEATHEN CYBER RITUALS IN THE BELGIAN MOORLANDS TO DESTROY MY ENEMIES)

