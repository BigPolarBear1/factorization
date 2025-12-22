Disclaimer: At no point did AI contribute anything to this research project. Copilot can't even calculate basic congruences without making mistakes. People who think AI can do novel math research are delusional.

#### To run debug.py" (Prints the linear and quadratic coefficients to solve for 0 in the integers, for use with my paper):</br></br>

To run: python3 debug.py -keysize 12

This basically creates a system of quadratics. Solving them mod p is easy. But there is only one root solution (the factor of N) which solves the system for 0 for any mod p (aka solves it in the integers). Figuring out how to exactly do this quickly is still an ongoing area of research for me. And if a polynomial time algorithm for factorization exists, it is likely done by solving this system of quadratics. Finding a polyomial time algorithm is my ultimate goal, as this would make progress toward solving p = np as well. 

#### To run from folder "QS_variant" (Standard SIQS with our number theory as backend):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 220 -base 20_000 -debug 1 -lin_size 10_000_000 -quad_size 1</br></br>

See below for an improved way of performing what this PoC does.. I'll delete this Proof of Concept once the PoC for Improved_QS_Variant matures a little more<br><br>

Note: With a large enough -base and lin_size this PoC will find smooths for 110 digits. Albeit very slowly, but this is a highly unoptimized cython PoC. However, to push beyond that into novel terroritory for Quadratic Sieve-based algorithms we need to use quadratic coefficients and p-adic lifting, and that is what the PoC below (Improved_QS_Variant) will be for. 

#### To run from folder "Improved_QS_Variant" (Building smooths from the ground up... see paper on clues to finish the code):</br></br>

To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 14 -base 20 -debug 1 -lin_size 100_000 -quad_size 1</br></br>

Coming soon...

Update: Took a break today. I'm going to finish this tomorrow though. I know exactly how to do it.... just in time for christmas. Its always funny seeing threat intel folks cry about having to work during holidays and shit. Thats why I used to drop 0days, bc I hated those people and seeing their stupid nerd tears was funny as hell. Also they are a bunch of self-righteous losers. They arn't hackers. Real hackers dont work for big tech, real hackers peddle their goods on the dark web so their exploits get used to hack the FBI hahahaha.

I really despise belgium. This country has always been a dead end for me.

Fuck it. Time to covers some ground today. I need to get out of this horrible country. My life in belgium has always been shit. European people in general. They are just snobs. Life in european society is very much like highschool, where if you dont fit in, you are just relentlessly bullied and excluded. Its just like that. I honestly have 0 respect for europe, I even despise europe. Sadly there isnt many good alternatives in the world anymore. Everything is just going to shit. Time to finish my work, maybe once they realize over in China what I've got... I can find something over there... just go somewhere I'm respected for my work and capabilities. Surely isnt in the west. I deeply despise western infosec. Western infosec isnt about skills. Its about who you know. Which is exactly why countries like China are winning and the west is losing.
