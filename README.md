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

Update: Let me add a full numerical example to the paper. And I'll finish the code in the coming days after that. Even if I dont complete my work by christmas, I can still ruin people's new year celebrations by breaking the internet. Everyone will be just as miserable as I am. Hahaha. 

Update: One sec... I think the proper way to do it is to restrict to a single quadratic coefficient. To keep that quadratic coefficient small. Let me edit some stuff in the paper.

Update: Added a proper example to the paper. Let me add to the PoC some lifting for when we check the hashmap... that should further filter out results we cant combine...

Update: I'll go for a run. Actually you know what I'm going to do tomorrow. Just sieve at quadratic coefficient z = 1. So iterate x and y for 1\*x^2-yx+N .. make sure every factor has a valid mapping in the precalculated hashmap (eventually we should sieve using the hashmap directly)  .... if it has add it as smooth... perform gaussian elimination over GF(2) ... and then we try to Chinese remainder these quadratics together and with some luck we'll have a square when taking the discriminant and we take the GCD... that should work often enough... and if not.. I'm sure I can figure it out. And then I just take it from there.. keep improving it. I know I'm almost there... and I know even if people know I'm right, they would never tell me unless I set a factorization record... because western infosec is a bunch of losers. Could have just paid me years ago when I started this and prevented everything. But my entire career has been dealing with this kind of bullshit.. had it been anyone else doing this, they would have been treated a lot better. People are instead literally gaslighting me into believing I'm mentally ill. I'm not. I'm so fucking pissed at everything people have done, nothing can salvage this. Unless they pay me and my former manager somewhere close to a billion.. aside from that.. I will end up in China or some non-western country after I set a factorization record, because I will never forgive this.
