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

Update: A realization just hit me. Ergh. 

So with the paper, we can precalculate zx^2+N = 0 mod p.
We can create a 2d sieve interval where the rows are the primes and the columns the roots (and we extend all primes up to some sieve interval size).

Now the thing is.. we can do vector addition... which is very fast.. add all rows together... and see if we have a column that appears in a lot of primes (we can use log values and add them together). 
Then do those calculations to shift all my solutions to the next quadratic coefficient and repeat... all just fast vector operations.

Wait, what? I literally cant think of a single reason why this wouldn't work. Oh my god. Let me rush out a PoC by Christmas. If that fails miserably.. I'll just continue what I was doing these last few days. But its so simple and elegant that I have to drop everything and atleast try it. Writing the code is super simple since Ive worked out all the math for a long time now..

God damnit. How could I not have seen this earlier? Its so simple, yet elegant. And I have all the math to do it now. God damnit! Happy christmas assholes. Factorization is coming!

OMG. WHY DIDNT I SEE THIS EARLIER. I can literally set it up, so everything can just be done with vector operations. EVERYTHING. FUCK!!!! 

https://youtu.be/2zpCOYkdvTQ

TWO DAYS. THAT IS ALL I NEED. Happy holidays assholes. 

Update: Remember that scene from lord of the rings, when the orcs finally crawl out of the mud pits and start beating the war drums to march on their enemies? Thats me today after 4 hours of sleep. I may just try to use CUDA for python this time around, writing SIMD optimized code kinda sucks. But let me just start by writing a python script with numpy arrays first. Anyway... first I have to go to my lawyer, because of the US gov harassing me. Thanks for the lawyer fees btw. My price to not go to China after I break factorization is now 2 billion for me and my former manager, plus a huge castle somewhere super remote so I dont ever have to interact with humans ever again. It keeps going up the more I'm harassed and my former manager doesn't see justice for what microsoft did to him.
