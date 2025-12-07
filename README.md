Update: It would appear the americans are continueing to legally harass me because pete hegseth is a fucking loser and coward who has to hide behind the law. So much for all his pathetic bravado. Fucking pussy. If anyone wants to help cover legal costs email: big_polar_bear1@proton.me ... I am also considering fleeing europe, if anyone can grant me asylum also reach out, thanks.

DISCLAIMER: AT NO POINT IN MY 2.5 YEARS OF RESEARCH DID I USE AI OR HAS AI IN ANY WAY HELPED WITH MY RESEARCH. I'M PUTTING THIS HERE BECAUSE I KNOW PEOPLE WILL ACCUSE ME OF IT. 

Note: I am antifa, leader of the cipherpunk, fuck the FBI department. 

Note2: Once Improved_QS_Variant is finished, I will rewrite most portions of the paper completely from start. I keep finding myself "maturing" out of these papers.. because I keep learning and advancing with my math skills. There's so much stuff in that paper that I just want to completely rewrite because even to me it looks way too amateuristic now.

Note3: I am utterly broke so if anyone wants to make donations to keep this research going a little longer: big_polar_bear1@proton.me email me... thanks.

Note4: If what I'm trying to do currently doesn't work, I will completely abandon the quadratic sieve type of approach. Because using the number theory I worked out, another approach would be is to directly find a factor of N. Since using the factors of N solves the quadratic for 0 for any mod m when the correct linear and quadratic coefficients are used. So this will be the next avenue I will explore. 

#### To run debug.py" (Prints the linear and quadratic coefficients to solve for 0 in the integers, for use with my paper):</br></br>

To run: python3 debug.py -keysize 12

This basically creates a system of quadratics. Solving them mod p is easy. But there is only one root solution (the factor of N) which solves the system for 0 for any mod p (aka solves it in the integers). Figuring out how to exactly do this quickly is still an ongoing area of research for me. And if a polynomial time algorithm for factorization exists, it is likely done by solving this system of quadratics. Finding a polyomial time algorithm is my ultimate goal, as this would make progress toward solving p = np as well. 

#### To run from folder "QS_variant" (Standard SIQS with our number theory as backend):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 220 -base 20_000 -debug 1 -lin_size 10_000_000 -quad_size 1</br></br>

See below for an improved way of performing what this PoC does.. I'll delete this Proof of Concept once the PoC for Improved_QS_Variant matures a little more<br><br>

Note: With a large enough -base and lin_size this PoC will find smooths for 110 digits. Albeit very slowly, but this is a highly unoptimized cython PoC. However, to push beyond that into novel terroritory for Quadratic Sieve-based algorithms we need to use quadratic coefficients and p-adic lifting, and that is what the PoC below (Improved_QS_Variant) will be for. 

#### To run from folder "Improved_QS_Variant" (Implements more of my number theory and attempts to succeed with fewer smooths by using p-adic lifting):</br></br>

To build: python3 setup.py build_ext --inplace</br>
To run:python3 run_qs.py -keysize 120 -base 1000 -sbase 500 -debug 1 -lin_size 10_000 -quad_size 10_000   (note: Takes about 90 seconds for 120 bit.. its bottlenecking in sieve().. I need to improve that next)</br></br>

Note: I will remove this folder eventually.. I used this code as a start to construct NFS_Variant_WIP

#### To run from folder "NFS_Variant_WIP" (Implements more of my number theory and attempts to succeed with fewer smooths by using p-adic lifting):</br></br>

To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -key 4387 -base 30 -sbase 30 -debug 1 -lin_size 1_000 -quad_size 1_000   

This PoC is simply a slight modification of Improved_QS_Variant where we add the modulus to y<sub>1</sub>. This allows us to have better control over the size of the generated smooth (although in the code we should probably subtract the modulus instead so we shrink the smooths if they are already positive.. but I first want to get it to work like this before I worry about that).

In extract_factors() we now have to figure out the math to adjust this result. I think a similar setup to NFS using square roots over finite fields should now be possible. We can also calculate y<sub>0</sub> by taking the derivative. 

Note: I may get arrest on Monday because the americans have been pressuring the belgian police to harass me. If I go dark after monday... you know what happened. And I promise you, they will not get me alive. These americans they know I'm closing on with my work... they are frantically doing anything they can to stop me, so a transgender person doesn't win. It's nearly there now.... just figure out how you can use zx^2+y1\*x-n instead of zx^2-n like the uploaded PoC does... how to manipulate and even lift coefficients is all described in the paper. And QS_Variant also has some coefficient related code. I would prefer that I finish this math project myself... but it appears like people are desperately trying to stop me.




