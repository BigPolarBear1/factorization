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

#### To run from folder "NFS_Variant_WIP" (Implements more of my number theory and attempts to succeed with fewer smooths by using p-adic lifting):</br></br>

To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -key 4387 -base 30 -sbase 30 -debug 1 -lin_size 1_000 -quad_size 1_000   

Update: I got NFS_Variant_WIP working with the linear coefficient. Now I need to figure out how I can add the modulus to the linear coefficient and still succeed at the linear algebra step. This will allow superior control over the smooth bit length. After that we should be finished with this math project after nearly 3 years... I have one week to finish this as I may be arrested on the 16th of december since the americans are using the belgian judicial system to harass me and I have to go into the police office on the 16th.

Update: Oh yea. It is possible! If we add the modulus to the linear coefficient y0. Then we must simply reconstruct the quadratic for y1 to generate the same output. IT IS LIKE NFS! I knew it! I knew I was right! Expect disclosure this week. 

UPDATE: YES. EUREKA! THIS IS 100% IT! I KNOW IT! If we add the modulus to y0, then this means the derivate (aka y1) will be non-zero. But in order to take a derivative, we need to find the correct root, since we only know zx^2 is square after the linear algebra step, but we dont know the root... which we need if we have non-zero linear coefficient. And ofcourse, this root isnt garantueed to simply be the square root of zx^2. But with everything I got, NOW I CAN TAKE THE SQUARE ROOT OVER A FINITE FIELD! I GOT IT! HAH! If I ever find out the belgian government knew I was on the right track with my math and yet still kept harassing me with the police, I garantuee you, I will flee this shit country. They'll have to kill me to prevent that from happening. I'm not a fucking dog. Fucking nazis. I'm a polar bear, and this week, I'm going to devour this shit world.

UPDATE: I expanded the paper a little more. Now only one final thing is missing in the paper. Showing how to reconstruct the root if the derivative of the polynomials we multiply together is non-zero. And that's it.. after that I added a few lines of code to the PoC and we're done (atleast all the math will be done.. after that its just optimizing the PoC).
