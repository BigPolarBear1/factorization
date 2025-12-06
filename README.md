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
To run: python3 run_qs.py -keysize 100 -base 500 -sbase 500 -debug 1 -lin_size 10_000 -quad_size 5_000</br></br>

Update: Alright. I added some interval code. </br></br>

To do:</br></br>

1. The calculations in sieve() need to be sped up many many times more so that we can us a much bigger quad_size parameter.</br>
2. Currently it only uses -sbase... which is the small factor base.. but we should also use -base, the large factor base to saturate the intervals with large squares.</br></br>

What this code does:</br></br>

It will roughly generate smooths that are of keysize/2. Normally with normal SIQS you cannot do this (unless you use a very small sieve interval). We achieve this by using multiple quadratic coefficients.</br>
You want to keep -lin_size relatively small and really increase quad_size, but I will need to optimize the calculations in sieve() first. </br>
Additionally, we really really need to saturate our sieve intervals with large squares. Since they wont increase the required amount of smooths and they will help to chip away at the bit length.</br></br>

Update: I've also quickly added lifting for the primes in the -sbase</br>
Update: Added large prime variation and cleaned up dead code</br>

NOTE: The americans are using the courts to move against me. I may get formally arrested on Monday. If this happens and my github repo dissappears.. then you can put 2 and 2 together...  my situation is very dire here in Europe. I've faced persistent unemployment and legal harassment under pressure of the americans. They know I will succeed one day and that I am getting closer as time passes. They are trying everything they can to break me mentally. 
