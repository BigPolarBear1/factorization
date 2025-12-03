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
To run: python3 run_qs.py -keysize 80 -base 10_000 -sbase 1_000 -debug 1 -quad_size 10_000</br></br>

Update: I will upload an improved version later this week. The "overall" idea of the uploaded PoC is fine. But rather then finding a large square with a single sieve interval to reduce the size of the quadratic coefficient, we need a recursive approach instead. I.e if we define the "center" as the point where the smooth candidate is smallest for a given root, then we must simply find a square as close to that center as possible... divide the quadratic coefficient by the square and multiply the root with the square root... and then we repeat this. Ofcourse while this can be represented as recursion, I don't want to actually write a recursive function... as recursive functions kind of suck... so I just do a iterative approach where I track the state at each depth in a list (aka simulate recursion)... I really hate programming finnicky shit like this and I'm having a hard time focusing these last few days.. I'll try to finish it soon. This will work once its done.. its just a pain to code. Anyway.. I told myself.. I will keep working on this fulltime till the end of the year before focusing on making money again... so I still got time to finish this project.. 4 weeks is plenty to finish this code and finish the paper. 

Update: I'll go for a run. I've started coding these improvements... I'll try to upload either tonight or tomorrow. I've implemented some high level logic already.. all thats left is a few lines of code to push and pop the state to a list as we go down or up in depth. I'll finish it when I get back from running. I'll have to quicly write down a numerical example by hand to make sure I'm tracking the state correctly. 
