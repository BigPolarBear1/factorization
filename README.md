Update: It would appear the americans are continueing to legally harass me because pete hegseth is a fucking loser and coward who has to hide behind the law. So much for all his pathetic bravado. Fucking pussy. If anyone wants to help cover legal costs email: big_polar_bear1@proton.me ... I am also considering fleeing europe, if anyone can grant me asylum also reach out, thanks.

DISCLAIMER: AT NO POINT IN MY 2.5 YEARS OF RESEARCH DID I USE AI OR HAS AI IN ANY WAY HELPED WITH MY RESEARCH. I'M PUTTING THIS HERE BECAUSE I KNOW PEOPLE WILL ACCUSE ME OF IT. EITHER THAT OR SOMEONE WILL CLAIM HE PUBLISHED BEFORE I DID WHILE PLAGIARIZING WORK FROM SOME OF MY EARLIER RELEASED PAPERS. THIS WORK HAS BEEN COMPLETED FOR 99% FOR A LONG TIME NOW, IT JUST TOOK ME AGES TO MAKE THESE FINAL CONNECTIONS BECAUSE I HAVE NO MATH EDUCATION, I HAVE NOBODY TO ASK FOR HELP. I HAVE TO GRIND IT OUT THE HARD WAY. THERE IS NO WAY PEOPLE WILL LET A TRANS PERSON WIN. EITHER WAY. I KNOW MY WORK WILL BE DONE THIS WEEK, SO GOOD LUCK

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
To run: python3 run_qs.py -keysize 14 -base 100 -sbase 11 -debug 1 -lin_size 100 -quad_size 1000

Alright nearly there. To achieve superior control over the size of smooth candidates we need to add or subtract the modulus to the linear coefficient. Right now the PoC will fail if after the linear algebra step the k != 0 (aka the modulus has been added to the linear coefficient). To still be able to take the GCD after we add the modulus to the linear coefficient, the first step we must achieve is by using Chinese remainder and lifting of roots and coefficients, find roots and coefficients that produce the same value as the product of smooths. Once we have that we'll be able to take the GCD, since we'll also know the coefficients for both y0 and y1. And that's that. A lot simpler then NFS. I actually think this is the correct way to do it. NFS actually overcomplicates this entire setup by using number fields. :) Fixed PoC will be released this week and full paper too. 


Update: oh yea... actually, the only thing that matters is that we can create this x^2-n*z, this is why not adding the modulus to the linear coefficient works, because if the derivative is 0 then that means there will be a quadratic generating the same result without linear coefficient.. so when we do add the modulus to the linear coefficient, we must meet the conditions where we can manipulate the root to get rid of that linear coefficient. I figured out these mechanics before. Let me add a new chapter to the paper before chapter 8 to build up to what I'm currently trying to do.
