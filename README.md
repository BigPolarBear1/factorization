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

#### To run from folder "Improved_QS_Variant" (Tries to fruther reduce the bitlength of smooths by combining SIQS's approach of reducing bitlengths with our method of calculating moduli):</br></br>

To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -key 4387 -base 6 -debug 1 -lin_size 100_000 -quad_size 1</br></br>

If you use the above command you can generate exactly what I'm describing in chapter 8 of my paper.
If is now also generating the hashmap with quadratic coefficient to y0/y1 mappings.

Next I will show you how to pull small quadratic coefficients from this. This is the only thing missing in my code. Since it already generates values that are always smooth as long as we do modular division using factors from the factor base. I will make sure my work is finished before tuesday or atleast leave detailed instructions before tuesday as I'm expected to go into the police office with my lawyer on tuesday. And I may end up in jail, since the Americans have been weaponizing the Belgian judicial system against me relentlessly.

Note: To finish this, you use the residue map to figure out how often you need to do modular division to get to a small quadratic coefficient (modulo some large modulus). This works because using this residue map we can work modulo m rather than modulo n. I will try to implement the necessairy calculations before I have to go into the police office on tuesday. Once that is done.. thats it.. you would just need to optimize the code and you're good.... I likely wont have time to do much optimizations before tuesday, but i'm doing what I can.

Update: I'll add the exact math of how k in zx^2-n*k changes as we perform modular division and we end up with a wrap around due to z not being divisible in the integers, to the PoC and paper tonight. I already worked out the math, its not difficult. Then tomorrow I'll write the code to use that hashmap... and monday I'll prepare everything should I get arrested on tuesday. Its hilarious, these fucking americans know I'm going to win, so this is their only way to try and stop me. Desperate fuckers. On monday, I will send my work to every non-western contact I have with instructions on how to finish my work. Good luck with that. I know what you people are doing. Go to hell. Fuck belgium most of all. I disown my own country. Do not ever refer to me as belgian. I am not belgian. Being called belgian is a fucking insult. I despise my own country.

Update: added the math for that to the PoC... now let me add it to the paper aswell.
