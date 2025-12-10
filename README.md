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
To run: python3 run_qs.py -keysize 100 -base 1000 -debug 1 -lin_size 100_000 -quad_size 1</br></br>

I've started refactoring. Still lots of work to be done.

There is a bunch of bugs in the uploaded PoC due to quickly putting it together copy pasting from some of my olds PoCs.
What this PoC does:

1. First it generates a large quadratic coefficient made up from factors of the factor base (this is basically what SIQS does).
2. We generate moduli for the sieve interval step size, such that the smooth candidate is divisible by the modulus.

All still very easy and basic. However, remember when I was earlier experimenting with large squares to reduce the size of smooth candidates? So since the quadratic coefficient adds an offset to the parabola and I now actually understand the mechanics at play here now... I'm feeling confident I can figure it out now. We must use the quadratic coefficient to shift a very large square into our sieve_interval range. I know this is do-able. 

Tomorrow I will thus also add p-adic lifting again, and then I will begin working out the math to use this quadratic coefficient to get a large square in range of the sieve_interval (with large square I mean, a square almost as big as N so the smooth_candidate is extremely small.. a technique like that would the only way to actually beat existing factorization methods)

Update: I'll figure it out tomorrow. This setup makes it a lot easier... now we dont need to worry about that quadratic coefficient like we had to before. So tomorrow I need to re-order some things. The outline would look like this:

1. Generate a large square modulus. Close to the size of N
2. Generate the root for this modulus for quadratic coefficient = 1
3. Multiply the quadratic coefficient with factors from the factor base (we know the mechanics of this now, easy)
4. Figure out how to gain fine grain control over the size of the root so we can get it setup to generate very small smooth candidates when divided by the very large square.

I know this is possible. It cant be so hard. I'm sure equipped with everything I figured out today, I can finally figure out this final step tomorrow.

UPDATE: NO NO NO. WE DONT GENERATE A SINGLE LARGE SQUARE. WE GENERATE SMALL SQUARES!!!!! THEN WE DO TRIAL FACTORIZATION ON THE ROOT, IF THE INVERSE OF A FACTOR ALSO FACTORS OVER THE FACTOR BASE, WE CAN MULTIPLY THE QUADRATIC COEFFICIENT AND BE GARANTUEED A SIZE REDUCTION OF THE ROOT. THEN WE JUST NEED TO SUCCEED AT THIS FOR A COUPLE OF SQUARES, MULTIPLY THEM TOGETHER AND BOOM BADABOOM DONE. LOL! EZ!! This ends tomorrow assholes. Prepares yourselves (good luck patching in a single day, guess you cant really prepare yourselves, hahaha, fuck you losers).
