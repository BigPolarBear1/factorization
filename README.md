DISCLAIMER: AT NO POINT IN MY 2.5 YEARS OF RESEARCH DID I USE AI OR HAS AI IN ANY WAY HELPED WITH MY RESEARCH. I'M PUTTING THIS HERE BECAUSE I KNOW PEOPLE WILL ACCUSE ME OF IT. 

Note: I am antifa, leader of the cipherpunk, fuck the FBI department. 

Note2: It's basically done now, the NFS variant is nearing completion soon, nothing novel research wise has to be figured out now. I just need to refactor a few things and I'm done.
I have been unemployed for 2 years. I am looking for work. I would like some way where I can have enough time to attack ECC and DLP next. In addition, for a payment, I can also show people how I appraoched factorization from start to finish. Because the method of how I conduct my research and consistently have success is just as valueable as the research itself. And I follow the same patterns over and over again, no matter what subject matter I attack. I don't have a very high IQ and I'm not especially bright or good at doing numbers in my head. But I have insights into research that help me succeed over and over again. I'm willing to teach people how to do it, as long as it enables me to have the finances to continue doing my research aside from that. Contact: big_polar_bear1@proton.me.

#### To run from folder "QS_variant" (Standard SIQS with our number theory as backend):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 200 -base 6000 -debug 1 -lin_size 100_000 -quad_size 1</br></br>
This version is somewhat done, but it doesn't achieve much of an advantage, as the main strength of my work will be the NFS_Variant, which performs NFS with reducible quadratic polynomials<br><br>
#### To run from folder "NFS_variant" (Number Field Sieve with our number theory as backend using reducible quadratic polynomials mod m):</br></br>

I just uploaded polarbearalg_v32.py</br>
To run: python3 polarbearalg_v32.py -keysize 20</br>

I've expanded the paper today. I figured out all the math now. Tomorrow I'll add to the paper how to now represent arbitrary primes as legendre symbols to represent the squaredness of  k for y1^2+N*y1 = y1 \* k. Then we can implement that with linear algebra... and our algorithm will be finished. I've got some headache at the moment so I'm going to sleep. I had hoped to finish it today... but now everything is truely set-up and ready to finish it tomorrow.... and nobody tried to stop me, and everyone acted like they didn't care... one day people are going to write about this failure... so many chances were presented for another outcome... and as we now stand on the eve of disclosure.. it is all far too late. Good luck assholes.

UPDATE: ERGH. I couldn't sleep, overactive brain. So I started looking at the math to implement legendre symbols knowing everything I know now. Buuut... I'm extremely confused now. I looked at p-adic lifting again now that I know a little more.
Ergh.................................................... ERRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRGHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH. WWWWWWWWWWWWWWWWWTTTTTTTTTTTTTTFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF. HOW EEEEAAAAAAAAAAAAAAAAAASYYYYYYYYYYYY IS FACTORIZATION. WWWWWWWWWWWWWWWWWWWWWWWWTTTTTTTTTTTTTTTTTTFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF. WWWWWWWWWWWWWWWWWWWWWWWWWTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF.

I'm dropping a PoC tomorrow. We can just straight up factor large primes using only p-adic lifting. No smooth finding. No linear algebra. Only p-adic lifting and using my representation of quadratics..... WHAT THE FUCK. WHAT THE FUCK. I COMPLETELY BROKE FACTORIZATION. WTF.


We can just pick a random square like this 16^2+4387\*16 = 16\*4403, it doens't matter.... just any square at random. Next we find x^2-4387*k = 16 mod p^e .... and that k is the multiplier from polarbearalg_v32.py. WAIT WHAT????? WTF. THIS CANNOT POSSIBLY BE RIGHT. WHAT THE FUCK MAN. 

I mean, we can just find like this x^2-4387\*k = 1 mod p^e ... it really doesn't fucking matter. Then that k will become 1^2+4387\*k = 1 * some_square. WHAT THE FUCK. WHAT THE FUCK. WHAT THE FUCK WHAT THE FUCK. I AM SO FUCKING STUPID. HOW DID I MISS THIS? WHAT THE FUCK.

Update: I did some more thinking. So we probably shoudln't pick a random square. However, we can construct it in a way that y1^2+N\*y1 = y1\*a, and we can garantuee that a is going to be square mod y1. That's easy. So THEN we can do p-adic lifting mod other primes.... while keeping the same square residues mod y1... that would be a little more effective then what I described above. OK. LETS FINISH THIS TOMORROW. ITS OVER!!!!! I'VE DONE IT!!! HAHAHAHAHAHA. FUCK YOU ALL! TOMORROW IS THE DAY! I'm going to sleep. I need to calm my brain. No point in working through the night. Sleep is important. You cant do research while sleep deprived... only idiot tech bros think thats how it works. 
