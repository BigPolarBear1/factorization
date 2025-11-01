DISCLAIMER: AT NO POINT IN MY 2.5 YEARS OF RESEARCH DID I USE AI OR HAS AI IN ANY WAY HELPED WITH MY RESEARCH. I'M PUTTING THIS HERE BECAUSE I KNOW PEOPLE WILL ACCUSE ME OF IT. 

Note: I am antifa, leader of the cipherpunk, fuck the FBI department. 

Note2: It's basically done now, the NFS variant is nearing completion soon, nothing novel research wise has to be figured out now. I just need to refactor a few things and I'm done.
I have been unemployed for 2 years. I am looking for work. I would like some way where I can have enough time to attack ECC and DLP next. In addition, for a payment, I can also show people how I appraoched factorization from start to finish. Because the method of how I conduct my research and consistently have success is just as valueable as the research itself. And I follow the same patterns over and over again, no matter what subject matter I attack. I don't have a very high IQ and I'm not especially bright or good at doing numbers in my head. But I have insights into research that help me succeed over and over again. I'm willing to teach people how to do it, as long as it enables me to have the finances to continue doing my research aside from that. Contact: big_polar_bear1@proton.me.

#### To run from folder "QS_variant" (Standard SIQS with our number theory as backend):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 200 -base 6000 -debug 1 -lin_size 100_000 -quad_size 1</br></br>
This version is somewhat done, but it doesn't achieve much of an advantage, as the main strength of my work will be the NFS_Variant, which performs NFS with reducible quadratic polynomials<br><br>
#### To run from folder "NFS_variant" (Number Field Sieve with our number theory as backend using reducible quadratic polynomials mod m):</br></br>
To run: python3 polarbearalg_debug.py -key 4387 (for use with the final chapter of the paper... but ignore this and use PoC below until I refactor the paper one final time)

I just uploaded polarbearalg_v25.py</br>
To run: python3 polarbearalg_v25.py -keysize 14</br>

This is just for the math I added to the paper today (v25). It shows how if we find a square with our polynomial while using square moduli, then that also tells us y1 will be a multiple of that square root. And y0 will be the square (or N minus the square.. both are congruent). This math works perfectly for the square found at quadratic coefficient 1 ... I need to add a few more calculations to generalize this.. let me work out the math for this. But this is much better .... this is now well and truely becoming a NFS variant. 

Update: wait I got it, this also tells us about the divisors of the quadratic coefficient... ergh... soooooooooooooooooo close. Let me take a short break first. My goal for this weekend is to figure out the exact math to make this work for the square relation found at any quadratic coefficient, because the uploaded PoC only finds the one from quadratic coefficient 1. But I know that we can also reverse engineer what the factors of the quadratic coefficient should be with this approach.. I just did the math for that. I just need to think clearly for a few hours and figure out exactly how I'm going to pull all of this together into something coherent. 
