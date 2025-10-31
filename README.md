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

I just uploaded polarbearalg_v24.py</br>
To run: python3 polarbearalg_v24.py -keysize 14</br>

I finally understand that quadratic coefficient!!!! 
So what we need to have is for zx^2 to be a square. And if it produces a square with zx^2-N, then we have a square relation. This quadratic coefficient can be used to compress these multiples of N down to just a single N in a sense. 

ok so zx^2 has to be square. So we find z values (quadratic coefficients) that are square. Then zx^2-n also has to be square. However, if z mod p^2 is square, then we can simply add p^2 to the quadratic coefficient. Until we generate a very small value.. and we are garantueed to have a smooth. THIS IS THE 2D SIEVING I WAS LOOKING FOR! THIS IS IT! It's not a perfect analogue of NFS, but its some hybrid. Holy fuck, I got to move fast on this, CIA assasination incoming for sure now (not that they stand a chance, I would kill them with my bear hands if they tried).

I GOT IT! EUREKA! FINALLY. HOLY FUCK. GOD DAMNIT. IT'S SO EASY AND SIMPLE. WHY DID I MISS THIS? WTF. OK TOMORROW... I WILL WORK HARD. I WILL SIMPLY MODIFY THE QS_Variant AND IMPLEMENT WHAT I LEARNED HERE. :)
