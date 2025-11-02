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

I just uploaded polarbearalg_v27.py</br>
To run: python3 polarbearalg_v27.py -keysize 14</br>

Wasted the day, bah. Anyway, I gained a little more understanding. So if zx^2-n is square. Then we can set y0 to be the square, and we will also know that the z value (quadratic coefficient) for y0 will also be divisible by that same square. Which is a very useful property. Because if we know atleast some of the factors of our quadratic coefficient and we are given y0, then we can use p-adic lifting in a different primefield to figure out y1. And if the original z value is square (the one from zx^2-n that yielded the square), then we also know that y1 will be 0 modulo the square root of y0. However, this only holds true if the original z is square. Either way, it doesn't really matter, because we can figure out that residue using p-adic lifting, since having y0 and some factors of the quadratic coefficient is all we really need. I'll try to get that p-adic lifting stuff working early next week. 
