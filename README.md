DISCLAIMER: AT NO POINT IN MY 2.5 YEARS OF RESEARCH DID I USE AI OR HAS AI IN ANY WAY HELPED WITH MY RESEARCH. I'M PUTTING THIS HERE BECAUSE I KNOW PEOPLE WILL ACCUSE ME OF IT. 

Note: I am antifa, leader of the cipherpunk, fuck the FBI department. 

Note2: It's basically done now, the NFS variant is nearing completion soon, nothing novel research wise has to be figured out now. I just need to refactor a few things and I'm done.
I have been unemployed for 2 years. I am looking for work. I would like some way where I can have enough time to attack ECC and DLP next. In addition, for a payment, I can also show people how I appraoched factorization from start to finish. Because the method of how I conduct my research and consistently have success is just as valueable as the research itself. And I follow the same patterns over and over again, no matter what subject matter I attack. I don't have a very high IQ and I'm not especially bright or good at doing numbers in my head. But I have insights into research that help me succeed over and over again. I'm willing to teach people how to do it, as long as it enables me to have the finances to continue doing my research aside from that. Contact: big_polar_bear1@proton.me.

#### To run from folder "QS_variant" (Standard SIQS with our number theory as backend):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 200 -base 6000 -debug 1 -lin_size 100_000 -quad_size 1</br></br>
This version is somewhat done, but it doesn't achieve much of an advantage, as the main strength of my work will be the NFS_Variant, which performs NFS with reducible quadratic polynomials<br><br>
#### To run from folder "NFS_variant" (Number Field Sieve with our number theory as backend using reducible quadratic polynomials mod m):</br></br>

I just uploaded polarbearalg_v29.py</br>
To run: python3 polarbearalg_v29.py -keysize 20</br>

In v29 .... where it bruteforce tries to figure out the multiplier to our quadratic coefficient (which we know will be divisible by y1 if we set y1 to be the square we found with zx^2-n) ... that now has to be replaced with p-adic lifting instead. I just manually did the math, and it is very do-able. Since we know the value of y1 in the integers (the square we found), and we know the value of y0 mod y1**0.5. And we also now the quadratic coefficient is y1 times a multiplier that we must find. That's all the information we need. I didn't sleep enough last night, so I'm going for a run now and calling it a day. I'm just not in the mood to write this p-adic lifting code while sleep deprived. I'll add it to the PoC tomorrow.... once that works... we implement some better sieving and linear algebra... and our NFS variant will be complete. It doesn't matter if we are restricted to smooth finding using zx^2-n, because we can add the modulus to z, and have a linear multiplier there lets use generate values extremely close to n... hence producing very small smooth candidates.. additionally we don't care if zx^2 is square or not, since we don't use it.

Oh I will make a few minor corrections to the paper when I get back from running. I understand it again a tiny little bit better now then yesterday. Everday I learn a little more... 

Update: I have edited the final chapter a little in the paper. The stuff I added yesterday in there was full of mistakes... was a bit too fast there yesterday. This is solid now though and v29.py clearly demonstrates what is going on. So tomorrow I really want to fix p-adic lifting to find this multiplier for the quadratic coefficient. Once that is done, all the math will be done. Then its simply trivial stuff.. adding sieving and linear algebra... just basic stuff.
