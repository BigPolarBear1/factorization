DISCLAIMER: AT NO POINT IN MY 2.5 YEARS OF RESEARCH DID I USE AI OR HAS AI IN ANY WAY HELPED WITH MY RESEARCH. I'M PUTTING THIS HERE BECAUSE I KNOW PEOPLE WILL ACCUSE ME OF IT. 

Note: I am antifa, leader of the cipherpunk, fuck the FBI department. 

Note2: It's basically done now, the NFS variant is nearing completion soon, nothing novel research wise has to be figured out now. I just need to refactor a few things and I'm done.
I have been unemployed for 2 years. I am looking for work. I would like some way where I can have enough time to attack ECC and DLP next. In addition, for a payment, I can also show people how I appraoched factorization from start to finish. Because the method of how I conduct my research and consistently have success is just as valueable as the research itself. And I follow the same patterns over and over again, no matter what subject matter I attack. I don't have a very high IQ and I'm not especially bright or good at doing numbers in my head. But I have insights into research that help me succeed over and over again. I'm willing to teach people how to do it, as long as it enables me to have the finances to continue doing my research aside from that. Contact: big_polar_bear1@proton.me.

#### To run from folder "QS_variant" (Standard SIQS with our number theory as backend):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 200 -base 6000 -debug 1 -lin_size 100_000 -quad_size 1</br></br>
This version is somewhat done, but it doesn't achieve much of an advantage, as the main strength of my work will be the NFS_Variant, which performs NFS with reducible quadratic polynomials<br><br>
#### To run from folder "NFS_variant" (Number Field Sieve with our number theory as backend using reducible quadratic polynomials mod m):</br></br>
To run: python3 polarbearalg_debug.py -key 4387 

I just uploaded polarbearalg_v15.py</br>
To run: python3 polarbearalg_v15.py -keysize 14</br>

Update: Ergh. v15 has some strange behavior. It's able to consistently yield the factorization with two squares that are not congruent mod N. And very abstractly I have a sense for what is happening. The interesting thing is, this behavior also persist when I modulo reduce the entire discriminant. There is something happening there. 

And this is also something I noticed earlier.
So if we have a square mod M, it's roots mod M have a significant chance of yielding the factorization when taking the gcd. And I understand how that may happen, since we calculate these things from possible factors mod p. (see first chapters in paper). And I do remember these becoming more sparse as N becomes bigger.... which prevented me from just bruteforce taking the GCD on root combinations mod M. But I may be able to achieve it with linear algebra... hmm. I got to dig in tomorrow. There is something here I can use.

I really need to get to the bottom of this tomorrow.. I have this feeling that figuring this out, will also show me how to correctly implement this linear algebra step.  Lets see. Anyway... I need to sleep i guess. Days are too short.
