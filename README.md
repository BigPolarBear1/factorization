DISCLAIMER: AT NO POINT IN MY 2.5 YEARS OF RESEARCH DID I USE AI OR HAS AI IN ANY WAY HELPED WITH MY RESEARCH. I'M PUTTING THIS HERE BECAUSE I KNOW PEOPLE WILL ACCUSE ME OF IT. 

Note: I am antifa, leader of the cipherpunk, fuck the FBI department. 

Note2: It's basically done now, the NFS variant is nearing completion soon, nothing novel research wise has to be figured out now. I just need to refactor a few things and I'm done.
I have been unemployed for 2 years. I am looking for work. I would like some way where I can have enough time to attack ECC and DLP next. In addition, for a payment, I can also show people how I appraoched factorization from start to finish. Because the method of how I conduct my research and consistently have success is just as valueable as the research itself. And I follow the same patterns over and over again, no matter what subject matter I attack. I don't have a very high IQ and I'm not especially bright or good at doing numbers in my head. But I have insights into research that help me succeed over and over again. I'm willing to teach people how to do it, as long as it enables me to have the finances to continue doing my research aside from that. Contact: big_polar_bear1@proton.me.

#### To run from folder "QS_variant" (Standard SIQS with our number theory as backend):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 200 -base 6000 -debug 1 -lin_size 100_000 -quad_size 1</br></br>
This version is somewhat done, but it doesn't achieve much of an advantage, as the main strength of my work will be the NFS_Variant, which performs NFS with reducible quadratic polynomials<br><br>
#### To run from folder "NFS_variant" (Number Field Sieve with our number theory as backend using reducible quadratic polynomials mod m):</br></br>

I just uploaded polarbearalg_v28.py</br>
To run: python3 polarbearalg_v28.py -keysize 20</br>

I also updated the paper a bit more today. So v28, if zx^2-n is a square, and x^2-nz completely factors over our square and z, then we can set y0 to be the square and the quadratic coefficient for this new square will very likely be a multiple of the same square, hence allowing us to find y1 and get our square relation mod N. This succeeds at a fairly consistent rate if these conditions are met. However, this will succeed 100% of the time if our initial z value is a square... but we do not want to restrict ourselves to having a square there as it defeats the advantage of having a linear modifier there (meaning we can generate very small smooth candidates... much smaller then any QS variant is able to). 

It demonstrates the core ideas and math already. I need to do proper sieving now and add linear algebra.. this should definitely be done now very shortly... people would have known I'm doing everything right... I hope you all go to hell for treating me like this despite knowing I'm doing everything right. 

