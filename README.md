DISCLAIMER: AT NO POINT IN MY 2.5 YEARS OF RESEARCH DID I USE AI OR HAS AI IN ANY WAY HELPED WITH MY RESEARCH. I'M PUTTING THIS HERE BECAUSE I KNOW PEOPLE WILL ACCUSE ME OF IT. 

Note: I am antifa, leader of the cipherpunk, fuck the FBI department. 

Note2: It's basically done now, the NFS variant is nearing completion soon, nothing novel research wise has to be figured out now. I just need to refactor a few things and I'm done.
I have been unemployed for 2 years. I am looking for work. I would like some way where I can have enough time to attack ECC and DLP next. In addition, for a payment, I can also show people how I appraoched factorization from start to finish. Because the method of how I conduct my research and consistently have success is just as valueable as the research itself. And I follow the same patterns over and over again, no matter what subject matter I attack. I don't have a very high IQ and I'm not especially bright or good at doing numbers in my head. But I have insights into research that help me succeed over and over again. I'm willing to teach people how to do it, as long as it enables me to have the finances to continue doing my research aside from that. Contact: big_polar_bear1@proton.me.

#### To run from folder "QS_variant" (Standard SIQS with our number theory as backend):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 200 -base 6000 -debug 1 -lin_size 100_000 -quad_size 1</br></br>
This version is somewhat done, but it doesn't achieve much of an advantage, as the main strength of my work will be the NFS_Variant, which performs NFS with reducible quadratic polynomials<br><br>
#### To run from folder "NFS_variant" (Number Field Sieve with our number theory as backend using reducible quadratic polynomials mod m):</br></br>

I just uploaded polarbearalg_v31.py</br>
To run: python3 polarbearalg_v31.py -keysize 20</br>


I've added some legendre symbols to that loop that tries to find the correct multiplier for the quadratic coefficient. 
And actually thinking about this, I need to get these legendre symbols working in my linear algebra step, similar to what NFS does with the quadratic character base.
I was thinking before that I should do p-adic lifting here to find that multiplier... but since the correct solution can be represented by legendre symbols, I should probably just implement that in my linear algebra step... then we don't really need p-adic lifting to find that correct multiplier. So I guess the next step now is to calculate some examples by hand to figure out the exact implementation details of the linear algebra step and how that would work. I think v31 is good now, it shows all the math, and also how legendre symbols come into play. Now it's time to get serious and start building a proper algorithm.

