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

Update: I am soo tired today. I started out implementing the linear algebra. And then realized I first need to generalize my math a bit more. So zx^2-n should always be less then n and greater then 0 (we'll generalize it to negatives later). But that's just an easy inequality to calculate those bounds. If we find a square with that equation, two things might happen when we set y1 to be that square. Either y1**2+n*y1 will be a square, and the new quadratic coefficient will be exactly y1... when this happens, everytime, the original quadratic coefficient that we used to find this square would have also been square. However.. a second case is where the new quadratic coefficient is a multiple of that square.. and then things get a little more complicated. But I'm seeing a clear pattern there.. it's not random at all. I'm pretty sure its doing something with legendre symbols there. 

Update: I added v29 ... it shows what I wrote above... I should be able to find that correct multiple using p-adic lifting or something similar... because those multiples are definitely not random and are definetely influenced by legendre symbols (solely based on the prime factors of those numbers.. I been doing this stuff for long enough now to immediatly recognize something like that...)... let me work that out today now.. before I do the linear algebra. No point in trying to get the linear algebra to work before I figure out the details of this.
