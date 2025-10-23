Note: I am antifa, leader of the cipherpunk, fuck the FBI department. 

Note2: It's basically done now, the NFS variant is nearing completion soon, nothing novel research wise has to be figured out now. I just need to refactor a few things and I'm done.
I have been unemployed for 2 years. I am looking for work. I would like some way where I can have enough time to attack ECC and DLP next. In addition, for a payment, I can also show people how I appraoched factorization from start to finish. Because the method of how I conduct my research and consistently have success is just as valueable as the research itself. And I follow the same patterns over and over again, no matter what subject matter I attack. I don't have a very high IQ and I'm not especially bright or good at doing numbers in my head. But I have insights into research that help me succeed over and over again. I'm willing to teach people how to do it, as long as it enables me to have the finances to continue doing my research aside from that. Contact: big_polar_bear1@proton.me.

#### To run from folder "QS_variant" (Standard SIQS with our number theory as backend):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 200 -base 6000 -debug 1 -lin_size 100_000 -quad_size 1</br></br>
This version is somewhat done, but it doesn't achieve much of an advantage, as the main strength of my work will be the NFS_Variant, which performs NFS with reducible quadratic polynomials<br><br>
#### To run from folder "NFS_variant" (Number Field Sieve with our number theory as backend using reducible quadratic polynomials mod m):</br></br>
To run: python3 polarbearalg_debug.py -key 4387 

I just uploaded polarbearalg_10.py</br>
To run: python3 polarbearalg_10.py -key 4387  (or -keysize 16 to generate a random 16 bit modulus)</br>

Update: Oops, there is still a bug with the quadratic character base. Let me debug it and fix it. I see now how the polynomials can be used to find a square multiple of a smooth with the linear algebra step (and i've earlier figured out how to derive the factorization of N from those) ... I'm confident now it is do-able. Let me debug though what is happening. I'm sure it's something simple and easy to fix.

Update2: Oh nvm, I was able to generate this example (ignore the root and eq1 in the debug info btw... I need to fix that): 

[Debug Info]Adding to left mod N: 947 adding to right mod N: 1556 y0: 237066 z: 19 discriminant: -124355 qchar: [0, 0, 0, -1, 0, 0, 0, -1, -1, 0] root: 167931 eq1: 11</br>
[Debug Info]Adding to left mod N: 1941 adding to right mod N: 389 y0: 138216 z: 769 discriminant: -6093395 qchar: [0, 0, 0, -1, 0, 0, 0, -1, -1, 0] root: 207125 eq1: 17</br>
sqrt_right:  870485</br>
sqrt_left:  60021</br>
[i]Debug info: disc%N: 1038 (disc): 757744135225 jacobi_symbols: [1, 1, 1, 1, 1, 1, 1, 1, 1, 1] prod_left: 3602520441 prod_right: 4283702500 z: 619 ja: 1 y0: 60021</br>
[i]Found a square in the integers, the quadratic character base worked</br>
[SUCCESS]Factors are: 43 and 47</br>

Here we find a square multiple with the discriminant but distinct coefficients (coefficients that generate a square multiple with the discriminant but are not a multiple of eachother). So it will find these... let me refactor the PoC now... we need to work with these 0 solutions like QS_Variant does.... eventually I'll need to get it to work, such that it can assemble multiple results into a square multiple... but I see the number theory behind it now and how we can work mod m instead of purely in the integers. Because in the example above if we have discriminant: -124355 and discriminant: -6093395, since -6093395 \* 49 = -124355, we would have also found this relations if we had reduced everything mod -124355. 

UPDATE!!!!!!!!!!!!!! UPDATE!!!!!!!!!!!!!!!!!: I changed that last PoC to only use 0 solutions instead of all quadratic residues when I do create_hashmap. So it is now the same as in the QS_Variant. I then also changed a few things around.. and something cool is happening now. I'm doing the quadratic character base calculations the same as in v10, but I'm having one side of the congruence fully in the integers and the other one reduces mod m. And the linear congruence step is producing some really interesting results now... like... it's definitely working now. But I need to debug the details a little before I upload tomorrow. But it 100% is doing something right now, I'll get to the bottom of it tomorrow. PREPARE FOR FACTORIZATION! HAHAHAHHAAHHAHA. 
