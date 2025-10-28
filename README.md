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

I just uploaded polarbearalg_v19.py</br>
To run: python3 polarbearalg_v19.py -keysize 14</br>

Finally, with v19 we are now doing the RIGHT appraoch. I still need to implement linear algebra and proper sieving. But this demonstrates working with square moduli and 0 solutions. Where we look for a square mod m, and then just multiply the modulus with the square and take the GCD. Simple as that. IT IS SO FUCKING SIMPLE. And I also know how I can finish the algorithm now with linear algebra. It was just a matter of look at this problem from "the right angle".

Expect the full algorithm to drop this week. I have court ordered therapy tomorrow (bc I sent angry emails to the FBI and they are pieces of shit, come fight me losers), and I'll run 30k to and from therapy, so that will take up most of my day.. but I'll have some time to properly think about all the implementation details in my head while running.

Sad day for the NSA, FBI, CIA, DIA, and all the other cunts. The amerian cryptologic advantage, destroyed by a single polar bear. I am stronger then all of america combined, HAHAHAHAHAHAHAHAHHAHAHHAA.

I know this is it. This is the thing I overlooked for half a year. Using square moduli and finding a square mod m. It's so straight forward and simple... but I guess I just missed it for half a year. All my heroes, Galois, Euler, Fermat and Legendre would be laughing at me if they knew it took me half a year... but that's fine, I succeeded in the end. And now we're here. I hope you have a good day today western infosec, because the storm is coming now. And it's time to pay the price.

Ps: In v19... the discriminant im calculating there, that can now also be swapped out with polynomials of the shape zx^2-yx. Without appending +N... hence resulting in a small integer value. Its easy now to finish it. Its very similar to NFS. Its literally NFS with reducible quadratics.

The biggest irony of all this is that I tried so damn hard to sell my work.. yet I was forced down this path. And I think we all know why. There is no way I will stay in the west after all this. No way.

Anyway.. if I have an "accident" this week.. garantuee you it was cowardly americans.

