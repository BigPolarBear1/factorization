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

I just uploaded polarbearalg_v21.py</br>
To run: python3 polarbearalg_v21.py -keysize 14</br>

Right side of the congruence is now 0 solutions for y1 and square moduli.
Right side if perfect already.

Left side, in v19 we used the discriminant there. 
But I am now reducing it to the polynomial:

(y1 will always be 0)
x^2+0*x-(N%mod)\*z

Note that we also modulo reduce N. And for any z, we recalculate the root by multiplying it by z mod m. This way we can move z to (N%mod)\*z instead.
This will allow us to discard many of the bad solutions. Not all. But many. But with a quadratic character base, we can reduce this further.

I will now continue re-implementing the linear algebra step. It's simple now :). 

Another sad day for the NSA, FBI, CIA, DIA, fuck america and micosoft. There was so many chances to avoid this, but none were taking. Which doesn't suprise me with western leaders constantly crying about "dudes in dresses". How does it feel to be destroyed by a dude in a dress huh? Not so funny now? Fucking losers. Come fight me hahaha.

Somewhere in the US, a nazi cryptologist has to review the code I uploaded today... and he's going to have a really bad day today. I guess yesterday was a bad day too.. but with that discriminant further reduced... well.. life's a bitch, sucks to be a nazi.

Update: added some further improvements to v21. It will succeed now if the modulus is large enough. But since everything is properly aligned now the way it should be.. we should be able to assemble a sufficiently large enough modulus now with linear algebra. This should now set me up perfectly to finally get that linear algebra portion working.



