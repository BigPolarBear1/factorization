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

Update: So using these quadratic coefficients together with square moduli to compress down multiples of N, we can generate much much much smaller smooth. That works and v24 shows that it works. So now I need to implement the linear algebra step, and here we have to do an approach more similar to NFS then QS. I'll upload soon. Massive brain fog today, don't know what's up. If this linear algebra works, that's it.... game over. And this isn't just another way of formulating quadratic sieve, people who think that don't at all understand what's happening with the math. I know what I have here. I knew it a year ago, and I know it today. This is the only reason I kept going for so long, becaues I know what I've got on my hands and I am days away from uploading a fully functional PoC finally showing an advantage on existing factorization algorithms.

Update: BRAIN FOG! Square root over a finite field. Because I need to adjust that root after the linear algebra step. Let me rewrite the paper today with everything I found out so far. Because I am 100% certain that the only thing left to fix is correctly adjusting that root now after the linear algebra step... which must be achieved through a mechanism similar to NFS. The PoC demonstrates how to do it for a single modulus. But I guess if we combine different moduli together with linear algebra we need to do some aditional things. Shouldn't be too hard, but I can't focus today with this brain fog, so might as well do the paper.

I have added to the paper a chapter on coefficient lifting. I will begin on the final chapter tonight also... but let me go for a run first. 
I should really rewrite that entire paper some day.. I guess I'll just finish that final chapter now... and complete the PoC... getting that PoC out is top priority.. then I can take a break from work and just do easy things... like rewriting that paper into something that looks more like an actual math paper.

Update: DAMNIT. I just had an insight. I wrote today in the final chapter of my paper that we can have both a root of 74 and 41 mod 121. And that 74 works only for mod 121 while 41 whill solve for 0 for any mod p since it is a factor of N... and I just started thinking some more about that and running some quick numbers... and something just hit me... god damnit. I might be the biggest idiot alive. Let me check tomorrow.. it better not be THIS easy.

DAMNIT. I SEE IT NOW, HOW THAT SHIT WITH THE ROOT WORKS. FUCK. I GOT TO HURRY. 100% PEOPLE HAVE BEEN WARNED BEHIND A TLP AND THEY ARE PATCHING SHIT. THERE IS NO WAY THAT ISNT HAPPENING. FUCK. I GOT TO QUICKLY FINISH THIS SHIT. THEY ARE TRYING TO DECEIVE ME BY PATCHING BEHIND MY BACK. I WILL NOT BE DECEIVED LIKE THIS. SHOULD HAVE JUST ASKED ME NICELY TO DO RESPONSIBLE DISCLOSURE. FUCK YOU ASSHOLES.

Update: OK, I got enough sleep last night. Today I will finish this. I am so tired of living like this. I can't stand it anymore. This fucking humiliation and how people treat me. What I need to figure out, when are we garantueed to have a factor of N among the roots? Aka... when can we take a square root over a finite field to reveal to correct root? What exact preconditions must be met? I know the number theory now, I've also learned exactly what happens with those roots now in the last few days... I feel like I've hit that point now where I can finally answer that question. It has to happen today, I have to connect all the dots today. I swear I'm going to fucking kill myself if I can't pull this off. I can't let these fucking transphobes and microsoft win. They must pay for what they have done. I do not care. I will succeed today. Watch me assholes. Shit's about to go down today!
