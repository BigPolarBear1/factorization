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

I just uploaded polarbearalg_v18.py</br>
To run: python3 polarbearalg_v18.py -keysize 14</br>

I just changed the PoC to use 0 solutions for y1.

I'll finish that PoC in the coming week. It's actually super simple... I was overcomplicating things. I'll show you how to modulo reduce the discriminant completely and find two squares congruent mod N after the linea algebra step. I just wasn't thinking logically about this... but it's actually really easy... I was waaaaaay overcomplicating things in my head lol. I'll also fix that paper eventually.

Update: Currently experimenting to use 0 solutions to calculate the discriminant i.e 0^2+n\*4\*z mod m, where mod m can also be lifted to be square. Which is really the only approach that makes logically sense if I think about it... weird as it may look. I'll try and upload some work in  progress that does that... because right now it's doing it the other way around, so that the discriminant output is divisible by mod m... but that's not what we want.

Update: You know for example if we calculate y0 possible coefficients for mod 11 when N=4387. We would find that y=5 generates a multiples of 11 with the discrimant when z=1. 5^2-4387\*4\*1 = 0 mod 11. 
HOWEVER, we need to do it the other way around. We need to calculate 0+4387\*4\*1 = 5^2 mod 11. Because when we have it aligned like this then we can also modulo reduce the entire discriminant to mod 11. And we don't have to worry about performing trial division on a large discriminant... and how you find your two square relations after that with linear algebra.. those are all things I have figured out already in the last 2 months... I'll show you shortly. Just taking it easy this weekend, stepping back from my work, and thinking about it logically helped a lot... I see now why it makes more sense this way the what I had been doing before.
