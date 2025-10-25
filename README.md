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

I just uploaded polarbearalg_v12.py</br>
To run: python3 polarbearalg_v12.py -keysize 14</br>

Update: Ignore that PoC for now. I did some thinking, and I'm going off on a tangent here. I need to regroup. Tomorrow I'll take some example I'm familiar with, i.e 148^2 = 66^2 mod 4387, and modulo reduce both squares, try to figure out how I would assemble those mod m using linear algebra. There's no point in just trying random stuff. Literally, the only times I made progress with my math is when I stop doing chaotic shit, and just think about the problem for a moment. I know I almost have it... I just need to implement it correctly. 

Update: I was thinking some more.... and the way to do this, if I want to leverage what I learned earlier, doing the whole thing with square multiples of a smooth, aka, finding a square modulo a smooth... the way you would do it with this setup is, you would have your integer side on the left side of your congruence, and then on the right side you find a square modulo your integer side. That is probably how that is supposed to work. I reread the chapter on the quadratic character base in matthew brigg's paper on NFS. And they also seem to be finding a square like this.. I was misinterpreting what NFS did before. This type of setup would thus be most like NFS. Let me first generate some examples where the right side is a square mod the left side, and then adjust the coefficients and take the gcd....I got the be deliberate in my math, not just try random shit, and I need to verify my assumptioms first.

Update: You know, I'm a fucking idiot. In one or two months, I could have studied abstract algebra, and worked through some undergraduate course on youtube and read some books about it... and then I could have read matthew brigg's paper on NFS after that... and I would have probably instantly made all the connections to my own work. Instead I feel like I'm tying to decipher some alien scripture. When this is done... I want to atleast be able to understand and read graduate level books on number theory and abstract algebraa and work through those famous MIT lectures on linear algebra. I'm not going to do any more research until I do this. I am fucking shooting myself in the foot here by just yolo-ing it. I mean... I'm doing it, making progress, no matter how slow... but I think if I keep working like this, I am going to have a mental breakdown because I'm making things way too hard for myself.
