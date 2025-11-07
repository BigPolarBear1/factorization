DISCLAIMER: AT NO POINT IN MY 2.5 YEARS OF RESEARCH DID I USE AI OR HAS AI IN ANY WAY HELPED WITH MY RESEARCH. I'M PUTTING THIS HERE BECAUSE I KNOW PEOPLE WILL ACCUSE ME OF IT. 

Note: I am antifa, leader of the cipherpunk, fuck the FBI department. 

Note2: It's basically done now, the NFS variant is nearing completion soon, nothing novel research wise has to be figured out now. I just need to refactor a few things and I'm done.
I have been unemployed for 2 years. I am looking for work. I would like some way where I can have enough time to attack ECC and DLP next. In addition, for a payment, I can also show people how I appraoched factorization from start to finish. Because the method of how I conduct my research and consistently have success is just as valueable as the research itself. And I follow the same patterns over and over again, no matter what subject matter I attack. I don't have a very high IQ and I'm not especially bright or good at doing numbers in my head. But I have insights into research that help me succeed over and over again. I'm willing to teach people how to do it, as long as it enables me to have the finances to continue doing my research aside from that. Contact: big_polar_bear1@proton.me.

#### To run from folder "QS_variant" (Standard SIQS with our number theory as backend):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 200 -base 6000 -debug 1 -lin_size 100_000 -quad_size 1</br></br>
This version is somewhat done, but it doesn't achieve much of an advantage, as the main strength of my work will be the NFS_Variant, which performs NFS with reducible quadratic polynomials<br><br>
#### To run from folder "Improved_QS_Variant" (Sort of a QS/NFS hybrid that improves on Quadratic Sieve using the number theory I developped):</br></br>

I just uploaded polarbearalg_v34.py</br>
To run: python3 polarbearalg_v34.py -keysize 20</br>

I fixed the PoC..... lol. 
It was actually extremely simple. I was way too overcomplicating everything inside my head.
Its done! MY WORK IS COMPLETE! Check v34, it has the linear algebra working.

Update: So using zx^2-N, we can garantuee very small smooth candidates. Since we can add the modulus to either z or x. However, zx^2 must also be square, but there is a number of tricks we can use there to garantuee its factorization.Since x^2 is always square, when performing trail factorization, we only need to factorize z. And since we can add the modulus to x, we can keep that z value small. Its easy. I'll go for a run, take a break, drink some beer to celebrate having finally solved it. And I'll begin working on an optimized PoC this weekend. It will blow any public factorization algorithm completely out of the water... you will see. I'm happy it is finally over, but part of me is really cursing that the solution is this simple and I just didn't see it for this long. I literally cannot believe I didn't see it until now. Because now that I do see it, I'm like "WTF idiot". I guess that's usually the case... same with software bugs... even if they are obvious in retrospect... it still takes a healthy dose of luck of looking at it from the right angle to find it.
