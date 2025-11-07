DISCLAIMER: AT NO POINT IN MY 2.5 YEARS OF RESEARCH DID I USE AI OR HAS AI IN ANY WAY HELPED WITH MY RESEARCH. I'M PUTTING THIS HERE BECAUSE I KNOW PEOPLE WILL ACCUSE ME OF IT. 

Note: I am antifa, leader of the cipherpunk, fuck the FBI department. 

#### To run from folder "QS_variant" (Standard SIQS with our number theory as backend):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 200 -base 6000 -debug 1 -lin_size 100_000 -quad_size 1</br></br>

See below for an improved way of performing what this PoC does.. I'll delete this Proof of Concept once the PoC for Improved_QS_Variant matures a little more<br><br>
#### To run from folder "Improved_QS_Variant" (Sort of a QS/NFS hybrid that improves on Quadratic Sieve using the number theory I developped):</br></br>

I just uploaded polarbearalg_v36.py</br>
To run: python3 polarbearalg_v36.py -keysize 20</br>

I fixed the PoC..... lol. 
It was actually extremely simple. I was way too overcomplicating everything inside my head.
Its done! MY WORK IS COMPLETE! Check v36, it has the linear algebra working.

Update: So using zx^2-N, we can garantuee very small smooth candidates. Since we can add the modulus to either z or x. However, zx^2 must also be square, but there is a number of tricks we can use there to garantuee its factorization.Since x^2 is always square, when performing trail factorization, we only need to factorize z. And since we can add the modulus to x, we can keep that z value small. Its easy. I'll go for a run, take a break, drink some beer to celebrate having finally solved it. And I'll begin working on an optimized PoC this weekend. It will blow any public factorization algorithm completely out of the water... you will see. I'm happy it is finally over, but part of me is really cursing that the solution is this simple and I just didn't see it for this long. I literally cannot believe I didn't see it until now. Because now that I do see it, I'm like "WTF idiot". I guess that's usually the case... same with software bugs... even if they are obvious in retrospect... it still takes a healthy dose of luck of looking at it from the right angle to find it.

Ps: I don't know how my enemies will react now that I've solved it, but I do know it won't be pretty. But I can garantuee you, I won't go down without a fight, so come at me you fucking losers. I'll fight you all to the death. I'm not scared of the NSA, FBI, CIA, DIA.. all the fucking NATO wankers. Come then, come fight me pussies. Morrons. I piss on you all. You fucking losers. You forced me down this path. You are responsible for this. Because you just can't accept the reality of what is happening. I'll fight you all to the death. Fucking pussies. Cowards. And when you do come after me, bring some guns, you'll need them, like the cowards you all are. I'll end you people with my bear hands.

Pps: Just so you assholes know, I will move fast on this now this weekend. And I will make a lot of noise about this once the PoC hits novel performances. It's all far too late dumbasses. Morrons. Idiots. Don't come near because I will punch you morrons in your face. WHAT DID YOU THINK WAS GOING TO HAPPEN YOU STUPID MORRONS?????
