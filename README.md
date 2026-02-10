Disclaimer: At no point did AI contribute anything to this research project. Copilot can't even calculate basic congruences without making mistakes. People who think AI can do novel math research are delusional. I did not use AI for the paper, and I did not use AI for any of the code. I attacked this problem, from the very basics 3 years ago, slowly burrowing deeper and deeper. I know people will say lies because they will never let a person like me win. Not when it comes to something as important as this. Because that would literally go against all the bullshit propaganda the far right has been spreading for years. Also, this could have been avoided. Its a disaster of your own making. 

#### To run from folder "Improved_Sieving":</br></br>

Note: Experimental WORK IN PROGRESS.</br>
To build: python3 setup.py build_ext --inplace</br>
To run:  python3 run_qs.py -keysize 60 -base 500 -debug 1 -lin_size 10_000 -quad_size 10</br></br>

Update: Removed the jacobi symbols again. Didn't like it and I dont think its needed anyway. I think the smarter way is sieving zx+y by also using the quadratic coefficient z to change the residue of zx. Let me work out some math. The only other thing I can think of is to reduce the constant N in zx^2+xy-N to mod m ... just reduce everything mod m and then try to find a valid solution mod m. I'll mess around a bit more. Something is going to work eventually.

Update: Thinking about it some more.. I could just reduce the entire quadratic mod m.. and THEN calculating jacobi symbols for primes not part of the modulus actually makes sense. Why am I only thinking of this just now? Damnit. This could actually work. And that also solves all of the issues I'm currently having. And even if we have just a solution mod m, but because we can take jacobi symbols for other primes we can then probably use p-adic lifting to calculate the solution in the integers.  Damnit. I'm too slow. Its funny I cant find employment while people of mediocre skills and no willpower or perseverance are paid 6 figure salaries, but then again, people are straight up cunts. These mediocre losers at microsoft blamed my manager for giving me a promotion and then fired him, even though I was finding bugs in secure channel and openssl etc. And then they retroactively went back and actually downgraded some of my bugs.. I guess to justify the shit they just pulled. I feel such immense hate and bitterness. Those people should stay very far away from me, for their own safety.

Those events in 2023 and the hopelessness, social isolation and unemployement in the years after have permanently changed me as a person. And if I had the money, I would go live by myself in nature somewhere (preferably the arctic) and never interact with humans ever again. I would teach myself all the survival skills I need, and just quit the human world and be that mean fucker you dont want to come across. I have seen enough misery and ugliness in this shit world and I cant take it anymore.

I'll try tomorrow to sieve everything mod m.. thats something i've tried earlier in my research, but I think with the way I'm doing it now, I may be able to actually pull it off. I'm struggling hard mentally. I wish I could just take a break from work and just run ultra marathons every day. Like just be in the woods all day running. I especially like the woods when its cold outside. Reminds me of my hikes in the arctic and scotland. Simpler times. But I cant. I have no future, I'm not getting any younger, I need to keep grinding. Sometimes it feels like all this shit is crushing my sanity, breaking me and the only thing I know what to do is to keep going, because the alternative is to go into the woods and hang myself. This world is so fucking shit man. 

Update: Oh yea, I just remembered something, this thing about using the polynomial value as modulus that I had found a month ago:

So if N = 4387:

We create the following discriminants for polynomial values, 0, 65,128,189,248: </br>

148^2-4\*(4387+248) = 58^2 and 66^2 = 58^2 mod 248</br>
148^2-4\*(4387+189) = 60^2 and 66^2 = 60^2 mod 189</br>
148^2-4\*(4387+128) = 62^2 and 66^2 = 62^2 mod 128</br>
148^2-4\*(4387+65) = 64^2 and 66^2 = 64^2 mod 65</br>
148^2-4\*(4387) = 66^2</br>

This is actually quite interesting. There is definitely a way to work mod m ... and do something similar to NFS. Fuck it, lets go hard on this tomorrow. So close to figuring it out now. I'm not going to quit right before a big breakthrough. lol. I know that is what people want.. but they can just fucking pay me for my work if they want me to stop putting shit on github.

Update: Oh, something just lit up in my head. Let me get some progress on this tomorrow now. I'm pretty sure I know how it can be done. Using the polynomial value as modulus, that is indeed the correct intuition.. because modulo the polynomial value, the discriminant is aligned with N. Yea yea yea... easy, I got this. 

#### (Outdated, check Improved_Sieving instead) To run from folder "CUDA_QS_variant":</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run:  python3 run_qs.py -keysize 240 -base 100_000 -debug 1 -lin_size 100_000_000 -quad_size 100</br></br>
 
Prerequisites: </br>
-Python (tested on 3.13)</br>
-Numpy (tested on 1.26.2)</br>
-Sympy</br>
-cupy-cuda13x</br>
-cython</br>
-setuptools</br>
-h5py</br>
(please open an issues here if something doesn't work)</br></br>

Additionally cuda support must be enabled. I did this on wsl2 (easy to setup), since it gets a lot harder to access the GPU on a virtual machine.

Update: While making this PoC, I suddenly had an idea, linking back to some things I tried earlier last spring. See below.

#### To run debug.py" (Prints the linear and quadratic coefficients to solve for 0 in the integers, for use with my paper):</br></br>

To run: python3 debug.py -keysize 12

This basically creates a system of quadratics. Solving them mod p is easy. But there is only one root solution (the factor of N) which solves the system for 0 for any mod p (aka solves it in the integers). Figuring out how to exactly do this quickly is still an ongoing area of research for me. And if a polynomial time algorithm for factorization exists, it is likely done by solving this system of quadratics. Finding a polyomial time algorithm is my ultimate goal, as this would make progress toward solving p = np as well. 

</br></br></br></br></br></br></br></br></br></br></br></br></br></br></br></br></br></br></br></br></br></br></br></br></br>
-----------------------------------------------------------------------------------------
Contact: big_polar_bear1@proton.me , I am looking for employment. Mainly to continue my math research. I am willing to teach my methodologies also. In addition I am willing to renounce my Belgian citizenship, because people in the west must have known my math was correct from the start, yet I was treated like absolute shit for years. I dont ever again want to call myself Belgian. As far as I can see, Belgium is a country of spineless cowards. Imagine doing this to one of your own citizens lol.. Europe could have had this algorithm and math to themselves. But this is pretty much inline with most of my life in Europe, never had oppurtunities here.

I dont know what lies people will tell to take away this achievement from me. But they are just that, lies by people who cannot accept that I did this. Even though I sacrificed 3 years of my life, working fulltime on this. When I started this project, I had everything. Now I have nothing, and I'm sure as hell not going to let lesser humans take my achievements away from me either. Losers going to be mad because they dont have the courage themselves to put everything on the line when all the odds are stacked against them.

I was thinking, you know, when I hiked northern scandinavia, those wide open landscapes, its weird, I always felt at home there. Like I had lived a thousand lifes there. Ofcourse one could imagine Ice Age Europe having looked a lot like Lapland does today. I wonder if genetic encoding, just like the fear of snakes and other dangerous animals, also encodes familiar places and places of safety. I wish I could wander that place for the rest of my life. I have no more ambitions in the human world. I want to go wander that place, and live like my ice age ancestors. I love the tundra. Nothing can hurt me there. It is strange to wander a place, I had never been before, yet where everything feels like I've been there before. Perhaps nostalgia for simpler times, hunting woolly mammoths and eating berries and living in the magnificience of unspoiled nature. Why must I have been born in this century, when there was thousands of years, of simpler lifes I could have been born into.  *sighs* Oh well, 3 am, I should sleep.
