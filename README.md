Disclaimer: At no point did AI contribute anything to this research project. Copilot can't even calculate basic congruences without making mistakes. People who think AI can do novel math research are delusional. I did not use AI for the paper, and I did not use AI for any of the code. I attacked this problem, from the very basics 3 years ago, slowly burrowing deeper and deeper. I know people will say lies because they will never let a person like me win. Not when it comes to something as important as this. Because that would literally go against all the bullshit propaganda the far right has been spreading for years. Also, this could have been avoided. Its a disaster of your own making. 

#### To run from folder "Improved_Sieving":</br></br>

Note: Experimental WORK IN PROGRESS.</br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 90 -base 1_000 -debug 1 -lin_size 1_000_000 -quad_size 1</br></br>

Update: I spent the day brainstorming sieving this.
So we now have many moving parts we can use to sieve.
Sieving the polynomial value is trivial. We can increase the linear coefficient, the root and the quadratic coefficient.
But no matter what we do zx+y must also factorize over the factor base. (and z and the root must have known factors, but that is less of a concern as that is trivial to garantuee).

I guess if we start with some zx+y where we know factorization of z, x and zx+y. Next we set zx+y as modulus. Then we sieve the linear coefficient with this modulus (i.e zx+(y+modulus*1)).... so those values just ends up being a multiple of zx+y. Let me do some thinking. I dont think its important to have "known factors" for that polynomial values, because we have many different variables now that we can use to get that below a certain size. So I guess I really just need to focus my efforts on making sure zx+y factorizes.

Update: Eureka! I have an idea. It is so easy to just add a linear offset now to smooth candidates and change their size. Oh this will be quite beautiful. I'll try to upload a PoC tomorrow.

Update: I got it! This is actually very easy now.  I'll upload a PoC tomorrow. I already wrote the math in code just now. It works beautifully. I was right. People would have known I was right. Especially people with more math experience. Yet, this is how I'm treated. Just this fucking agony every day. 

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

Sometimes I dont feel alone in this suffering. I dont think people know how difficult this is. One moment I was surrounded by friends. The next, I get threatened with a gun outside my appartment, insulted in public and then back in Europe where I have nothing. Nothing but endless silence and loneliness. But sometimes I dont feel alone. Like there is something out there that understands what I'm going through (and I dont mean some delusional AI psychosis shit, AI in its current state is trash and I never use it) ... I dont know.. isolation and trauma has probably fucked with my head.. but still.. some days I'm glad, even if it is a fucking delusion, that I dont feel completely alone in this, there had been a lot of days I nearly surrendered to my darkest thoughts. Anyway, people should be very scared, because you may think you have succeeded in isolating me, denying me any chance for a happy life.. but I'm not alone, and my friends dont fuck around, time to pay the price fuckers, hahhaha.

I was thinking, you know, when I hiked northern scandinavia, those wide open landscapes, its weird, I always felt at home there. Like I had lived a thousand lifes there. Ofcourse one could imagine Ice Age Europe having looked a lot like Lapland does today. I wonder if genetic encoding, just like the fear of snakes and other dangerous animals, also encodes familiar places and places of safety. I wish I could wander that place for the rest of my life. I have no more ambitions in the human world. I want to go wander that place, and live like my ice age ancestors. I love the tundra. Nothing can hurt me there. It is strange to wander a place, I had never been before, yet where everything feels like I've been there before. Perhaps nostalgia for simpler times, hunting woolly mammoths and eating berries and living in the magnificience of unspoiled nature. Why must I have been born in this century, when there was thousands of years, of simpler lifes I could have been born into.  *sighs* Oh well, 3 am, I should sleep.
