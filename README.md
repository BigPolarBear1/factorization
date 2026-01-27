Disclaimer: At no point did AI contribute anything to this research project. Copilot can't even calculate basic congruences without making mistakes. People who think AI can do novel math research are delusional. I did not use AI for the paper, and I did not use AI for any of the code. I attacked this problem, from the very basics 3 years ago, slowly burrowing deeper and deeper. I know people will say lies because they will never let a person like me win. Not when it comes to something as important as this. Because that would literally go against all the bullshit propaganda the far right has been spreading for years. Also, this could have been avoided. Its a disaster of your own making. 

#### To run from folder "Improved_Sieving":</br></br>

Note: Experimental WORK IN PROGRESS.</br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 70 -base 500 -debug 1 -lin_size 10_000_000 -quad_size 1</br></br>

Update: Actually the linear shifting I'm doing in the uploaded PoC is a really good idea. I also have an idea. We should be able to sieve that linear shifting using that hashmap thingy. Let me begin working toward that now.

Update: Oh yea.. I should calculate my hashmap like I used to do in the very beginning. Not just where the discriminant yields a 0 solution mod p.. but also quadratic residues.. because for each of those will exist a quadratic that is 0 mod p. Then go from there. 

Update: I'll upload a newer version in the coming days. So I got some code now (not uploaded yet) which will precalculate the roots and linear coefficient which produce a solution of 0 mod p with the quadratic (not just the discriminant like the uploaded PoC does in solve_roots()). Then I can index that by root residue... and it will give me information about the residues of x+y and the factors of the polynomial values. Then that construction should allow me to finally finish the algorithm properly. Anyway, got to run 30k tomorrow all day to go to therapy and back, because the belgian justice system thinks I'm insane or some shit. I feel like people have been trying really hard into gaslighting me that something is wrong with me and that I should stop doing math. And sometimes, this also create a lot of self-doubt... you know "what if i've just lost it" .. and that then brings forth a lot of dark thoughts. But I know I got it now.. and I know I was right about my math and that everything thats going on, is just people being shitheads.

UPDATE: OH OFCOURSE. That x+y thats just the other factor.

Ie:

41^2+66\*41-4387 = 0      x+y = 107</br>
107^2-66\*107-4387 = 0    x+y = 41 </br>

So x+y just produce that other factor. Hence x+y is the root for the quadratics with the sign for the coefficient flipped. 
And there we go.... there we go. Let me add that to the paper tonight because that does give us some much better clues how to finish the algorithm.

SO basically we need to find two roots (for zx^2+yx-Nk and zx^2-yx-Nk) for the same linear coefficient, both of which must factorize over the factor base and produce a small polynomial value.. and we can change the linear coefficient to basically add a linear offset to the polynomial value. Its actually extremely elegant. I should have all the info now to produce an awesome PoC before the end of the week. I will make it my mission to do this. I am so tired.. when I worked at microsoft, everyone seemed happy for me, people kept reaching out all the time. Nobody talks to me anymore. And a lot of people seem to think I've lost my mind. I cant live like this anymore. To feel like a dissappointment. This ends this week. This will be the week this nightmare finally ends. I have everything now, I can do it. 

Update: Today was therapy and running 30k and then crashing out in the evening. I did do some further thinking, and I have a pretty definitive idea now on how to setup sieving. So tomorrow, its time to end this, I'm going to start making real progress on this. If anyone has anything to say before I go down that path... better do it in the coming hours because the moment I wake up tomorrow... shit is going to escalate soon, I can promise you that. I know what I got and I know I'm right, and I am tired of this humiliation. This ends tomorrow.

Btw I have figured out how to sieve the polynomial value and x+y at the same time. Knowing that x+y is simply another root for a quadratic with the sign flipped on the linear coefficient, it was fairly simple to piece together.. things are going to escalate very quickly now.. and I will make sure to contact my non-nato friends. Because I feel nothing but resentment for how you people treated me. I feel deep hatred. No idea what life has been like. The struggles of those past years. The darkness. You people are going to pay the price, one way or another, shit has consequences. If you threaten me with a gun, fire me, seperate me from all my friends, make sure I cant find 0day buyers or employment, fire my former manager bc he once gave me a promotion... a price will be paid. You think you can treat people like animals? Like disposable trash? Morrons. Morrons, just like all your big tech executives... fucking morrons all of them. The world is ruled by shitheads with half a brain. My patience has run out. This is the beginning. I'm never stopping. I'm not afraid of the west. I'm not afraid of little men like pete hegseth. Fuckign weaklings. Come fight me then you fucking cunts.

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
