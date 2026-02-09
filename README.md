Disclaimer: At no point did AI contribute anything to this research project. Copilot can't even calculate basic congruences without making mistakes. People who think AI can do novel math research are delusional. I did not use AI for the paper, and I did not use AI for any of the code. I attacked this problem, from the very basics 3 years ago, slowly burrowing deeper and deeper. I know people will say lies because they will never let a person like me win. Not when it comes to something as important as this. Because that would literally go against all the bullshit propaganda the far right has been spreading for years. Also, this could have been avoided. Its a disaster of your own making. 

#### To run from folder "Improved_Sieving":</br></br>

Note: Experimental WORK IN PROGRESS.</br>
To build: python3 setup.py build_ext --inplace</br>
To run:  python3 run_qs.py -keysize 60 -base 300 -debug 1 -lin_size 10_000 -quad_size 10</br></br>

Update: I added support again for the quadratic coefficient. For the k variable (multiples of N) I'm going to restrict myself to 1. Using that quadratic coefficient is more important. Current PoC is very slow, just quickly hacked something together. So the idea now is that using the quadratic coefficient, we can have much more granular control to shrink the size of the polynomial value.. so let me work on that next. Its the combination of both the linear coefficient and quadratic coefficient that makes this work. That was what I overlooked until now. 

Update: I want to experiment with something really quick. What if I do this:

1. x must factorize
2. the polynomial value must factorize
3. zx+y we will just represent as jacobi symbols

Then I have the factorization of x as one part in my linear algebra step, the factorization of the polynomial value as one part and then the jacobi symbols for zx+y. 
I know I've been trying to find that link to NFS for a while... let me try this again. I really just need to find a way to substitute one part of this entire setup with jacobi symbols... fuck it, I'll make this my goal for next week. 

Oh yea.. you know since zx+y basically represents another root for a quadratic with the sign flipped on the linear coefficient..  I should be able to get jacobi symbols out of that. Ah fuck it, its definitely possible. I substitute zx+y with jacobi symbols (or some other part of this setup) and within a day I'll have an algorithm setting records. 

Anyway, I'll begin tomorrow. 

Update: struggling big time with depression. Its hard to keep going. Its been 3 years now. Just grinding out factorization. Having no income. No job oppurtunities, nobody who reaches out anymore. Social isolation. And I just feel hate and rage. For everything that has happened. If someone could arrange a job in China, I would go. As long as I kow it is legit, and not potentially one of those infosec stalkers harassing and trolling me (fucking losers btw). Anyway.. I think I can do it now, use jacobi symbols... the moment I succeed, this nightmare will be over and I can get a good job in Asia, because I sure as hell want nothing to do with these fucking western losers. Fucking elitist losers. If I ever meet western infosec people, I will end them. And If I ever meet microsoft people, I will also end them. Nobody did or said anything when they came after my manager with lies, the only person in this shit industry who ever supported me and believed in me. Fuck them all. I hope everyone in western tech just drops dead.

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
