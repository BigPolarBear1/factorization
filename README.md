Disclaimer: At no point did AI contribute anything to this research project. Copilot can't even calculate basic congruences without making mistakes. People who think AI can do novel math research are delusional. I did not use AI for the paper, and I did not use AI for any of the code. I attacked this problem, from the very basics 3 years ago, slowly burrowing deeper and deeper. I know people will say lies because they will never let a person like me win. Not when it comes to something as important as this. Because that would literally go against all the bullshit propaganda the far right has been spreading for years. Also, this could have been avoided. Its a disaster of your own making. 

#### To run debug.py" (Prints the linear and quadratic coefficients to solve for 0 in the integers, for use with my paper):</br></br>

To run: python3 debug.py -keysize 12

This basically creates a system of quadratics. Solving them mod p is easy. But there is only one root solution (the factor of N) which solves the system for 0 for any mod p (aka solves it in the integers). Figuring out how to exactly do this quickly is still an ongoing area of research for me. And if a polynomial time algorithm for factorization exists, it is likely done by solving this system of quadratics. Finding a polyomial time algorithm is my ultimate goal, as this would make progress toward solving p = np as well. 

#### (Outdated, check nfs_variant instead) To run from folder "CUDA_QS_variant":</br></br>
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

Update: While making this PoC, I suddenly had an idea, linking back to some things I tried earlier last spring. See below (I think this approach is interesting, but its hard to get an advantage.. however using these different multiples of N.. there is a much easier way to exploit this, hence see NFS_variant):

#### To run from folder "nfs_variant" (note: I'll rename this folder eventually, its really just improved way of sieving using quadratic polynomials):</br></br>

Note: Experimental WORK IN PROGRESS.</br>
To build: python3 setup.py build_ext --inplace</br>
To run:  python3 run_qs.py -keysize 30 -base 50 -debug 1 -lin_size 100 -quad_size 100</br></br>

I will implement proper sieving in the coming days and prepare the paper for release. This is my breakthrough. It finally happened. This setup will allow for much better sieving. You will see in the coming days...

Update: Barely got work done today. Just stressed out of my mind. You would think, when you make a historical breakthrough, that people would notice. This entire situation feels off. I know the west has always hated me, my entire career has been a shitshow, so I probably shoudln't be suprised I'm treated like this the moment I started doing cryptologic research (and succeeding). They probably know what I did, but they would rather watch the world burn then admit it. I figured out how I'm going to appraoch sieving tomorrow.. a way to sieve every part of the smooth candidate. I'm so fucking stressed. I probably should work day and night to finish this now, because the sooner people become aware of what I did... the safer it will be for me and everyone around me. People can behave like savages, just look at the americans, had I lived in any other country then a western country, I would have been murdered already. I have no doubt about that. Not that I'm afraid of monkey brained losers like pete hegseth. I'll fight those cunts.

-----------------------------------------------------------------------------------------
Contact: big_polar_bear1@proton.me , I am looking for employment. Mainly to continue my math research. I am willing to teach my methodologies also. In addition I am willing to renounce my Belgian citizenship, because people in the west must have known my math was correct from the start, yet I was treated like absolute shit for years. I dont ever again want to call myself Belgian. As far as I can see, Belgium is a country of spineless cowards. Imagine doing this to one of your own citizens lol.. Europe could have had this algorithm and math to themselves. But this is pretty much inline with most of my life in Europe, never had oppurtunities here.

I dont know what lies people will tell to take away this achievement from me. But they are just that, lies by people who cannot accept that I did this. Even though I sacrificed 3 years of my life, working fulltime on this. When I started this project, I had everything. Now I have nothing, and I'm sure as hell not going to let lesser humans take my achievements away from me either. Losers going to be mad because they dont have the courage themselves to put everything on the line when all the odds are stacked against them.

Sometimes I dont feel alone in this suffering. I dont think people know how difficult this is. One moment I was surrounded by friends. The next, I get threatened with a gun outside my appartment, insulted in public and then back in Europe where I have nothing. Nothing but endless silence and loneliness. But sometimes I dont feel alone. Like there is something out there that understands what I'm going through (and I dont mean some delusional AI psychosis shit, AI in its current state is trash and I never use it) ... I dont know.. isolation and trauma has probably fucked with my head.. but still.. some days I'm glad, even if it is a fucking delusion, that I dont feel completely alone in this, there had been a lot of days I nearly surrendered to my darkest thoughts. Anyway, people should be very scared, because you may think you have succeeded in isolating me, denying me any chance for a happy life.. but I'm not alone, and my friends dont fuck around, time to pay the price fuckers, hahhaha.
