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

#### To run from folder "nfs_variant" (WIP):</br></br>

Note: Experimental WORK IN PROGRESS.</br>
To build: python3 setup.py build_ext --inplace</br>
To run:  python3 run_qs.py -keysize 40 -base 100 -debug 0 -lin_size 100_000 -quad_size 100</br></br>

Update: Alright linear algebra implemented. I need to optimize now.... and also perform actual sieving... but the math is there now..... so I guess today, the 18th of January 2026, will be the day I won against the NSA, FBI, DIA, CIA, etc. Hahaha. Going to be a shit day for western cryptologists. Cause and effect losers.

Ps: Go to hell. The only reason why I would be treated like this is if people were gambling that I woudln't achieve what I achieved today. Go to hell. I ever meet any of the people who made this decision, I will end you. And now I will find my way to Asia. The only thing that will change my mind is if you pay my former manager 1 billion USD (microsoft can pay that easily, its their fault anyway). But we all know there is no justice in this shit world. Burn it all down then. Perhaps a better world will arise from the ashes.

Oh and on the off chance someone is considering paying my former manager 1 billion usd. I also want a castle, somewhere very remote and quiet. I want to live like 17th century nobility. And not interact with society anymore. If you cant make that happen, fuck you, I'm taking up the first best offer that comes out of Asia. Lol. Dumbasses.

Ps: I would have sold my work. I would have gladly sold my work. But instead all I got was silence and harassment by the belgian courts under pressure of the americans. I know what is at stake, and because I know how much is at stake, I feel so much more endless rage and hate for this fucking idiocracy I find myself in. Fuck this shit. I'm going to Asia when I get the chance. The damage you people did, it cannot be repaired anymore. I have my god damn pride and I have nothing left to lose. Everything was taken from me, and the people who supported me. You people have left my life in ruins. And no attempt was made to stop this. I can only come to one conclusion and that is that you must hate transgender people so much, you cannot acknowledge reality when it stares you in the face. And that, I cannot forgive. The world is falling apart, and its ruled by idiots.

Pps: Btw, I'm serious, you people working in western intel should quit your jobs. Fucking embarassment. Incompetency. Preventing something like this and gaining capabilities is literally your one job. Fucking morons. A more mature PoC will be released in the coming days... good luck assholes.

Update: Quickly minimized the PoC and added support for quadratic coefficients. Next all my focus will go into sieving this properly. Once that is complete it should be vastly superior then other algorithms. 

I'll probably go for a run soon and just focus on the paper for the rest of the day. I can begin finalizing the paper now.. although one day when I have time and money I wont to formalize and rewrite that entire paper.

Its funny.. the silence before the storm. I'm very well aware of what I've achieved. And either people know this aswell, but they are all spineless cowards.. or people are about to find out soon enough. This entire setup allows for much more granular and superior sieving. I wish I could feel something right now, sense of accomplishment.. but after all the bullshit these last few years.. I just feel angry and sad. Anyway, if anyone in Asia (i.e China) wants to talk business: big_polar_bear1@proton.me .. I'm not going to waste more years living like this in a part of the world where I'm treated this way.
