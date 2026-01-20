Disclaimer: At no point did AI contribute anything to this research project. Copilot can't even calculate basic congruences without making mistakes. People who think AI can do novel math research are delusional. I did not use AI for the paper, and I did not use AI for any of the code. I attacked this problem, from the very basics 3 years ago, slowly burrowing deeper and deeper. I know people will say lies because they will never let a person like me win. Not when it comes to something as important as this. Because that would literally go against all the bullshit propaganda the far right has been spreading for years. Also, this could have been avoided. Its a disaster of your own making. 

Contact: big_polar_bear1@proton.me , I am looking for employment. Mainly to continue my math research. I am willing to teach my methodologies also. In addition I am willing to renounce my Belgian citizenship, because people in the west must have known my math was correct from the start, yet I was treated like absolute shit for years. I dont ever again want to call myself Belgian. As far as I can see, Belgium is a country of spineless cowards. Imagine doing this to one of your own citizens lol.. Europe could have had this algorithm and math to themselves. But this is pretty much inline with most of my life in Europe, never had oppurtunities here.

I dont know what lies people will tell to take away this achievement from me. But they are just that, lies by people who cannot accept that I did this. Even though I sacrificed 3 years of my life, working fulltime on this. They will tell lies, they wont ever let a transgender person win, let alone one who dropped out of highschool. Its lies by lesser men (and women), who lack ambition and perseverance themselves. I did this. I did this through sacrifice most people cannot imagine. I gave everything for this and I have nothing left. And I piss on everyone who will attempt to deny me what I've achieved, you pathetic weaklings.

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

Update: I'm insanely stressed today. I made a breakthrough, and its sitting right here uploaded on github. True, the code needs some more polish, I'm just really stressed. I worked out the math to find similar smooths.. and its actually really really really easy now. There's a lot of fortunate things about sieving this way that line up nicely. Also notice how the factorization of the constant in the discriminant is partially decided by the factorization of the root of the quadratic polynomial we're sieving with..  everything is just lining up perfectly for superior sieving. I made a breakthrough, I did it! Just stressed man. Anyway.. two days since disclosure, nobody seems to know yet. Either that or people are told to ignore my work and not speak about it. 

Update: God damnit. I really cant finish this code. Because I FUCKING HAVE A BREAKTHROUGH. THIS IS SUPERIOR TO ANY OTHER TYPE OF SIEVING. I HAVE A FUCKING BREAKTHROUGH AND THE WORLD SEEMS BLIND TO IT. And I can proof it, by just finishing the code. But I'm so stressed I feel like fucking throwing up.  I've done it. I've really done it. I trippled checked the math today, that fucking factorization of the root and root plus the linear coefficient is what determines the factorization of that constant in the discriminant.. I know exactly what this means. This is the best way of sieving ever invented. Nothing fucking comes close. I made a breakthrough, and the world seems to be sleeping. I'm so stressed man.

I guess I'm also just nervous that someone else is just going to claim credit for my work. Because anyone with more reach then me, can just copy paste my work and claim it theirs. I should finish the PoC and Paper this week and make some noise about my work. I'm really stressed out of my mind right now. I'm so stressed, I literally feel my fight or flight triggering, even though I'm in bed with my laptop and there is no danger. Yet I'm fucking freezing up and I really just want to close my laptop and not think about math. How does one finish something like this, while being fully aware of its consequences? I mean.. nobody tried to stop me, or ask me to wait.. or offered to buy my work when that was still an option.. its not like another choice was available. But its going to absolutely ruin my life, I know that much.

Update: I've added some more code. You can see how the factors of the root show up in the final smooth. The rest of the factors are determined by the quadratic coefficient, the quadratic polynomial value, and the root + the linear coefficient. So basically, ever factor in our smooth, we know EXACTLY, where this comes from and we CAN sieve for it. GUYS, DO YOU KNOW WHAT THAT MEANS? DO YOU?? Jesus christ, I cant believe nobody tried to stop me or negotiate with me lol.
