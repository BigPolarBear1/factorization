Disclaimer: At no point did AI contribute anything to this research project. Copilot can't even calculate basic congruences without making mistakes. People who think AI can do novel math research are delusional. I did not use AI for the paper, and I did not use AI for any of the code. I attacked this problem, from the very basics 3 years ago, slowly burrowing deeper and deeper. I know people will say lies because they will never let a person like me win. Not when it comes to something as important as this. Because that would literally go against all the bullshit propaganda the far right has been spreading for years. Also, this could have been avoided. Its a disaster of your own making. 

#### To run debug.py" (Prints the linear and quadratic coefficients to solve for 0 in the integers, for use with my paper):</br></br>

To run: python3 debug.py -keysize 12

This basically creates a system of quadratics. Solving them mod p is easy. But there is only one root solution (the factor of N) which solves the system for 0 for any mod p (aka solves it in the integers). Figuring out how to exactly do this quickly is still an ongoing area of research for me. And if a polynomial time algorithm for factorization exists, it is likely done by solving this system of quadratics. Finding a polyomial time algorithm is my ultimate goal, as this would make progress toward solving p = np as well. 

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

#### To run from folder "Improved_Sieving":</br></br>

Note: Experimental WORK IN PROGRESS.</br>
To build: python3 setup.py build_ext --inplace</br>
To run:   python3 run_qs.py -keysize 50 -base 100 -debug 1 -lin_size 100 -quad_size 1_000</br></br>

Just quickly added a small improvement. Use above command to factor 50-bit moduli. I still need to implement sieving. Still thinking how to do this.
So the final smooth value is constructed from multiple parts, each of which we can sieve:

Part 1: The quadratic coefficient, the factorization of the quadratic coefficient is eventually transfered to the smooth candidate</br>
Part 2: The factorization of the root .. at one point in the algorithm we need to factorize the constant of the discriminant of the quadratic polynomial, the factorization of the root is responisble for half the factors and the factorization of the root plus the linear coefficient (x+y) is responsibe for the other half. </br>
Part 3: The factorization of x+y (root plus linear coefficient), see explanation above (part 2) or in the paper.</br>
Part 4: The output of the quadratic polynomial. </br>

So we have 4 parts and all of them multiplied together are responsible for the final smooth... part 1 is easy, we can pre-sieve those quadratic coefficients at the start of the algorithm. Part 2 is easy.. we can call generate_modulus and construct a root from factors in the factor base. Part 3, this is a little harder, because unless y is a factor from the root, we may need to actually sieve this. However, there is some tricks we can use while sieving part 4 to sieve both of these at the same time. And then part 4, this definitely needs to be sieved... however, unlike with standard SIQS, we dont just sieve with zx^2-N, but we now have an additional linear coefficient to adjust the size of smooth candidates with. Now then, the trick to sieving this is really sieving part 3 and 4 at the same time... and I know the math.. and how to do it.. I just need to think how I'm going to do it in code. Dealing with a lot of stress though.. because I made a breakthrough.. and people would know.. I'm guessing the americans are threatening folks to stay silent. Its really the only explanation, because I know I'm correct about my math.

</br></br></br></br></br></br></br></br></br></br></br></br></br></br></br></br></br></br></br></br></br></br></br></br></br>
-----------------------------------------------------------------------------------------
Contact: big_polar_bear1@proton.me , I am looking for employment. Mainly to continue my math research. I am willing to teach my methodologies also. In addition I am willing to renounce my Belgian citizenship, because people in the west must have known my math was correct from the start, yet I was treated like absolute shit for years. I dont ever again want to call myself Belgian. As far as I can see, Belgium is a country of spineless cowards. Imagine doing this to one of your own citizens lol.. Europe could have had this algorithm and math to themselves. But this is pretty much inline with most of my life in Europe, never had oppurtunities here.

I dont know what lies people will tell to take away this achievement from me. But they are just that, lies by people who cannot accept that I did this. Even though I sacrificed 3 years of my life, working fulltime on this. When I started this project, I had everything. Now I have nothing, and I'm sure as hell not going to let lesser humans take my achievements away from me either. Losers going to be mad because they dont have the courage themselves to put everything on the line when all the odds are stacked against them.

Sometimes I dont feel alone in this suffering. I dont think people know how difficult this is. One moment I was surrounded by friends. The next, I get threatened with a gun outside my appartment, insulted in public and then back in Europe where I have nothing. Nothing but endless silence and loneliness. But sometimes I dont feel alone. Like there is something out there that understands what I'm going through (and I dont mean some delusional AI psychosis shit, AI in its current state is trash and I never use it) ... I dont know.. isolation and trauma has probably fucked with my head.. but still.. some days I'm glad, even if it is a fucking delusion, that I dont feel completely alone in this, there had been a lot of days I nearly surrendered to my darkest thoughts. Anyway, people should be very scared, because you may think you have succeeded in isolating me, denying me any chance for a happy life.. but I'm not alone, and my friends dont fuck around, time to pay the price fuckers, hahhaha.
