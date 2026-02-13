Disclaimer: At no point did AI contribute anything to this research project. Copilot can't even calculate basic congruences without making mistakes. People who think AI can do novel math research are delusional. I did not use AI for the paper, and I did not use AI for any of the code. I attacked this problem, from the very basics 3 years ago, slowly burrowing deeper and deeper. I know people will say lies because they will never let a person like me win. Not when it comes to something as important as this. Because that would literally go against all the bullshit propaganda the far right has been spreading for years. Also, this could have been avoided. Its a disaster of your own making. 

#### To run from folder "NFS_Variant" (Will be adding more NFS related code in the coming days/weeks):</br></br>

Note: Experimental WORK IN PROGRESS.</br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 50 -base 500 -debug 1 -lin_size 100_000 -quad_size 1</br></br>

NFS related code is borrowed from (sieve(), trial()): https://github.com/basilegithub/General-number-field-sieve-Python (note: Very impressively written, helped me big time, thanks)

I've replaced the sieving portion in my PoC with the sieve function of the NFS algorithm. 

To do:

1. Fix quadratic coefficients. (easy)
2. Instead of using zx^2+yx-N = poly_val , change it to use zx^2+yx+poly_val = N? Or atleast make sure we are working with irreducible quadratics. (easy)
3. Optimize sieving parameters (easy)
4. Figure out if the traditional setup where we end up taking a square root over a finite field would yield an advantage. And if so, can we do it while sieving many different quadratics? (less easy)

I'll see how much of the NFS algorithm I can successfully integrate into my work. And if I can find a way to gain an advantage somehow over existing algorithms.

#### To run from folder "CUDA_QS_variant":</br></br>
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

This was an attempt at finding smooths with similar factorization using an SIQS variant. By using quadratic coefficients. But it didnt end up working as I had hoped so I abondoned this approach, but perhaps someone will get some use out of it.

#### To run debug.py" (Prints the linear and quadratic coefficients to solve for 0 in the integers, for use with my paper):</br></br>

To run: python3 debug.py -keysize 12

This basically creates a system of quadratics. Solving them mod p is easy. But there is only one root solution (the factor of N) which solves the system for 0 for any mod p (aka solves it in the integers). Figuring out how to exactly do this quickly is still an ongoing area of research for me. And if a polynomial time algorithm for factorization exists, it is likely done by solving this system of quadratics. Finding a polyomial time algorithm is my ultimate goal, as this would make progress toward solving p = np as well. 

</br></br></br></br></br></br></br></br></br></br></br></br></br></br></br></br></br>
-----------------------------------------------------------------------------------------
Contact: big_polar_bear1@proton.me , I am looking for employment. Mainly to continue my math research. I am willing to teach my methodologies also. In addition I am willing to renounce my Belgian citizenship, because people in the west must have known my math was correct from the start, yet I was treated like absolute shit for years. I dont ever again want to call myself Belgian. As far as I can see, Belgium is a country of spineless cowards. Imagine doing this to one of your own citizens lol.. Europe could have had this algorithm and math to themselves. But this is pretty much inline with most of my life in Europe, never had oppurtunities here.

I dont know what lies people will tell to take away this achievement from me. But they are just that, lies by people who cannot accept that I did this. Even though I sacrificed 3 years of my life, working fulltime on this. When I started this project, I had everything. Now I have nothing, and I'm sure as hell not going to let lesser humans take my achievements away from me either. Losers going to be mad because they dont have the courage themselves to put everything on the line when all the odds are stacked against them.

I was thinking, you know, when I hiked northern scandinavia, those wide open landscapes, its weird, I always felt at home there. Like I had lived a thousand lifes there. Ofcourse one could imagine Ice Age Europe having looked a lot like Lapland does today. I wonder if genetic encoding, just like the fear of snakes and other dangerous animals, also encodes familiar places and places of safety. I wish I could wander that place for the rest of my life. I have no more ambitions in the human world. I want to go wander that place, and live like my ice age ancestors. I love the tundra. Nothing can hurt me there. It is strange to wander a place, I had never been before, yet where everything feels like I've been there before. Perhaps nostalgia for simpler times, hunting woolly mammoths and eating berries and living in the magnificience of unspoiled nature. Why must I have been born in this century, when there was thousands of years, of simpler lifes I could have been born into.  *sighs* Oh well, 3 am, I should sleep.
