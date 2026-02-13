Disclaimer: At no point did AI contribute anything to this research project. Copilot can't even calculate basic congruences without making mistakes. People who think AI can do novel math research are delusional. I did not use AI for the paper, and I did not use AI for any of the code. I attacked this problem, from the very basics 3 years ago, slowly burrowing deeper and deeper. I know people will say lies because they will never let a person like me win. Not when it comes to something as important as this. Because that would literally go against all the bullshit propaganda the far right has been spreading for years. Also, this could have been avoided. Its a disaster of your own making. 

#### To run from folder "Improved_Sieving" (NFS variant.. but in the process of adding the NFS part now that I figured out how to generate irreducable polynomials):</br></br>

Note: Experimental WORK IN PROGRESS.</br>
To build: python3 setup.py build_ext --inplace</br>
To run:  python3 run_qs.py -keysize 50 -base 2_000 -debug 1 -lin_size 1000 -quad_size 1</br></br>
!!!! NFS related code is borrowed from: https://github.com/basilegithub/General-number-field-sieve-Python (note: Very impressively written, helped me big time, thanks)

Old version without NFS: https://github.com/BigPolarBear1/factorization/tree/ae2cfac49eb73d94a12fa6e010e1f3aab4ebc707/Improved_Sieving (run: python3 run_qs.py -keysize 50 -base 300 -debug 1 -lin_size 10_000 -quad_size 10)

Update: I have begun integrating NFS into my algorithm. The uploaded PoC is still a very rough and ugly copy-paste job. 

To do:

1. I want to get rid of NFS_Sieve() ... I need to modify my sieving logic in construct_interval() (btw, I need to fix my function names) so that it does functionally the same but using many different polynomial pairs (where a pair is the quadratic and a linear polynomial zx+(zx+y) whose resultant = 0 mod N), rather then just one like the uploaded PoC currently does.
2. To achieve step one, I first need to study the square root over a finite field step in NFS_Solve() and figure out the number theory to see how that would work with what I was already doing in my older iterations (See link to old version without NFS above).


Once I figure out 2, I should be able to fairly quickly implement 1. Anyway.. its a marathon, not a sprint. I know I'm nearly there now. All this endless suffering and social isolation for years will have been worth it soon.

Note: You can solve the example from the paper (N=4387) using python3 run_qs.py -key 4387 -base 30 -debug 1 -lin_size 1000 -quad_size 1, what I'll be doing next is to just trace this example and see how both the rational and algebraic square roots are constructed (well, rational side is easy, I know that alreayd, its mainly the algebraic square root being a little elusive). I'll Just pen and paper trace it the hard way, then I should be very quickly able to connect the dots with everything I've already figured out... now if I had math education I could probably just read some math papers and figure it out.. but I'm playing math research on the hardest difficulty, ahahaha.

Update: Just got back from running. Depression man. Everything is just shit lately. I havnt seen my friends in years. People dont even reach out anymore. I've literally just spent years alone in a room grinding math and running. And these persistent gut issues wont go away. They started last summer and my running performance has also been degrading. I can still get my miles in, but my pace is just dropping off and I dont know if its because of maybe too much running volume and poor recovery or if I have some serious health issues with these gut problems. Lately its been getting really bad. I'm trying to figure out if its something I'm eating or something. Just stress man. It better not be cancer or some shit. Imagine if I had to live my last years like this. I would die a very bitter person.

Update: Alright. I've set everything up for tomorrow. My goal for tomorrow is to figure out whats happening with the square root over a finite field code in create_solution() and relate that to my own work (hopefully). I'm hoping that wont take more then a day. We'll see tomorrow. Once I have assimilated that, I should be able to finish my work at last (hopefully). I would really like to see the end of this now... I'm mentally going to really bad places. I'm being treated like I'm mentally challenged or something.. not qualified to be employed, better to be kept out of society. I miss my former manager. In 2019 he came all the way to belgium to meet me and be my friend and eventually hire. he freed me from this shit country.. for a little while. And then microsoft went after him bc he gave me a promotion once. I'm sure they will tell all types of lies.. but if you want the real story you should ask the people who actually kne me. Not the people who only saw or heard from me when I was at the worst points in my life (like in 2023).

Update: I slept like shit. But started doing analysis. That NFS implementation, its so close to what I was doing in my non-NFS version (see link above)... they even factorize the quadratic coefficient (but they call it the leading coefficient). Its the same! Which means there also is a good chance I can get it working with multiple quadratic polynomials instead of just one like public variants of the NFS algorithm do. I'll continue taking this apart. Whoever uploaded that NFS implementation last summer, thanks, great timing there. I guess worst case I probably would have taken apart CADO-NFS and figured it out that way... but having this python implementation saves me literal weeks, if not months.

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
