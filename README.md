Disclaimer: At no point did AI contribute anything to this research project. Copilot can't even calculate basic congruences without making mistakes. People who think AI can do novel math research are delusional. I did not use AI for the paper, and I did not use AI for any of the code. I attacked this problem, from the very basics 3 years ago, slowly burrowing deeper and deeper. I know people will say lies because they will never let a person like me win. Not when it comes to something as important as this. Because that would literally go against all the bullshit propaganda the far right has been spreading for years. Also, this could have been avoided. Its a disaster of your own making. 

#### To run from folder "Improved_Sieving":</br></br>

Note: Experimental WORK IN PROGRESS.</br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 70 -base 500 -debug 1 -lin_size 100_000 -quad_size 1</br></br>

Update: Bah, I'm going to switch it up. Just two small roots. Then we only need to worry about the factorization of the polynomial value. And we just work with the precalculated hashmap. Just pull smooths out of it. I know it can be done, i've done it previously. Plus I could probably try and find smooths with similar factorization this way. 

So just calculate an enormous factor base and then be strategic about the smooth factorization we try to find, so we can succeed with very few smooths. I know this works. I think this is probably the best way to approach it.

Update: The more I think about this, the more I realize, that yes, finding smooths with a specific factorization is a lot more feasible now. Let me write a function for this tomorrow. And then the size of the factor base matters a whole lot less.

Update: I was messing with some numbers. Observe this sequence:

175^2+1\*175-4387\*7 = 91</br>
175^2+26\*175-4387\*8 = 79</br>
175^2+51\*175-4387\*9 = 67</br>
175^2+76\*175-4387\*10 = 55</br>
175^2+101\*175-4387\*11 = 43</br>
175^2+126\*175-4387\*12 = 31</br>

So using k (see chapter 8 in the paper) we can achieve some really good control over the polynomial values. Oh yea, I think I have an idea now. 

Update: Oh yea. That k variable is very useful. Damn. Should have seen that earlier. The thing is, this will actually let us sieve the polynomial value and x+y properly. Damnit. I got to move fast on this tomorrow. 

Update: Awful day. Been having gut problems for months now. And this last week its really bad. Hard to focus on work like this. I'm wondering if its the fructose in this sports drink mix. I think I'm just going to radically change my diet. No more sports drinks. Just salt pills for electrolytes and (real) honey for carbs when I run longer distances. I'm starting to think the fitness industry is a scam. Anyway.. for the PoC.. you need to use this k variable (multiples of N) to change the polynomial value. This way you can also increase the linear coefficient in such a way that x+y keeps a consistent factorization. This was exactly what I needed. I'll upload a more mature PoC this weekend. 

Gays, keep that MQ drone flying at altitude bc the sound of its propeller is very distinct. Noobs. You're welcome. 

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
