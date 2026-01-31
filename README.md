Disclaimer: At no point did AI contribute anything to this research project. Copilot can't even calculate basic congruences without making mistakes. People who think AI can do novel math research are delusional. I did not use AI for the paper, and I did not use AI for any of the code. I attacked this problem, from the very basics 3 years ago, slowly burrowing deeper and deeper. I know people will say lies because they will never let a person like me win. Not when it comes to something as important as this. Because that would literally go against all the bullshit propaganda the far right has been spreading for years. Also, this could have been avoided. Its a disaster of your own making. 

#### To run from folder "Improved_Sieving":</br></br>

Note: Experimental WORK IN PROGRESS.</br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 70 -base 500 -debug 1 -lin_size 100_000 -quad_size 1</br></br>

Update: Quickly added support for moduli in the polynomial value.  Next, for every root and linear coefficient combination where also x+y factorizes.. we should then also enumerate the k variable now. This will allow us to generate many near identical smooth candidates really quickly. 

So remember, the final smooth factorization is made up from x, x+y and the polynomial value (ignoring quadratic coefficient z) .... so if we iterate k (multiples of N) we can add to the linear coefficient a value, such that the factorization of x+y stays largely the same and we then use the k variable to shrink the polynomial (the amount by which it shrinks is dependent on N%x and the increase in the linear coefficient). So thats going to garantuee that all 3 parts, x,x+y and the polynomial value have consistent factorizations and we can also partially garantuee the factorization of the polynomial value and even use the hashmap to aid with that if we precalculate it for multiple quadratic coefficient (which just encodes that k variable).

Next version I upload will be a first rough draft of that.. so spit out smooth candidates with similar factorization quicly by iterating the k variable in the quadratic zx^2+yx-Nk

Update: I had an idea just now. So the easiest way to approach it, is by doing something like this:

3^2+10\*3-4387  = -4,348				x+y = 13</br>
3^2+(10+13\*10)\*3-4387 = -3958			x+y = 13\*11</br>
3^2+(10+13\*120)\*3-4387 = 332			x+y = 13\*11\*11</br>
3^2+(10+13\*1330)\*3-4387\*12 = -735		x+y = 13\*11\*11\*11</br>
3^2+(10+13\*14640)\*3-4387*130 = 689		x+y = 13\*11\*11\*11\*11</br>

So a setup like this will just repeadtly keep adding the same factor (11 here) to x+y. We can keep repeating this until the end of the universe basically. And since we dont care about the factorization of the k variable .. we can also keep increasing that until the end of the universe to keep the polynomial values small. The only thing missing from this setup is garantueeing that the polynomial values have similar factorization, but we already know how to do that using that hashmap. So that way.. we're going to end up with a bunch of smooths with similar factorization... like... extremely similar factorization. Which means we'll have an enourmous chance of finishing early. 



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
