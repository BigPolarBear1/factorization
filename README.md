Disclaimer: At no point did AI contribute anything to this research project. Copilot can't even calculate basic congruences without making mistakes. People who think AI can do novel math research are delusional. I did not use AI for the paper, and I did not use AI for any of the code. I attacked this problem, from the very basics 3 years ago, slowly burrowing deeper and deeper. I know people will say lies because they will never let a person like me win. Not when it comes to something as important as this. Because that would literally go against all the bullshit propaganda the far right has been spreading for years. Also, this could have been avoided. Its a disaster of your own making. 

#### To run from folder "Improved_Sieving":</br></br>

Note: Experimental WORK IN PROGRESS.</br>
To build: python3 setup.py build_ext --inplace</br>
To run:   python3 run_qs.py -keysize 30 -base 30 -debug 1 -lin_size 20 -quad_size 10</br></br>

Just run it on 30 bit using the command above. It should find smooths whose factorization is now mostly determined by x+y. Which in this PoC we are just setting to be lin_size, hence most smooth factorizations will have primes between 1 and 20 if using lin_size 20. This enables us to succeed with fairly few smooths. ofcourse choosing which x+y to check is an area for improvement, as once we find a smooth, we should probably set x+y to be its odd exponent factors and check for those. But anyway, this is a first extremely rough draft. Now we begin optimizing and improving. So in theory, the size of the factor base does not matter at all anymore now. Hence speedingup these calculations is what I'll focus next on. Once that is as good as it can be, I'll improve the selection of which x+y values to check for based on smooths we have found so far. 

This should have been my last dip in performance... now upwards and onwards... 

In addition there also probably is some linear algebra tricks to help select x+y values and k (multiples of N) values... but first things first.. building that factor base quickly..

Update: I'll add hensel's lifting lemma to solve_roots() tomorrow. Thats my goal for tomorrow (thursday). That should see a dramatic increase in computing the factor base. I havnt yet implemented hensel's lifting lemma for polynomials with non-zero linear coefficient... but it shouldn't be very hard, plus there is code examples online that do just that. I think once that is implemented, that part will probably be as optimized as can be math wise. Then after that I'll see coding wise what can be done there. Oh and also computing the results for different k values, I can probably derive those results like I do in CUDA_QS_Variant.. I dont need to do everything all over again there.

Update: Added hensel's lifting and removed the hashmap in favor of 1d lists. Building hmap is still very slow. But it doesn't matter that much for now. Let me now make sure we select linear coefficients such that the polynomial value is as small as possible. Because there is no reason to do it random like we do now, when we can use it to generate small polynomial values..

Update: Hmm. Let me do some thinking. The thing is, that k value (multiple of N for the constant), we dont care about its factorization. Hence we can use that to add or subtract N as much as we want, until eternity. Thats a useful tool. I should center my strategy around that. Then whats left is selecting an x and x+y, with a factorization that is know, which has enough hits in the precalculated datastructure (hmap) to garantuee a smooth. Let me do some thinking. Because this x and x+y can be huge values, we dont care, as long as they factorize.. as long for some k value, they have enough hits in hmap. And this k value can be anything and the factorization doesn't matter. Common. Its on the tip of my tongue now. Give me a day or two to work this out.

Update: You know, what I should do is find an x+y value, which has a large occurance in hmap, regardless of the k value. Because we can increase x to get the polynomial value closer to 0 if k is very large. Or visa-versa find an x with a large occurance and then increase x+y, which does give us more control. But either way, thats definitely a better approach then just the bruteforce way the code is doing and only taking into consideration k value < quad_size.. thats not really playing into the strengths of all the number theory we're using. Anyway, going to run in the dark in the woods. 

Update: If we have a singular root, then lifting seems to fail. Its not an urgent fix, it just means we'll be missing some solutions here and there. Anyway.. I'm an idiot, I know how to approach this PoC.

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
