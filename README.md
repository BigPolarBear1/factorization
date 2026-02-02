Disclaimer: At no point did AI contribute anything to this research project. Copilot can't even calculate basic congruences without making mistakes. People who think AI can do novel math research are delusional. I did not use AI for the paper, and I did not use AI for any of the code. I attacked this problem, from the very basics 3 years ago, slowly burrowing deeper and deeper. I know people will say lies because they will never let a person like me win. Not when it comes to something as important as this. Because that would literally go against all the bullshit propaganda the far right has been spreading for years. Also, this could have been avoided. Its a disaster of your own making. 

#### To run from folder "Improved_Sieving":</br></br>

Note: Experimental WORK IN PROGRESS.</br>
To build: python3 setup.py build_ext --inplace</br>
To run:  python3 run_qs.py -keysize 60 -base 500 -debug 1 -lin_size 100_000 -quad_size 1</br></br>

Update: I added a linear congruence that will calculate linear coefficients for polynomial values that are multiples of the modulus. So the only thing that is left missing in the PoC is using this hashmap so that the x+y multiple we are generating has the correct factorization (it will add the factorization of s2 to it, the result of the linear congruence). And right now it is using the factorization of s2%n (line 948), but this can be s2-n or s2-n\*2 or s2-n\*3 and so on.. and if any of those factorize, we can use it, and even more, we can use the hashmap to get the factorization we want there.. but that is what i'll begin implementing now. I dont know man, it is a surreal feeling, knowing the scope of the disaster you are creating.. yet, all I wanted was for a single person to take my work serious. But nobody does... so I have no choice. Good luck. Will upload in the coming days. 

Update: I was thinking, if I delete that line s2%=n and s also factorizes over the factor base.. then we get near identical smooth candidates. Then its just a matter of getting that s value to factorize. Hmm. I can figure that out for sure. That shouldnt be in the realm of the impossible. And in fact.... we can use much larger factor basis to get that s to factorize.. since once it factorizes, if my brain at 12:20am is understanding this correctly, should spit out near identical smooth candidates, so we wouldnt be constrained by using a very large factor base, since we'll have a super big chance to succeed with much fewer smooths. Something like that. Man, if that is it, people are going to be freaking out, Ill probably get droned striked. Hahahaha. If I die tonight, it was the american cowards. I'm not afraid. I did everything right. I tried to avoid this outcome. The fbi knew I was working on this way back in january 2024. 
I'm not a cyber terrorist, but I'm also not a cunt if you disrespect me like this, so be it, I press on.

Update: Well almost correct. So actually you need to subtract some multiple N from "s" so that the result has "enough" factors in common with s. Bc we need unique factorizations. Ok ok. I'm seeing how to do it now. Easy.

Update: You know, this will actually be pretty great. Tomorrow my work will come to a conclusion. Its pretty easy now. Man, the amount of mental anguish these last 3 years.. just because people are morrons. Imagine if someone, who has proven to be really persistent and has done plenty of creative research in the past, imagine someone like that proclaiming they are working on factorization and you decide its a good idea to treat that person like complete trash. Thats not strategy, thats being a fucking morron. And fuck belgium most of all. I come back to this? Years of unemployment and getting harassed by the courts. lol. Fuck this. 

Update: Just did some quick tests. And yes, this 100% works. Go to hell. People would have known. More experience math folks would have known had any of them looked at my work. Everything, from NOW until I begin setting factorization is records, its just implementation details. I have figured out all the pieces of the puzzle. Now it is optimization work and implementation details. 

Update: Added another loop at line 960 in construct_interval() .. now this s2 variable at line 962, this is the only thing introducing new factors. Hence we need to minimize introducing new factors here by sieving this and using a modulus. Its slow as hell right now bc its just bruteforce iterating those s2 values. I could just setup a sieve interval, easy. But I'm also thinking if I can just straight up calculate this using the hashmap. Because that would probably be more effective. Anyway.. just implementation details. All the math is there. Everything is there. This is as close to disclosure as you can get. And people would know. Now that I have worked out the exact strategy.. that crappy uploaded PoC is quickly going to become faster and faster. Plus the whole thing is that it is going to generate smooths with similar factorizations... hence we can use much larger factor bases. Whatever. I will keep grinding, and if nobody takes my work serious, then so be it. All that tells me is that this world deserves this disaster, if this is how you treat people.

There is probably also some simplifications that can be done in the PoC as there is quite a bit of steps involved. let me have a look. But everything is there. All the math is there. It proves everything. And faster versions will come online now in the coming days as I simplify and improve everything.

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
