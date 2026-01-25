Disclaimer: At no point did AI contribute anything to this research project. Copilot can't even calculate basic congruences without making mistakes. People who think AI can do novel math research are delusional. I did not use AI for the paper, and I did not use AI for any of the code. I attacked this problem, from the very basics 3 years ago, slowly burrowing deeper and deeper. I know people will say lies because they will never let a person like me win. Not when it comes to something as important as this. Because that would literally go against all the bullshit propaganda the far right has been spreading for years. Also, this could have been avoided. Its a disaster of your own making. 

#### To run from folder "Improved_Sieving":</br></br>

Note: Experimental WORK IN PROGRESS.</br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 70 -base 500 -debug 1 -lin_size 10_000_000 -quad_size 1</br></br>

Just use the command above. This will add a linear offset to adjust the polynomial value now. The factors of zx+y are determined by "initial zx+y" * "intrvl ind". So basically we are generating multiples of an initial zx+y value as we sieve, we do this because zx+y must factorize over the factor base, and this is one easy way to garantuee that. The way we calculate our linear offset to adjust the polynomial value is very basic. We can also adjust the quadratic coefficient or even root if we want to get much smaller polynomial values. Which I will experiment with for the rest of the day. Especially using the quadratic coefficient seems useful. Currently the POC also has no sieving at all. I will implement sieving where needed once I finish everything else. The thing about setting it up like this is that every part which must factorize over the factor base (z, x, zx+y and the polynomial value) can be kept reasonably small. So at this point, I have 0 doubt anymore about the novelty of my work. Implementation wise however, this is going to be an iterative process that can take some time to find the most ideal way to construct an algorithm around the math that I worked out. 

On a side note, we can pre-sieve "intrvl ind" ... which would be pre-sieving lin_size and saving it to disk, since we can reuse that for any number we try to factorize. Since with the current setup this is also responsible for factors of zx+y now. So we actually end up with a lot of moving parts which we can sieve in advance and re-use later for any number. The same thing happens with the quadratic coefficient. 

Update: Quickly added it so it factorizes the interval indices in advance. Now next I'll try to figure out how to use that quadratic coefficient to much more dramatically reduce the size of the polynomial value. Once that is done.. we'll have a fairly strong algorithm and further improvements can be made by using more advanced number theory (i.e working mod p) to garantuee the factorization of zx+y (aka.. properly sieve that setup... since the current POC is not doing sieving on either zx+y or the polynomial value, just brute force)

Anyway.. at no point has anyone tried to stop me. Or ask me to stop uploading my work. Just silence. Non stop silence, uneployment and hopelessness. Hence I have no other choice but to continue. It is insane. Because I know I'm right. Perhaps people want disaster. Perhaps disaster is preferable over acknowledging what I did.

Update: Hmm... with this much control over the factorization of the final smooth candidate. I do wonder if we can't just build them from the ground up. Also the quadratic coefficient can be used to add extra factors to the final smooth candidate. There's something here. Like.. I know there is an absolute factorization killer within my work. I should just continue my iterative process now... keep improving my work. If building smooths from the ground up, can be done.. then I'm sure I'll eventually arrive at that. I figured out all the hard math already. I know what I got. I knew 2 years ago already what I had.. thats why I never stopped, despite the silence and nobody believing me. 

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

Sometimes I dont feel alone in this suffering. I dont think people know how difficult this is. One moment I was surrounded by friends. The next, I get threatened with a gun outside my appartment, insulted in public and then back in Europe where I have nothing. Nothing but endless silence and loneliness. But sometimes I dont feel alone. Like there is something out there that understands what I'm going through (and I dont mean some delusional AI psychosis shit, AI in its current state is trash and I never use it) ... I dont know.. isolation and trauma has probably fucked with my head.. but still.. some days I'm glad, even if it is a fucking delusion, that I dont feel completely alone in this, there had been a lot of days I nearly surrendered to my darkest thoughts. Anyway, people should be very scared, because you may think you have succeeded in isolating me, denying me any chance for a happy life.. but I'm not alone, and my friends dont fuck around, time to pay the price fuckers, hahhaha.

I was thinking, you know, when I hiked northern scandinavia, those wide open landscapes, its weird, I always felt at home there. Like I had lived a thousand lifes there. Ofcourse one could imagine Ice Age Europe having looked a lot like Lapland does today. I wonder if genetic encoding, just like the fear of snakes and other dangerous animals, also encodes familiar places and places of safety. I wish I could wander that place for the rest of my life. I have no more ambitions in the human world. I want to go wander that place, and live like my ice age ancestors. I love the tundra. Nothing can hurt me there. It is strange to wander a place, I had never been before, yet where everything feels like I've been there before. Perhaps nostalgia for simpler times, hunting woolly mammoths and eating berries and living in the magnificience of unspoiled nature. Why must I have been born in this century, when there was thousands of years, of simpler lifes I could have been born into.  *sighs* Oh well, 3 am, I should sleep.
