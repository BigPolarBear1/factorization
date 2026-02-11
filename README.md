Disclaimer: At no point did AI contribute anything to this research project. Copilot can't even calculate basic congruences without making mistakes. People who think AI can do novel math research are delusional. I did not use AI for the paper, and I did not use AI for any of the code. I attacked this problem, from the very basics 3 years ago, slowly burrowing deeper and deeper. I know people will say lies because they will never let a person like me win. Not when it comes to something as important as this. Because that would literally go against all the bullshit propaganda the far right has been spreading for years. Also, this could have been avoided. Its a disaster of your own making. 

#### To run from folder "Improved_Sieving" (NFS variant.. but in the process of adding the NFS part now that I figured out how to generate irreducable polynomials):</br></br>

Note: Experimental WORK IN PROGRESS.</br>
To build: python3 setup.py build_ext --inplace</br>
To run:  python3 run_qs.py -keysize 60 -base 500 -debug 1 -lin_size 10_000 -quad_size 10</br></br>

Update: I was grinding some more numbers.

For example if we have 20^2+33\*20-4387 = -3327
We can rewrite this as: 20^2+33\*20+3327 = 4387
Now the constant in the discriminant is calculated like this: 20^2+33\*20 = 1060

And we see 33^2+4\*(4387-3327) = 73^2 where 4387-3327 = 1060 (obviously). 

Now a polynomial of this format: x^2+33\*x+3327 is irreducible and we know x=20 produces N (4387).
Since we have an irreducible polynomial we can apply number field sieve. 
And as a matter of fact, this should also help predicting the factorization of polynomial values. 

Ergh. 

I hate myself right now. I should have seen this a year ago already. I absolutely hate myself right now. Time to get moving on this.

Update: Quickly added it to the bottom of chapter 8 in the paper. I already have some python code for the NFS algorithm.. so I'm going to move fast on this. All I wanted was for someone to take me serious. Nobody did. Even now, just silence. I cant find a job. I cant get anywhere in life. I'm going to finish my work now.. and when all this is over, I hope you people can sleep at night.

You know, since I have a fairly good grasp of all the math I've spent 3 years figuring out.. I am also fairly confident that as I go into my final chapter of my research project, I can probably get the NFS algorithm to work using multiple quadratic polynomials.. so that we can switch polynomials if finding smooths becomes an issue. Ah, it will be easy. This will basically be my victory lap. Alright gays. You all knew this was coming, I hope you're ready now. Next time you see activity here, it will be with NFS added to the existing code. Or maybe I will give my non-nato friends a heads-up... it's fifty-fifty depending on my mood. I'm not in a great mood lately. Good luck shitheads.

ps: Even if, on the offchance, adding NFS to my code now doesn't result in a breakthrough algorithm... I started with 0 math knowledge 3 years ago. I dropped out of highschool and I never had math education beyond that. With the link to NFS finally made.. I have gone through the whole process now.. there is few people who understand the underlying structure of the factorization problem as well as I do. Sure, some people may understand QS, or even NFS... but they havnt spent years taking everything apart, studying the patterns, making sense of that underlying structure. I have. Even on the slim chance, I dont achieve my breakthrough before the end of the month, there is few people in the world better positioned then I am to end up making a breakthrough eventually. And you are all morrons. Imagine not taking advantage of that with war erupting everywhere. Instead you are all like deer staring into the fucking headlights. And I dont care anymore, because you people have made me absolutely miserable and perhaps disaster is what you people deserve.

Time to get started.. I'll just plug in the NFS algorithm very roughly, just use my code to generate irreducable polynomials and after that I'll start simplifying and streamlining everything. A problem I had before, was when looking at NFS, I found it hard to relate it to my own work because I just didn't had made enough progress yet with my own work. But now I have.. so I can finally put NFS in the context of my own work and hopefully that is now going to allow for a breakthrough. We'll see. I'm optimistic. And even if I dont have a breakthrough.. its another step forward and the breakthrough will eventually come. People spent many years learning math, I've only been doing it for 3 years. I literally couldnt do basic algebra 3 years ago. I'm not plateau-ing anytime soon skill-wise. If people underrestimate me.. then thats their problem. Like microsoft trying to portray me as incompetent and firing my manager. When I have my breakthrough, I'm going to make sure everyone knows about the shit microsoft did. And I'm sure they will tell lies or use out-of-context shit to smear me, I dont care, the people who knew me there, they know what really happened.

Update: I've begun. Will take a few days, I have to study a few topics and read some papers before I integrate my work into NFS. Apparently Kleinjung's polynomial selection method is the best method for selecting polynomials... but I kind of want to take apart how that works first and see if there's anything I can relate to my work there. I understand the classical approach of using just one polynomial. But using two of them like with Kleinjung's method is interesting, so I need to figure out why that works.

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
