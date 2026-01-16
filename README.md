Disclaimer: At no point did AI contribute anything to this research project. Copilot can't even calculate basic congruences without making mistakes. People who think AI can do novel math research are delusional. I did not use AI for the paper, and I did not use AI for any of the code. I attacked this problem, from the very basics 3 years ago, slowly burrowing deeper and deeper. I know people will say lies because they will never let a person like me win. Not when it comes to something as important as this. Because that would literally go against all the bullshit propaganda the far right has been spreading for years. Also, this could have been avoided. Its a disaster of your own making. 

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

#### To run from folder "nfs_variant" (WIP):</br></br>

Note: Experimental WORK IN PROGRESS.</br>
To build: python3 setup.py build_ext --inplace</br>
To run:   python3 run_qs.py -keysize 14 -base 50 -debug 1 -lin_size 1000 -quad_size 1</br></br>

Update: Oops, I've re-uploaded what I had yesterday, what I was trying earlier today wasnt quite the right direction. I just had a realization though. So the discriminant with an offset in the constant..the easiest way to get rid of that is to have that offset be the modulus. That works fine if the offset/modulus is larger then N. That would probably work... hmm. Let me explore that idea some more.


Update: I GOT IT!

So if N = 4387: 

148^2-4\*(4387+248) = 58^2 and 66^2 = 58^2 mod 248</br>
148^2-4\*(4387+189) = 60^2 and 66^2 = 60^2 mod 189</br>
148^2-4\*(4387+128) = 62^2 and 66^2 = 62^2 mod 128</br>
148^2-4\*(4387+65) = 64^2  and 66^2 = 64^2 mod 65</br>
148^2-4\*(4387) = 66^2     </br>

Another example:

1068^2-4\*(4387\*65-35) = 12^2 and 2^2 = 12^2 mod 35</br>
1068^2-4\*(4387\*65-24) = 10^2 and 2^2 = 10^2 mod 24</br>
1068^2-4\*(4387\*65-15) = 8^2 and 2^2 = 8^2 mod 15</br>
1068^2-4\*(4387\*65-8) = 6^2 and 2^2 = 6^2 mod 8</br>
1068^2-4\*(4387\*65-3) = 4^2 and 2^2 = 4^2 mod 3</br>
1068^2-4\*(4387\*65) = 2^2</br>

THATS MY LINK TO NFS COMPLETED! I see how to do it now!!!! I will notify my Chinese friends so they can destroy america, ahahahahahahaha. 

Update: Anyway, I'll upload a PoC tomorrow that abuses this. Because thats how we can use square finding and linear algebra together with what I'm doing now. Stressed. I keep having these really bad gut problems for months now. I'm actually starting to think I may have a serious health issue. Probably would have seen a doctor already if my life wasn't so absolutely shit. Either way, I stopped taking creatine a few days ago.. I'm hoping these gut issues were just related to water retention or something.. we'll see.

Update: Digging a better deeper into what I'm seeing.. I updated the above number sequence... so that is a very useful pattern. Let me move toward completing a PoC now.

Update: I'll finish the PoC this weekend.. I've added what I've found to the paper also. Just these gut issues are stressing me out so much. Legit starting to think I may have fucking cancer or some shit. And tomorrow I have to run 30k to and from therapy again. Maybe I'll complain to my therapist about these fucking health issues lol. I dont know, if I do have some serious issue like cancer, just let me die man. I'm fucking tired of this life anyway.

Update: Today with that offset in the discriminant figured out. I know now this was the final thing. What I've found is a much better way to do NFS. And its so absurd that everyone ignores what I'm doing. I dont even feel happy today. Despite knowing I've done it. A full PoC is a matter of days now. The road is open for me to finish this now. I'll email my Chinese friends. Fuck the west for treating me like this. Wouldn't even be suprised if these western fuckers gave me cancer or some shit and hoped I would die before finishing my work.

You know what would be funny, imagine if some highschool drop-out bug hunter, who 3 years prior couldn't even do basic algebra suddenly started completely destroying the factorization problem. I wonder what people would think. I guess people would be seriously spooked. But maybe I know people would be seriously spooked, maybe I enjoy spooking people, because my life is shit and making people scared is the only joy I have left. Anyway guys, One foot in the present, one in the future, we march on. People about to get spooked a whole lot more. Hahaahaha. Who knows what is going on, anything could be going on. 
