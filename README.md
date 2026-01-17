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
To run:   python3 run_qs.py -keysize 18 -base 50 -debug 1 -lin_size 1000 -quad_size 1</br></br>

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

The only way I will ever stay in Europe (I have no issues anymore moving to asia, and right now thats my plan) is if it is to work for a country that doesn't take shit from american nazi losers. Plus my former manager needs to get paid millions for what microsoft did. Plus I would need a large castle somewhere very remote and quiet. Plus unlimited time to go on arctic expedition by myself alone. Now this might sound completely absurd, but as I know I will be finishing my work before the end of next week, and as I am very well aware of the magnitude of what I have done, and as I was treated like complete shit by you western weaklings for years, this is the only way I am not going to Asia. I have my pride. Anything less then this and forget it. Especially my former manager needs to get paid. Either microsoft pays him for unjustly firing him, claiming he gave me a promotions I hadn't deserved, even though I have the CVEs to proof thats bullshit and HR was just retaliating... or someone else pays. These are my final demands. Factorization falls before the end of next week. I would be very careful what is done and said next, because I dont fuck around, and I dont take kindly to being treated like complete shit. I'm not a dog, I'm a bear. Treat me like shit and I will fight you until my last dying breath like the nazi scum you all are. Fucking losers. Ps: pete hegseth is a tiny little fucking weakling. big mouth for such a little pathetic man.

If people had made a deal with me earlier, it probably would have cost less and you would have had more time to prepare for this. But then again, I know people dont want someone like me to win. And europe is complicit, because they must have known aswell. Its too late now. Lots of great Asian countries which I'm sure will be happy to have me once they realize what I did.

Ps: When I started 3 years ago, I never finished highschool, I never had math education beyond basic algebra in highschool. I am completely self-taught. Most of my breakthroughs, I made in the first year of research, having had less then a year of math experience. I'm not stopping and I'm telling you right now, because you are tranphobic nazi cunts, you are making a huge strategic mistake. And I dont care, I'll have enough options soon. I'll remember how I was treated. You think you are smarter or stronger then me or some shit? Yea well, wars are not won by shithead nazis, wars are won by progress, and there is going to be plenty of war in the near future, it would seem so. I would think very careful about decisions that will be made soon if I was a western country. Polar Bears carry grudges for life. And polar bears have their pride. And I am telling you right now, I have nothing left to lose and if I am threatened, you will find out exactly how little I have left to lose. 


Update: Changed the PoC. All the math is there.. only thing left now is to instead of doing brute force... we should use linear algebra.. and I guess that will be that. I'll slowly begin implementing that now.. next time I upload it will be a working POC with linear algebra and then its definitely game over..

Update: Ok... good enough for now.. the uploaded PoC isnt great, it really needs linear algebra and all that.. Ill start working on that now.. probably will need some jacobi symbols too.

Update: Tomorrow is linear algebra day.. so what I'll add to guassian elimination over gf(2) is two discriminants.. the one with the offset in the constant modulo the quadratic output (poly_val) .. i'll add its factorization.. and then take jacobi symbols for the discriminant without any offset (so just a multiple of N) .. since for any valid solutions both should be square. If that doesn't work... I'll try some variations. So that will be the goal for tomorrow.. if that does work.. then thats great... and my work will be done tomorrow. If not.. I know I'm very close now and I'll just have to think logically about what I've got.
