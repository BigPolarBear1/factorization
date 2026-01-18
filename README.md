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
To run:   python3 run_qs.py -keysize 30 -base 50 -debug 1 -lin_size 1000 -quad_size 1</br></br>

Update: Almost there. Not quite yet ready to implement linear algebra. I need to spent a few more days taking apart what I've found so far. 
There is something identical going on as to what number field sieve is doing. But its fairly complicated. Pretty I have all the pieces now.. I just need a few focused days studying it so I can implement the linear algebra.

Update: Rolled back to a few days ago. I see the link with number field sieve, but I just got to dig in now. Let me figure this out.

Update: I changed it again to use y<sub>1</sub> instead of y<sub>0</sub> ... now I'm going to go ahead and figure out how to pull off linear algebra...

I guess if we have something like this:</br>

66^2+4*(4387+1749) = 170^2 = 148^2 mod 1749 </br>
66^2+4*(4387+1248) = 164^2 = 148^2 mod 1248 </br>
66^2+4*(4387+1085) = 162^2 = 148^2 mod 1085 </br>
66^2+4\*(4387+453) = 154^2 = 148^2 mod 453</br>
66^2+4\*(4387+300) = 152^2 = 148^2 mod 300</br>
66^2+4\*(4387+149) = 150^2 = 148^2 mod 149</br>
66^2+4\*4387 = 148^2</br></br>

453-153=300-151=149-149=0</br></br>

Reducing the square and constant eventually hits a multiple of N as constant.</br></br>

But this is not always the case:</br></br>

2^2+4\*(4387\*4379+15) = 8766^2</br>
2^2+4\*(4387\*4379-8750) = 8764^2</br>
2^2+4\*(4387\*4379-17513) = 8762^2</br>
2^2+4\*(4387\*4379-26274) = 8760^2</br>
etc..</br></br>

Ofcourse then we would need to find a square such that that offset to N is divisible by N. But that is the same as solving a quadratic congruence, which we cant do without knowledge of the factorization of N and in addition this would basically change that k variable.. aka the quadratic coefficient. 
However: 66^2+4*(4387+1749) = 170^2 = 148^2 mod 1749 

Because 1749 mod 11 = 0, so we know that it is atleast aligned to a multiple of N mod 11. And ofcourse 170 = 148 mod 11. 
So if we use the factorization of that constant as the modulus..  then that linear coefficient should also correctly map to the quadratic coefficient (z\*k)

I see how it works. I'm struggling hard with depression. Its so incredibly hard to think clearly this way. I really need vacation or a break, but how? I have no future, I have no income, no employment, I havnt seen my friends in years since what microsoft did. The math isn't hard. Its easy. Its just finding my focus with all the other bullshit thats the hardest. If my life wasnt so shit, I would solve this within a day.
