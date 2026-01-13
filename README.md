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
To run:  python3 run_qs.py -keysize 20 -base 50 -debug 1 -lin_size 1000 -quad_size 1</br></br>

Update: Doing some more work today. I realized that I can use different quadratic coefficients as well... as long as the linear coefficient is the same. Everything else can be different. Let me see how to sieve this now.. all the math is there.. I just need to develop a proper sieving algorithm now.

UPDATE: After pondering every angle of attack for a while.. I think the best way to use these full quadratics is to sieve the discriminant... we want to arrive at quadratics where 2zx = y... because if that factors over the factor base, then also the discriminant factors over the factor base. So  ideally I now need to figure out how do I combine multiple quadratics where 2zx != y to one where 2zx = y. Using these full quadratics as a way to build smooths, I think is the correct appraoch. Its really the only angle of attack that I can see working. Let me figure that out now... almost there... almost. 

UPDATE: So let me begin figuring out how to use these quadratics to basically "sub-sieve" the discriminant. I know the discriminant will produce the same output as the quadratic if 2zx = y. So I need to figure out the math to combine quadratics so that the combination satifies the condition that 2zx = y and we produce a smooth. That would improve the quadratic sieve algorithm... sieving with the full quadratic gives us much more control.

I think I know how to "sub-sieve" that discriminant using quadratics with non-zero linear coefficients and where 2zx != y .. let me try to write something tomorrow.. 

UPDATE: Slow day today. I changed the PoC quickly to use y<sub>0</sub> instead of y<sub>1</sub> ... because the signs were giving me an headache. The modulus code should also be working correctly. And in addition I'm not sieving mod N anymore... we just use that quadratic coefficient to sieve multiples of N instead. Sieving this way also gives me a little less of a headache. I guess in the future I can always change it bac to mod N... but this will be a little easier to get the sieving portion implemented without the added complexity of working mod N.
So now to implement the sieving algorithm next... Once 2zx = y we know we'll have a smooth aswell with the discriminant (where the constant is N without an offset).

Anyway... going to go for a run.. fucking stressed man. I know I'm right about my math. I know what I got. And people must know... not everyone in infosec can be a complete morron. Ofcourse they think people like me are "DEI" and blah blah blah.. people are so fucking brain washed, they woudln't even believe it if it was staring them right in their eyes. Anyway, factorization will come to this world, even if it is the last thing I do in my life.
