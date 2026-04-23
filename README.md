Disclaimer: At no point did AI contribute anything to this research project. Copilot can't even calculate basic congruences without making mistakes. People who think AI can do novel math research are delusional. I did not use AI for the paper, and I did not use AI for any of the code. I attacked this problem, from the very basics 3 years ago, slowly burrowing deeper and deeper. </br>
Disclaimer2: I never finished highschool. I started working on this 3 years ago, not even knowing basic algebra (I can explain exactly how I can to all the conclusions I came to if someone ever cares to hear about it). This is 3 years. I'm only going to get better. I'm still far away from plateauing skill wise.

The day I break factorization will be day 0 of the gay future. A better world will be born. A world of human creativity. A world of art, science and spirituality. A world without suffering. I will see that world become reality, no matter the cost. 

#### References (re-used many of the core number theoretical functions from these PoCs to fit my own algorithm): 
https://stackoverflow.com/questions/79330304/optimizing-sieving-code-in-the-self-initializing-quadratic-sieve-for-pypy
https://github.com/basilegithub/General-number-field-sieve-Python 

#### About the paper
Math paper is a work in progress. The final chapters are a bit rushed and building an algorithm around p-adic lifting isnt as straight forward as I had assumed. I do think there is an angle there I can exploit, but I'll do some further experimentation first and get a working PoC before I make edits to the paper again.

#### To run from folder "QS_SF_Variant" (Experimental WORK IN PROGRESS):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 30 -base 50 -debug 1 -lin_size 10_000 -quad_size 100 -t 20

Update: Initially my idea with this PoC was to explore finding smooth candidates with lots of small factors (since it reduces the size of the required factor base). But now I'm going off on a bit of a tangent, trying to figure out if there is any way we can take a square root mod m when the discriminant isn't square in Z. Because that would still be the most ideal option. If I cant get that to work I'll go back to exploring my initial idea here.

Update: Actually I have an idea for tomorrow. Calculating the divisor is very hard, it would be an area of research which would take time. However, if we can make it so our divisor value is very small, which ends up being the case if the discriminant is square.. then it does not matter. So tomorrow I should have a look at cases, for example when the discriminant is 3\*3\*2  .. aka.. nearly square. What happens with this divisor value if it is almost square? Because if I can calculate those cases.. then that would be an historical breakthrough and a major improvement on existing algorithms. Lets see tomorrow.

Update: I spent the day working on this. I think the only way to do it is as following: If the discriminant for example is 7\*7\*7 then we can calculate the singular root for that. THen we need to figure out the divisor. But, and here is the crucial part.. rather then looking at solutions in a different prime field.. we look at solutions in 7\*7\*7\*7 instead. Because that then simplifies all this residue math. I dont know if I can see another way to do it. I guess if this fails I'll go back and try to optimize what I had before. JUst quickly finding large concentrations of small factors isn't a bad strategy either.. I just havnt yet tried hard enough to implement an optimized PoC for it. Although I am fairly sure if the discriminant is nearly square.. I should be able to do something with it.. it's definitely an angle I'll keep exploring for a while.

Update: You know what. Tomorrow I'm just going to go with the original idea I had here. Just a quadratic sieve variant trying to get many small factors in the discriminant so we can work with small factor bases. Get a super optimized PoC for that published. I can explore all these other things later. Its time I produce something I can show people instead of chasing perfection. 

#### To run from folder "Coefficient_Sieve" (Experimental WORK IN PROGRESS):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 30 -base 40 -debug 1 -lin_size 5000 -quad_size 100 -d 2

This is kind of NFS with second degree polynomials.. this needs to be expanded to higher degrees eventually..

Update: I quickly changed something so it takes a square root over a single large prime. Next we need to explore to different routs:

1. Can we succeed when the discriminant is not square? And if so, how and when?
2. Can we gain an advantage using higher degree polynomials?

#### To run from folder "CUDA_QS_variant":</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run: To run:  python3 run_qs.py -keysize 240 -base 100_000 -debug 1 -lin_size 100_000_000 -quad_size 100</br></br>
 
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

This was an attempt at finding smooths with similar factorization using an SIQS variant. By using quadratic coefficients. But it didnt end up working as I had hoped so I abondoned this approach, but perhaps someone will get some use out of it.

#### To run debug.py" (Prints the linear and quadratic coefficients to solve for 0 in the integers, for use with my paper):</br></br>

To run: python3 debug.py -keysize 12

This basically creates a system of quadratics. Solving them mod p is easy. But there is only one root solution (the factor of N) which solves the system for 0 for any mod p (aka solves it in the integers). Figuring out how to exactly do this quickly is still an ongoing area of research for me. And if a polynomial time algorithm for factorization exists, it is likely done by solving this system of quadratics. Finding a polyomial time algorithm is my ultimate goal, as this would make progress toward solving p = np as well. 


