Disclaimer: At no point did AI contribute anything to this research project. Copilot can't even calculate basic congruences without making mistakes. People who think AI can do novel math research are delusional. I did not use AI for the paper, and I did not use AI for any of the code. I attacked this problem, from the very basics 3 years ago, slowly burrowing deeper and deeper. I did experiment from time to time with AI, but it always failed to be useful at math. It is useful for writing simple and short python functions for simple things.. but that is the extend of its usefulness. Saying this here because I am closing in fast now and I know people will claim I just used AI or something bc of my non-math background and being transgender. And I'm honestly tired of the trolls feeding my work into AI and then attacking my work using the garbarge analysis AI makes of it.</br>

Disclaimer2: I never finished highschool. I started working on this 3 years ago, not even knowing basic algebra (I can explain exactly how I can to all the conclusions I came to if someone ever cares to hear about it). This is 3 years. I'm only going to get better. I'm still far away from plateauing skill wise. 

The day I break factorization will be day 0 of the gay future. A better world will be born. A world of human creativity. A world of art, science and spirituality. A world without suffering. I will see that world become reality, no matter the cost. 

Been unemployed for years without income, still looking for work: big_polar_bear1@proton.me ... able to relocate to anywhere in the world except the US.

#### References (re-used many of the core number theoretical functions from these PoCs to fit my own algorithm): 
https://stackoverflow.com/questions/79330304/optimizing-sieving-code-in-the-self-initializing-quadratic-sieve-for-pypy
https://github.com/basilegithub/General-number-field-sieve-Python 

#### About the paper
Math paper is a work in progress. The final chapters are a bit rushed and building an algorithm around p-adic lifting isnt as straight forward as I had assumed. I do think there is an angle there I can exploit, but I'll do some further experimentation first and get a working PoC before I make edits to the paper again.

#### To run from folder "polysieve" WORK IN PROGRES...extremely early version:</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run:    python3 run_qs.py -keysize 100 -base 500 -debug 0 -lin_size 10_000 -quad_size 1

You need to use the set-up from find_same and just calculate residues and use p-adic lifting. There has to exist some multiple of N (part of which can also be encoded by the leading coefficient) where f(x) will be reducible in Z... meaning the discriminant will be square. In find_same those coefficients and roots are bounded by co_ind .. so you would only consider residues within that bound.. as not to have to deal with residue explosion because of CRT or lifting. So we can just take a few primes, keep lifting them.. and if a solution exists in Z.. within these bounds, then it will show up.

Update: Alright.. .this PoC is OKAY-ish for now. What I actually want to do in find_same() is closer to Coefficient_Sieve.. if a bounded set of coefficients produces a discriminant that is a square residue in enough fields (Z/p) ... at some point the discriminant is garantueed to be square. Within the condition that the coefficients are bounded. So we basically want to do exactly this.. but we need to properly define the bounds for our coefficients. Let me work out that math while in Iceland.... 

Probably should do something with those leading coefficients of f(x) and i(x) ... I'll work it out in Iceland.

Update: You know what also hits me... the more factors of a B-smooth are square, the less ground an algo like coefficient_sieve (working with quadratics and quadratic residues) has to cover to make it completely square. It's like an SIQS type reduction.. but for coefficient sieve... and in theory if we can use linear algebra to produce a nearly square B-smooth... then that should become even easier to make square. I'm not going to be able to finish this math before Iceland.. I will do it while there.

I am still looking for employment: big_polar_bear1@proton.me , belgian citizenship, cant travel to the US anymore. Maybe after the next administration with a lot of convincing (which has to include payment to my former manager bc of what microsoft did. He was the only person in this shit industry to have ever believed in me and he lost everything bc of it). I'm also preparing a course on rings and modulus to work through while in Iceland so I can master the number field approach of NFS. Once I master that material.. I'll be right there at the cutting edge.. and if I dont have all the chess pieces already.. I will have them soon. Either I die or I succeed, because I'm sure as hell not going to abandon this project for as long as I am alive and every year I get better and better and better. Couldn't do basic algebra 3 years ago. 

#### To run from folder "Coefficient_Sieve" (For use with the paper):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run:  python3 run_qs.py -keysize 40 -base 50 -debug 1 -lin_size 10_000 -quad_size 100</br>

Just demonstrates the math from the paper using quadratics. For educational purposes. And rather then taking a square root over a large prime we can also just calculate the discriminant. But this demonstrates the interesting relation between these quadratics and the factors of N.

#### To run from folder "CUDA_QS_variant" (Failed Experiment):</br></br>
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


