Disclaimer: At no point did AI contribute anything to this research project. Copilot can't even calculate basic congruences without making mistakes. People who think AI can do novel math research are delusional. I did not use AI for the paper, and I did not use AI for any of the code. I attacked this problem, from the very basics 3 years ago, slowly burrowing deeper and deeper. </br>
Disclaimer2: I never finished highschool. I started working on this 3 years ago, not even knowing basic algebra (I can explain exactly how I can to all the conclusions I came to if someone ever cares to hear about it). This is 3 years. I'm only going to get better. I'm still far away from plateauing skill wise.

The day I break factorization will be day 0 of the gay future. A better world will be born. A world of human creativity. A world of art, science and spirituality. A world without suffering. I will see that world become reality, no matter the cost. 

#### References (re-used many of the core number theoretical functions from these PoCs to fit my own algorithm): 
https://stackoverflow.com/questions/79330304/optimizing-sieving-code-in-the-self-initializing-quadratic-sieve-for-pypy
https://github.com/basilegithub/General-number-field-sieve-Python 

#### About the paper
Math paper is a work in progress. The final chapters are a bit rushed and building an algorithm around p-adic lifting isnt as straight forward as I had assumed. I do think there is an angle there I can exploit, but I'll do some further experimentation first and get a working PoC before I make edits to the paper again.

#### To run from folder "binomial_sieve":</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run:  python3 run_qs.py -keysize 35 -base 60 -debug 1 -lin_size 1000 -quad_size 100

Sieves with binomial expansions... 

Update: I have changed some things so it can sieve with non-zero derivatives.. but sieving degrees > 2 needs to be refined and in addition we need to implement residue sieving. 

Update: What I'll do next is simple... do an approach similar to CUDA_QS_Variant and if we find a smooth candidate using second degree polynomials... just expand to the 4th degree and sieve.. and so on. That should hopefully allow us to find smooth candidates with similar factorization... and succeed at the linear algebra step with much fewer smooth candidates.. bc traditional QS variants they cant punch past that 110 digit ceiling due to the fact of how quickly the required factor base size grows.. but this approach would solve that. And if that works... I'll pubish and take some time to seriously study math to begin the journey towards a polytime algorithm. That is my dream. Nobody can take my dream away. No matter how hard they try.

Update: I'll add a sieve interval for the second degree tomorrow. Then once a smooth candidates is found we can go hunting for similar factorizations on both sides of the congruence using binomial expansions.. thats the benefit of that setup..because we can ensure both sides of the congruence have similar factorizations... hopefully I can finish a function that goes hunting for similar factorization at higher degrees towards the weekend. If it works.. that concludes my work and I would have a prove-able breakthrough.. it would be good to see the end. After that I should go to a doctor and have these gut issues checked out.. does not seem to go away.. at this point its either a bad case of SIBO or colon cancer lol, just persistent localized pain for well over 6 months now.. its annoying bc I cant think straight anymore with the constant pain. If its cancer, the so be it, I'll go into the woods and take my life. Tired of this bullshit anyway. 

Update: Quickly added a sieve interval. So next we will have our mainloop exclusively using second degree polynomials.. once we hit a smooth candidate we then use binomial expansions to higher (even) degrees.. since they come from the same binomial term and as root we use multiples of this binomial term, it means atleast one side of the congruence will have similar factorizations.. for the other side we need to use precalculated residue maps.. but unlike normal quadratic sieve.. since we can sieve with polynomials having non-zero derivatives.. we have much better control over this. Anyway... I'll release this toward the end of the week. In addition.. we could also use moduli... like in my SIQS variant.. which doubles the speed, but since that complicates the number theory a fair bit, we'll implement that once everything else is done.

Update: Refactored some things. So we start sieving with linear polynomials... collects all the factors for b-smooths founds.. and then does binomial expansions to try and find similar factorizations. Fairly straightforward set-up ... but I now need to fix a lot of shit and do a lot of improvements. Idea of this is that we can make the amount of B-smooths required independent from the factorbase size... 

UPDATE: EUREKA! Sieving those multiples of N, I see it. Its often hitting very similar factorizations, and not just with factors I'm sieving for.. I see whats going on now.. let me improve the PoC so we zero in on this behavior..

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


