Disclaimer: At no point did AI contribute anything to this research project. Copilot can't even calculate basic congruences without making mistakes. People who think AI can do novel math research are delusional. I did not use AI for the paper, and I did not use AI for any of the code. I attacked this problem, from the very basics 3 years ago, slowly burrowing deeper and deeper. </br>
Disclaimer2: I never finished highschool. I started working on this 3 years ago, not even knowing basic algebra (I can explain exactly how I can to all the conclusions I came to if someone ever cares to hear about it). This is 3 years. I'm only going to get better. I'm still far away from plateauing skill wise.

The day I break factorization will be day 0 of the gay future. A better world will be born. A world of human creativity. A world of art, science and spirituality. A world without suffering. I will see that world become reality, no matter the cost. 

Been unemployed for years without income, still looking for work: big_polar_bear1@proton.me ... able to relocate to anywhere in the world except the US. Depression man. Just really intense suicidal thoughts today. Everything that was good in my life has been gone for years now. Just been living a life without dignity, stuck in a tiny room, no income for years. Perhaps when it is an exceptionally nice day, I'll just do it, get it over with. Fck this bullshit. What a world filled with morrons.

#### References (re-used many of the core number theoretical functions from these PoCs to fit my own algorithm): 
https://stackoverflow.com/questions/79330304/optimizing-sieving-code-in-the-self-initializing-quadratic-sieve-for-pypy
https://github.com/basilegithub/General-number-field-sieve-Python 

#### About the paper
Math paper is a work in progress. The final chapters are a bit rushed and building an algorithm around p-adic lifting isnt as straight forward as I had assumed. I do think there is an angle there I can exploit, but I'll do some further experimentation first and get a working PoC before I make edits to the paper again.

#### To run from folder "binomial_sieve":</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run:  python3 run_qs.py -keysize 75 -base 500 -debug 1 -lin_size 1_000 -quad_size 1_000_000 -t 70

Update: Already, added it so that it creates a 2d interval. Now only the important part is still remaing. See To Do below.

To Do:

If you run the PoC, you will occasionally notice the PoC finding B-smooths that are a slight variant of each other (only 1 or 2 small factors are different) ... this happens a lot actually and these count as unique smooth candidates. The reason this is useful is because they cancel out eachothers large factors and give us just a handful of small factors for the linear algebra step. And if we can find a lot of these... we can potentially succeed at the linear algebra step much sooner.
Now, the current problem is, that these might fall outside our sieve interval.. but we can quickly check if they exist using our residue map. So thats what needs to be done now. 

A good example is for example this: 

Poly: [1, 4, -33032964601443606179686070] x: 4670 pval: -33032964601443606157858490 seen_primes: [-1, 2, 5, 17, 19, 89, 173, 317, 419, 467, 1949, 1949, 2819] k: 88730</br>
Poly: [1, 4, -33089552121103466747017438] x: 4674 pval: -33089552121103466725152466 seen_primes: [-1, 2, 17, 19, 89, 173, 317, 419, 1949, 1949, 2339, 2819] k: 88882</br>

These are near identical but have a handful of different factors. You will find cases like this all the time. And this behavior is quite important.

x1: 4670 = 2 x 5 x 467    </br>
x2: 4674 = 2 x 3 x 19 x 41 </br>
k1: 88730 = 2 x 5 x 19 x 467</br>
k2: 88882 = 2 x 19 x 2339</br>

Now if we look at the factorizations of the root and multiples of N, k ... we see that the factors that are different show up in these. So that gives us a clue how to help checking for these variations. These are unique because the non-unique ones are when you multiply a root by a factor and multiply k with the square of it. Which ends up adding a square to the polynomial value. These are non square factors being added. Hence these are completely valid for us to us during linear algebra. My work is nearing completion now.... perhaps finally the suffering will end... or perhaps it is only beginning.

Update: Eh, starting to get there.. but need to optimize everything now. And there will probably also be a benefit to adding p-adic lifting when we go into find_same2.

Update: After doing some testing and thinking today. I dont believe sieving is the correct approach for find_same2 ... just got to be more deliberate trying to find these near identical smooths by using the fact that those different factors show up in the k variable (the multiples of N). Anyway... it is the final thing that needs to be figured out somehow now.

Update: I guess if we have i.e x^2+b\*x-N\*k then multiplying root x and multiplying k basically shifts the entire thing. And which "shifts" result in near identical smooth candidates..thats just residue math..I should be able to figure that out. I'll go for a run.. let me give it a go when I get back. Either this will work or fck it.. I'm just going to end my life. Fcking tired of shit.

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


