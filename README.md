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
To run: python3 run_qs.py -keysize 80 -base 10_000 -debug 0 -lin_size 100 -quad_size 100

This PoC now requires as number of B-smooths to succeed around 1/4th the size of the factor base. Beyond a doubt proving my ideas. It can be reduced much more though. 
It achieves this by finding an initial B-smooth using square moduli, such that half the bitlength can be discarded as the factors of the square modulus are irrelevant as far as guassian elimination over Z/2 is concerned, and then trying to find more B-smooths that have factor overlap with the initial B-smooth.. but by using much smaller moduli and using the coefficients to generate smaller B-smooth candidates.. far fewer new factors are introduced.. because of how these factors end up canceling eachother out (see final chapter paper for an example) we can succeed much easier. There is no trickery involved here. It works like I had theorized it might since september already.

To do:

1. Need to presieve a bunch of numbers roughly within range of what these quadratic coefficients will end up being, so that we can quickly check if there is a leading coefficient nearby that factors over the factor base.
2. The paper demostrated using third degree polynomials.. this has the benefit of having more roots.. but I'm still unsure if it yield a performance increase. So let me just go ahead and implement the POC using purely quadratics and large quadratic coefficients. If I come across a scenario where having more roots per prime ends up as an advantage then I'll implement that later.

I'll add some pre-sieving logic next.

Update: Quickly added some bug fixes. Will go ahead and slowly start implementing pre-sieving values at a certain bit lenght now. Let me just start by pre-calculating residues and creating an interval, and then implement something more elaborate.

Update: Quickly added a sieve interval to make the PoC a little less crappy. If you run the above command, using a factor base of 10_000 primes (that is 10_000 primes, not primes up to 10_000.. important distinction), it will succeed at around 3000 B-smooths. Which is less then half.
Now the thing is... we can precompute numbers that factor over a small factor base within a certain range.. and then only check roots of a certain bitlength that we know will end up inside or near that range. That's how you factorize! People would have known I was right. Yet my life has been reduced to literally living like a fcking cockroach. I'm not sure what slander and lies they will weaponize against me when this finally gets out... I guess they are hoping by that time I'm dead and unable to speak for the truth.

Update: Shit day today. This heatwave is killing me. And I'm living in this tiny room that somehow seems to trap all the heat. Can hardly sleep. Let me try and implement some rudimentary pre-sieving tomorrow. And since this would be independent of N, these results can also be saved to disk and re-used. First I'll add the modulus to the root, so the root is always a specific bitsize.. whose optimal quadratic coefficient to generate the smallest possible polynomial value is always going to be inside or near the pre-sieved range. So I'll implement that logic first.. and then add some pre-sieving logic... easy enough.

Update: Made some progress toward pre-sieving.. will be ready to upload tomorrow. Need to change a few more small things and it will be done. Really wasnt too difficult.

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


