# factorization

The PoC under PoC_Files will easily factor above 200 bits. </br></br>

However, that one is basically just using the number theory from my paper to achieve default SIQS. </br></br>

in PoC_Files_Find_Similar_WIP I am attempting to sieve in the direction of the quadratic coefficient and linear coefficient. This way in theory we need much less smooths.</br></br>

For PoC_Files_Find_Similar_WIP:</br>
Build: python3 setup.py build_ext --inplace</br>
Run: python3 run_qs.py -keysize 30 -base 100 -debug 1 -lin_size 1_00 -quad_size 100_000</br></br>

You will see some output like this:</br></br>

Constructing quad interval took: 5.059455415990669</br>
Processing quad interval took: 0.16657884302549064</br>
Generating new modulus:  681</br>
Looking for similar smooths, amount needed: 6 initial smooth factors: {2, 227, 59, -1}</br>
local_factors:  [2, 59, 227]</br>
Looking for similar smooths, amount needed: 9 initial smooth factors: {2, 3, 227, 41, 79, 17, -1}</br>
local_factors:  [2, 3, 17, 41, 79, 227]</br>
Looking for similar smooths, amount needed: 8 initial smooth factors: {2, 227, 229, 7, 11, -1}</br>
local_factors:  [2, 7, 11, 227, 229]</br>
Looking for similar smooths, amount needed: 8 initial smooth factors: {3, 227, 7, 103, 43, -1}</br>
local_factors:  [3, 7, 43, 103, 227]</br>
Found a similar smooth {43, 2, 3, 227}</br>
Found a similar smooth {2, 3, 7, 43, -1}</br>
Found a similar smooth {2, 103, 7, -1}</br>
Found a similar smooth {3, -1}</br>
Found a similar smooth {3, -1}</br>
Found a similar smooth {3, 7}</br>
Found a similar smooth {227}</br>
sim_found:  7</br>
Running linear algebra step with #smooths:  11</br>
[SUCCESS]Factors are: 17203 and 18211</br>

If it finds a smooth, it will look in the direction of  the quadratic coefficient to find more similar smooths. And if it finds enough of them, it skips to the linear algebra step. You can succeed with much less smooths this way.
All the basic functionality is kind of there. But it won't work well yet.

I realized I need to address the following things:

1. Our very first step, righ after finding the initial smooth is performing p-adic lifting on the linear coefficients modulo the smooth factors (with odd exponents) using all possible linear coefficient solutions.
   Then we look at those solutions which lift very high while keeping a small coefficient (or dist1, dist2 values.. I need to think what to use as modulus.. but the idea is the lift them in advance so we can just find quadratic coefficients where they will cluster together, while minimizing the linear and quadrtic coefficient) and we find at which quadratic coefficients they occur. Ony then should we build the sieve interval (or possibly we could straight up calculate everything without sieve interval). The way it's doing it now is not going to work. Because we don't want to limit our selves to a small quad_sieve_size due to precalculating the primes that are quadratic residues for each. That's a wrong appraoch.

2. With that out of the way... it should be much better... then I need to fix a lot of other shit and shit code. Had some particulary bad brain fog when I wrote this find_similar code. I need to rework literally everything.  It's getting there. The basic building blocks to achieve it are there... it's just getting the implementation right now.

