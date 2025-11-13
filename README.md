DISCLAIMER: AT NO POINT IN MY 2.5 YEARS OF RESEARCH DID I USE AI OR HAS AI IN ANY WAY HELPED WITH MY RESEARCH. I'M PUTTING THIS HERE BECAUSE I KNOW PEOPLE WILL ACCUSE ME OF IT. 

Note: Remember, as this now comes to completion in the coming days. Both the FBI and Microsoft knew I was working on factorization when I had only been working on it for a few months and nobody else in the entire world knew I was working on it. I tried desperately to avoid the path of public disclosure as I would have much rather sold my work. But at every turn people made it impossible. This would not have happened 10 years ago. An oppurtunity presented itself, it was literally handed to you on a silver platter.. and instead you people all treated me like shit and gambled on me failing with my project. I ask you, where has your common sense gone? I hate what this will do to everyone around me, this is never what I wanted.

Note: I am antifa, leader of the cipherpunk, fuck the FBI department. 

#### To run from folder "QS_variant" (Standard SIQS with our number theory as backend):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 200 -base 6000 -debug 1 -lin_size 100_000 -quad_size 1</br></br>

See below for an improved way of performing what this PoC does.. I'll delete this Proof of Concept once the PoC for Improved_QS_Variant matures a little more<br><br>
#### To run from folder "Improved_QS_Variant" (Implements more of my number theory to also take advantage of quadratic coefficients):</br></br>

To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 160 -base 4000 -debug 1 -lin_size 100_000 -quad_size 100</br></br>

Update: I very quickly implemented the basic ideas that I developed in recent months using quadratics of the shape zx^2-N. I need to check for bugs still though and that entire code needs to be vastly improved. All of it is less then ideal right now. As I improve it, it should eventually overtake the original QS_Variant in terms of performance. Since we have much more fine grain control over the size of our smooths this way as we can now displace the location of our parabola created by the sieve interval with a linear multiplier (aka the quadratic coefficient).


##### #To do: 
1. Quadratic coefficient should not be squares, since then we generate the same parabola as x^2-N. We should check if restricting to prime quadratic coefficients vs non-square composites makes a difference in terms of smooth diversity (smooths that are not just going to create trivial factorizations).
2. Remove what I'm calling the "quad interval" just rely on jacobi symbols insteads...
3. When 2 is implemented, building the iN datastructure should also be restricted to primes found at quadratic coefficients 1. That will drastically improve the building time there..
4. The big ticket item will be to be "smart" about shifting the parabole, such that we can garantuee smaller smooths. I should rework that together with number 7. because right now as the quad co goes up, we drift toward bigger smooths.. which is not at all what we want.
5. Re-implement large-prime variant, since that is broken for now
6. We can reduce the required amount of smooth by reducing the size of the factor base used for the quadratic coefficients.. 
7. We also need to work on the whole "generate modulus" logic. This works fine for standard SIQS, but for us, the smooth size is determined by the root^2 \* quadratic coefficient. So we need to completely rework all of that.
8. Coefficient lifting to square moduli! But this will need to tie in with whatever we change in step 7.

------------------------------------------------------------------------------------------
##### #random rants below

UPDATE: I spent most of the day just taking this apart and pondering how to best appraoch it. I have settled on an approach. I'll upload tomorrow or in the weekend. We need to work with fairly small moduli, so we can use the modulus to shift the parabola by adding it to the quadratic coefficient. There will still be a "range" where we know we can get smaller smooth candidates by adding the modulus to the quadratic coefficient.. even if we can't factor all of the quadratic coefficients over the factor base.. as long as we can factor a few of them and have garantueed smaller smooths... that's a win. That's the way that should be done. The way the uploaded PoC does it is not good. It does create unique parabola, but using a single modulus, you need to shorten the sieve interval... otherwise your parabola runs from bitlen(N) to bitlen(N)+bitlen(quadratic coefficient). 

