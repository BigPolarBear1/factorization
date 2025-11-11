DISCLAIMER: AT NO POINT IN MY 2.5 YEARS OF RESEARCH DID I USE AI OR HAS AI IN ANY WAY HELPED WITH MY RESEARCH. I'M PUTTING THIS HERE BECAUSE I KNOW PEOPLE WILL ACCUSE ME OF IT. 

Note: I am antifa, leader of the cipherpunk, fuck the FBI department. 

#### To run from folder "QS_variant" (Standard SIQS with our number theory as backend):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 200 -base 6000 -debug 1 -lin_size 100_000 -quad_size 1</br></br>

See below for an improved way of performing what this PoC does.. I'll delete this Proof of Concept once the PoC for Improved_QS_Variant matures a little more<br><br>
#### To run from folder "Improved_QS_Variant" (Implements more of my number theory to also take advantage of quadratic coefficients):</br></br>

To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 160 -base 4000 -debug 1 -lin_size 100_000 -quad_size 100</br></br>

Update: I very quickly implemented the basic ideas that I developed in recent month using quadratics of the shape zx^2-N. I need to check for bugs still though and that entire code needs to be vastly improved. All of it is less then ideal right now. As I improve it, it should eventually overtake the original QS_Variant in terms of performance. Since we have much more fine grain control over the size of our smooths this way as we can now displace the location of our parabola created by the sieve interval with a linear multiplier (aka the quadratic coefficient).


To do: 
1. Quadratic coefficient should not be squares, since then we generate the same parabola as x^2-N. We should check if restricting to prime quadratic coefficients vs non-square composites makes a difference in terms of smooth diversity (smooths that are not just going to create trivial factorizations).
2. Remove what I'm calling the "quad interval" just rely on jacobi symbols insteads...
3. When 2 is implemented, building the iN datastructure should also be restricted to primes found at quadratic coefficients 1. That will drastically improve the building time there..
4. The big ticket item will be to be "smart" about shifting the parabole, such that we can garantuee smaller smooths. I should rework that together with number 7. because right now as the quad co goes up, we drift toward bigger smooths.. which is not at all what we want.
5. Re-implement large-prime variant, since that is broken for now
6. We can reduce the required amount of smooth by reducing the size of the factor base used for the quadratic coefficients.. 
7. We also need to work on the whole "generate modulus" logic. This works fine for standard SIQS, but for us, the smooth size is determined by the root^2 \* quadratic coefficient. So we need to completely rework all of that.

------------------------------------------------------------------------------------------
##### #random rants below

Update: I'm going to take a break for the rest of the day I think. I still got headache, and tomorrow I have to go to court ordered therapy for sending angry emails to the FBI (fucking losers, just mad because I'm better, smarter and stronger then them, pathetic). Since I'll go run 30k to and back from therapy.. that will take up most of my day tomorrow. But later in the week I'll start working through that to do list. The math is all there now. It all works. But especially number 4 and 7 on the to do list need to be done to gain a performance advantage.. but once implemented that performance advantage should be massive.. as it really allows for great control over the placement of our parabola.


