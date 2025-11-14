DISCLAIMER: AT NO POINT IN MY 2.5 YEARS OF RESEARCH DID I USE AI OR HAS AI IN ANY WAY HELPED WITH MY RESEARCH. I'M PUTTING THIS HERE BECAUSE I KNOW PEOPLE WILL ACCUSE ME OF IT. 

Note: Remember, as this now comes to completion in the coming days. Both the FBI and Microsoft knew I was working on factorization when I had only been working on it for a few months and nobody else in the entire world knew I was working on it. I tried desperately to avoid the path of public disclosure as I would have much rather sold my work. But at every turn people made it impossible. This would not have happened 10 years ago. An oppurtunity presented itself, it was literally handed to you on a silver platter.. and instead you people all treated me like shit and gambled on me failing with my project. I ask you, where has your common sense gone? I hate what this will do to everyone around me, this is never what I wanted.

Note: I am antifa, leader of the cipherpunk, fuck the FBI department. 

#### To run from folder "QS_variant" (Standard SIQS with our number theory as backend):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 200 -base 2000 -debug 1 -lin_size 100_000 -quad_size 1</br></br>

See below for an improved way of performing what this PoC does.. I'll delete this Proof of Concept once the PoC for Improved_QS_Variant matures a little more<br><br>
#### To run from folder "Improved_QS_Variant" (Implements more of my number theory to also take advantage of quadratic coefficients):</br></br>

To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 160 -base 4000 -debug 1 -lin_size 1_000_000 -quad_size 1_000</br></br>

Update: I very quickly implemented the basic ideas that I developed in recent months using quadratics of the shape zx^2-N. I need to check for bugs still though and that entire code needs to be vastly improved. All of it is less then ideal right now. As I improve it, it should eventually overtake the original QS_Variant in terms of performance. Since we have much more fine grain control over the size of our smooths this way as we can now displace the location of our parabola created by the sieve interval with a linear multiplier (aka the quadratic coefficient).


##### #To do: 
1. Quadratic coefficient should not be squares, since then we generate the same parabola as x^2-N. We should check if restricting to prime quadratic coefficients vs non-square composites makes a difference in terms of smooth diversity (smooths that are not just going to create trivial factorizations).
2. Remove what I'm calling the "quad interval" just rely on jacobi symbols insteads...
3. When 2 is implemented, building the iN datastructure should also be restricted to primes found at quadratic coefficients 1. That will drastically improve the building time there..
4. The big ticket item will be to be "smart" about shifting the parabole, such that we can garantuee smaller smooths. I should rework that together with number 7. because right now as the quad co goes up, we drift toward bigger smooths.. which is not at all what we want.
5. Re-implement large-prime variant, since that is broken for now
6. ~~We can reduce the required amount of smooth by reducing the size of the factor base used for the quadratic coefficients..~~
7. We also need to work on the whole "generate modulus" logic. This works fine for standard SIQS, but for us, the smooth size is determined by the root^2 \* quadratic coefficient. So we need to completely rework all of that.
8. Coefficient lifting to square moduli! But this will need to tie in with whatever we change in step 7.

------------------------------------------------------------------------------------------
##### #random rants below

UPDATE: I spent most of the day just taking this apart and pondering how to best appraoch it. I have settled on an approach. I'll upload tomorrow or in the weekend. We need to work with fairly small moduli, so we can use the modulus to shift the parabola by adding it to the quadratic coefficient. There will still be a "range" where we know we can get smaller smooth candidates by adding the modulus to the quadratic coefficient.. even if we can't factor all of the quadratic coefficients over the factor base.. as long as we can factor a few of them and have garantueed smaller smooths... that's a win. That's the way that should be done. The way the uploaded PoC does it is not good. It does create unique parabola, but using a single modulus, you need to shorten the sieve interval... otherwise your parabola runs from bitlen(N) to bitlen(N)+bitlen(quadratic coefficient). 

UPDATE: Just uploading my work for the day. Didn't do much today aside from making the factor base for the quadratic coefficient seperate. To reduce the amount of required smooths there. Tomorrow I need to change some stuff... use much smaller moduli (maybe not run the parabola from N to N, since thats quite a big range to cover) ... then we can use the modulus to move the parabola around by adding it to the quadratic coefficient.

UPDATE: You know what I'm going to do first thing tomorrow? Graph that parabola. So that when I start messing with the quadratic coefficient, I can see exactly on a graph how it changes. That's going to save me some headache from staring at numbers all day.

UPDATE: I just quickly worked out the math of how adding the modulus to the quadratic coefficient displaces the parabola. So that gives us exact math .... we could just have a smooth candidate and basically shift it in one direction. I.e if we have a large negative smooth, we can then add the modulus to the quadratic coefficient to make it a smaller (closer to 0) number. That's actually super straight forward. Plus I guess we can use negative quadratic coefficients to make larger positive numbers smaller. That's actually super neat. It's just a lot of nitpicky coding now... I'll just take it easy.. my mental state is feeling pretty bad these last few days. It's just everything that has happened, and everything that will happen once my work is complete. It's starting to feel like a real nightmare. The world is ran by morrons who do not care about anything but themselves. Absolute fucking morrons.  

UPDATE: Anyway, I guess the way I'm going to try and finish the PoC is in the process_interval function, if we have a smooth, just make it smaller with the quadratic coefficient (add or subtract the modulus). As long as we're able to factor the quadratic coefficient, that still ends up being a valid smooth we can add to our linear algebra step. I'll make that my goal to finish before the end of the week. Then just take it from there. 

UPDATE: It's really finnicky to implement this. I think we just find smooths using quadratic coefficient = 1.. then for each smooth candidate, we can reduce it in size by using the quadratic coefficient to shift the smooth candidate. However, we should probably shift using a subset of the modulus and not the entire modulus since that will result in very large shifts. So the finnicky part is a fast heuristic to decide which shifts are worth doing. I guess I need an heuristic that will allow me to shift, create an as small as possible smooth, without sacrificing to much of the modulus (known factors). With is just really doing math with bit lengths. I'll need to generate some toy examples and do the math. This can take a few days.. but after that we should be nearly done.

UPDATE: I guess as modulus we should consider any factors of the original smooth candidate, not just the modulus it was generated with. I wonder if we can use that to salvage smooth candidates with a lot of factors that arn't fully smooth.  It's just coming up with a good heuristic now. It's doable. Especially now that we dont have to worry about multiples of N and can flatten it out like this, because that makes the math so much more complicated otherwise.
