DISCLAIMER: AT NO POINT IN MY 2.5 YEARS OF RESEARCH DID I USE AI OR HAS AI IN ANY WAY HELPED WITH MY RESEARCH. I'M PUTTING THIS HERE BECAUSE I KNOW PEOPLE WILL ACCUSE ME OF IT. 

Note: Remember, as this now comes to completion in the coming days. Both the FBI and Microsoft knew I was working on factorization when I had only been working on it for a few months and nobody else in the entire world knew I was working on it. I tried desperately to avoid the path of public disclosure as I would have much rather sold my work. But at every turn people made it impossible. This would not have happened 10 years ago. An oppurtunity presented itself, it was literally handed to you on a silver platter.. and instead you people all treated me like shit and gambled on me failing with my project. I ask you, where has your common sense gone? I hate what this will do to everyone around me, this is never what I wanted.

Note: I am antifa, leader of the cipherpunk, fuck the FBI department. 

#### To run from folder "QS_variant" (Standard SIQS with our number theory as backend):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 200 -base 2000 -debug 1 -lin_size 100_000 -quad_size 1</br></br>

See below for an improved way of performing what this PoC does.. I'll delete this Proof of Concept once the PoC for Improved_QS_Variant matures a little more<br><br>
#### To run from folder "Improved_QS_Variant" (Implements more of my number theory to also take advantage of quadratic coefficients):</br></br>

To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 100 -base 1000 -debug 1 -lin_size 10_000 -quad_size 1</br></br>

Update: I have began implementing a strategy that would allow the algorithm to finish much sooner. See the find_similar() function. I need to improve it... just uploading my work in progress for the day before I tae a break now.

Strategy is a follows:

1. Use square moduli to find a smooth at quadratic coefficient  = 1
2. When we find a smooth, divide by the square modulus and using the factors with negative exponents, construct a new modulus. Because of using square moduli in step 1, the bit size of this one will be relatively low.
3. Now find smooths with these factors at different quadratic coefficients.

To do:

1. In step two, if the new modulus is too small, we can us p-adic lifting to increase the bit length
2. In find_similar, we should use a sieve interval
3. Once the sieve interval is implemented, we should mark the sieve interval with p^2 for primes that arn't in our modulus, so those end up being ignored in the linear algebra step... this will allow the linear algebra portion to potentially succeed much earlier.

I'll begin attacking that to do list tomorrow. But first, I need to clean up all of the code... it's becoming convoluted with a lot of unneeded code. Frankenstein code due to making too many edits.

The code I uploaded today is really bad. I need to change and add a lot of things in the coming days. But I know I'm finally on the right track and that my work's completion is imminent now. 

Update: Also quickly added the trick the generate permutations of linear coefficients in the find_similar() function. Next this needs a sieve interval.. that will be my goal for tomorrow. Then the algorithm should really start to ramp up....
