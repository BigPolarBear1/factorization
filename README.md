DISCLAIMER: AT NO POINT IN MY 2.5 YEARS OF RESEARCH DID I USE AI OR HAS AI IN ANY WAY HELPED WITH MY RESEARCH. I'M PUTTING THIS HERE BECAUSE I KNOW PEOPLE WILL ACCUSE ME OF IT. 

Note: Remember, as this now comes to completion in the coming days. Both the FBI and Microsoft knew I was working on factorization when I had only been working on it for a few months and nobody else in the entire world knew I was working on it. I tried desperately to avoid the path of public disclosure as I would have much rather sold my work. But at every turn people made it impossible. This would not have happened 10 years ago. An oppurtunity presented itself, it was literally handed to you on a silver platter.. and instead you people all treated me like shit and gambled on me failing with my project. I ask you, where has your common sense gone? I hate what this will do to everyone around me, this is never what I wanted.

Note: I am antifa, leader of the cipherpunk, fuck the FBI department. 

#### To run from folder "QS_variant" (Standard SIQS with our number theory as backend):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 200 -base 2000 -debug 1 -lin_size 100_000 -quad_size 1</br></br>

See below for an improved way of performing what this PoC does.. I'll delete this Proof of Concept once the PoC for Improved_QS_Variant matures a little more<br><br>
#### To run from folder "Improved_QS_Variant" (Implements more of my number theory to also take advantage of quadratic coefficients):</br></br>

To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 100 -base 1000 -debug 1 -lin_size 100_000 -quad_size 1</br></br>

Update: I began improving everything today. From testing the following now needs to be done:

So its really important that our new modulus has a specific bit length to create just the size of the sieve interval we want in find_similar.
Hence some type of generate_modulus function like we use at the start of the algorithm needs to be implemented. It is important that in find_similar we find smooths that contain at the very least the large prime factors with odd exponents. We can use the small prime factors to adjust the bit length of the modulus.. since these small prime factors will be very common during the run of the algorithm anyway.
Once that is done... everything will be much better. Then next we also just need to do a lot more p-adic lifting.. not just to exponent 2. Because we really need to capitalize on those cases with high exponents clustered together on the sieve interval.
