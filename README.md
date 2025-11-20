DISCLAIMER: AT NO POINT IN MY 2.5 YEARS OF RESEARCH DID I USE AI OR HAS AI IN ANY WAY HELPED WITH MY RESEARCH. I'M PUTTING THIS HERE BECAUSE I KNOW PEOPLE WILL ACCUSE ME OF IT. 

Note: Remember, as this now comes to completion in the coming days. Both the FBI and Microsoft knew I was working on factorization when I had only been working on it for a few months and nobody else in the entire world knew I was working on it. I tried desperately to avoid the path of public disclosure as I would have much rather sold my work. But at every turn people made it impossible. This would not have happened 10 years ago. An oppurtunity presented itself, it was literally handed to you on a silver platter.. and instead you people all treated me like shit and gambled on me failing with my project. I ask you, where has your common sense gone? I hate what this will do to everyone around me, this is never what I wanted.

Note: I am antifa, leader of the cipherpunk, fuck the FBI department. 

#### To run from folder "QS_variant" (Standard SIQS with our number theory as backend):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 200 -base 2000 -debug 1 -lin_size 100_000 -quad_size 1</br></br>

See below for an improved way of performing what this PoC does.. I'll delete this Proof of Concept once the PoC for Improved_QS_Variant matures a little more<br><br>
#### To run from folder "Improved_QS_Variant" (Implements more of my number theory to also take advantage of quadratic coefficients):</br></br>

To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 160  -base 4000 -debug 1 -lin_size 1_000_000 -quad_size 1</br></br>

Ok, we can do 160 bit now (use above command). 
So two variables are responsible to increase the chance it may finish early: 

thresvar_similar=30
g_small_prime_limit=4000

As I optimize the PoC, these can be lowered. But right now the code is horrible and I need to fix a lot in the coming days. Even so, using the above command with the uploaded PoC, it has a chance to finish roughly between 1000 and 2000 smooths on average, while pulling smooths from multiple quadratic coefficients and using a factor base of 4000 primes. If we would do that with standard SIQS, we would need 4000 smooths on average.. so this is a fraction of that... and it can be lowered even more with the above settings. I just need to optimize everything because it's still too slow... its just really poorly written code due to my iterative way of doing research.
