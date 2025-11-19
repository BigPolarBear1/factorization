DISCLAIMER: AT NO POINT IN MY 2.5 YEARS OF RESEARCH DID I USE AI OR HAS AI IN ANY WAY HELPED WITH MY RESEARCH. I'M PUTTING THIS HERE BECAUSE I KNOW PEOPLE WILL ACCUSE ME OF IT. 

Note: Remember, as this now comes to completion in the coming days. Both the FBI and Microsoft knew I was working on factorization when I had only been working on it for a few months and nobody else in the entire world knew I was working on it. I tried desperately to avoid the path of public disclosure as I would have much rather sold my work. But at every turn people made it impossible. This would not have happened 10 years ago. An oppurtunity presented itself, it was literally handed to you on a silver platter.. and instead you people all treated me like shit and gambled on me failing with my project. I ask you, where has your common sense gone? I hate what this will do to everyone around me, this is never what I wanted.

Note: I am antifa, leader of the cipherpunk, fuck the FBI department. 

#### To run from folder "QS_variant" (Standard SIQS with our number theory as backend):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 200 -base 2000 -debug 1 -lin_size 100_000 -quad_size 1</br></br>

See below for an improved way of performing what this PoC does.. I'll delete this Proof of Concept once the PoC for Improved_QS_Variant matures a little more<br><br>
#### To run from folder "Improved_QS_Variant" (Implements more of my number theory to also take advantage of quadratic coefficients):</br></br>

To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 80  -base 2000 -debug 1 -lin_size 100_000 -quad_size 1</br></br>

Update: I GOT IT!!!

Alright. It's working. Just run it on 80 bit keys using "python3 run_qs.py -keysize 80  -base 2000 -debug 1 -lin_size 100_000 -quad_size 1" for now. There is some issues I still need to fix. But PoC demonstrates its ability to finish with much fewer smooths. In theory we dont even need to call into find_similar. Just generating square moduli and marking the sieve interval with even exponents only works to generate smooths with very few odd exponent factors... hence increasing the odds of succeeding early. However, if we then also call into find similar.. then we find more smooths that include the large primes of the smooth we found.. which should help with very large numbers and large factor bases. So in a way, it also serves a role. Tomorrow I need to also iterate quadratic coefficients in the main loop (before calling into find_similar) ... instead of just coefficient 1, because an issue we have now is that eventually it ends up failing to generate a modulus, because using square moduli, it's harder to get to a target bit size.. but simply switching to the next quadratic coefficient if generating the modulus fails works too. A lot of optimization works now needs to be done.

BUT IT WORKS! Back in September I could only get this to work with very small numbers. Now that it works with big numbers aswell... I finally got something that's breaking new ground, I just need to get it faster now. And even if it is a little slower in the end, as long as it has a chance of succeeding way sooner... that's the important thing. Just having that chance of completing much sooner. Especially an advantage when dealing with very large numbers and very large factor basis.

PS: One last rant before I try to sleep. It's funny, you people are continueing this deception, treating me like shit, wasting literal years of my life. I would have sold my work. I tried to. You all knew I was working on this all the way back in 2023. Nothing was done. Instead I come back to Europe, and right from the start, I get treated just the same here like I was treated around 2018. Like nothing had ever changed. Morons. Everyone working in western intelligence should quit their fucking jobs. You're all pathetic. Literally, quit your jobs you morons. You had all the warnings in the world, yet you chose the one course of action that resulted in this. You had an oppurtunity presented to you on a silver platter. Yet you choose this course of action. I'll bash your fucking heads in if our paths should ever cross. Fucking dimwitted morrons. They dont hire the brightest people in that intelligence business of yours huh? Losers. Quit your jobs. You all suck.
