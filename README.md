DISCLAIMER: AT NO POINT IN MY 2.5 YEARS OF RESEARCH DID I USE AI OR HAS AI IN ANY WAY HELPED WITH MY RESEARCH. I'M PUTTING THIS HERE BECAUSE I KNOW PEOPLE WILL ACCUSE ME OF IT. 

Note: Remember, as this now comes to completion in the coming days. Both the FBI and Microsoft knew I was working on factorization when I had only been working on it for a few months and nobody else in the entire world knew I was working on it. I tried desperately to avoid the path of public disclosure as I would have much rather sold my work. But at every turn people made it impossible. This would not have happened 10 years ago. An oppurtunity presented itself, it was literally handed to you on a silver platter.. and instead you people all treated me like shit and gambled on me failing with my project. I ask you, where has your common sense gone? I hate what this will do to everyone around me, this is never what I wanted.

Note: I am antifa, leader of the cipherpunk, fuck the FBI department. 

#### To run from folder "QS_variant" (Standard SIQS with our number theory as backend):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 200 -base 2000 -debug 1 -lin_size 100_000 -quad_size 1</br></br>

See below for an improved way of performing what this PoC does.. I'll delete this Proof of Concept once the PoC for Improved_QS_Variant matures a little more<br><br>
#### To run from folder "Improved_QS_Variant" (Implements more of my number theory to also take advantage of quadratic coefficients):</br></br>

To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 100  -base 2000 -debug 1 -lin_size 100_000 -quad_size 1</br></br>

Update: I GOT IT!!!

So in the main sieve interval, if we only mark with even exponents and we use square moduli ... then we can generate smooths with incredibly few factors. Then when we use find_similar() this will increase the odds that the linear algebra step can succeed sooner. And the PoC demonstrates this already. I need to optimize everything now. Finally! I did it!

I'll go for a run now and call it a day. it works! It fucking works! Big day tomorrow as we begin optimizing all the code :). 

Update: Actually let me swap tomorrow's rest day with today and run tomorrow instead. Fucking got the flu or something. My entire body is feeling inflammed. I'll just eat a lot of protein today and I'm sure tomorrow I'll be recovered. Hahahaha. I am invincible. I refuse to accept the limits of this mortal body. Anyway, I'll still take a break from work. It works. I had my PoC finish way earlier with way less smooths, at 100 bits (back in september I could only pull this off with very small numbers) .. so it works. It proves my ideas. It proves that I was right. So today I feast and do fuck all for the rest of the day. I think I deserved it.
