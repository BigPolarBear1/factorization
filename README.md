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

Nearly done now. The whole find_similar where we try to find smooths with similar factorization as the original smooth... that now needs to be optimized so we can go over more quadratic coefficients fasters and so we can also lower the thresvar_similar variable.. which is the tolerance for how many factors not part of the original smooth to allow. Its a bit too lose in the uploaded PoC now, but I need to optimize that code first. One step towards optimizing that is to also mark p^2 in the original sieve intervals... since that will potentially allow for fewer negative exponent factors.. meaning our new modulus in find_similar will be smaller... which is important for creating a sieve interval there that creates a good parabola. I'll begin optimizations tomorrow. Really the goal is to lower that thresvar_similar variable down to atleast 30 (from 40) and be able to find smooths there.. that will then really increase the chance of succesful factorization way earlier with much fewer smooths.

note to self: I shouldn't restrict myself to p^2. I should keep lifting for as long as the resulting distance is within the sieve interval. But going to need a bit of refactoring first.. 

Update.. going to sleep. Tomorrow I may be able to achieve a performance advantage over default SIQS finally.. just a lot of rough edges to fix now. I dont know how I feel about. I dont think I have feelings left after the last 3 years.

The one thing Ill never understand is why my own country is taking part in this obvious deception while there is war in eastern europe. The world is ruled by morrons. I know people must have seen/noticed my work. I know what I got. The only way this can be explained is that people rathee throw away potential advantages of strategic importance then to deal with someone like me. The lunacy of it. Morrons. Im never talking to anyone anymore once this is finished. Want nothing to do anymore with people. Fucking despise you all 

Anyway, doesnt surprise me. Belgium has always been a dead end for me. Its an awful country. There really arnt great countries left in the world. Everything has gone to shit.
