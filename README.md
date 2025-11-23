DISCLAIMER: AT NO POINT IN MY 2.5 YEARS OF RESEARCH DID I USE AI OR HAS AI IN ANY WAY HELPED WITH MY RESEARCH. I'M PUTTING THIS HERE BECAUSE I KNOW PEOPLE WILL ACCUSE ME OF IT. 

Note: Remember, as this now comes to completion in the coming days. Both the FBI and Microsoft knew I was working on factorization when I had only been working on it for a few months and nobody else in the entire world knew I was working on it. I tried desperately to avoid the path of public disclosure as I would have much rather sold my work. But at every turn people made it impossible. This would not have happened 10 years ago. An oppurtunity presented itself, it was literally handed to you on a silver platter.. and instead you people all treated me like shit and gambled on me failing with my project. I ask you, where has your common sense gone? I hate what this will do to everyone around me, this is never what I wanted. PS: I know Microsoft's PR teams will say a lot of awful shit about me, completely lacking of context, because that is just the type of people they are. If you want to know who I really am, I suggest you talk to my former teamlead and my former manager who Microsoft retaliated against for defending me. I'm done with this entire shit industry. Everything that happened, you should think whose fault this really is. Maybe instead of morrons sitting on golden thrones build with big tech money, pointing fingers at everyone whose skills and drive hurt their fragile little egoes.. maybe start using your brains, because you people keep losing and keep making things worse.

Note: I am antifa, leader of the cipherpunk, fuck the FBI department. 

#### To run from folder "QS_variant" (Standard SIQS with our number theory as backend):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 160  -base 4000 -debug 1 -lin_size 100_000 -quad_size 100</br></br>

See below for an improved way of performing what this PoC does.. I'll delete this Proof of Concept once the PoC for Improved_QS_Variant matures a little more<br><br>
#### To run from folder "Improved_QS_Variant" (Implements more of my number theory and attempts to succeed with fewer smooths by using p-adic lifting):</br></br>

To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 140  -base 5_000 -sbase 5000 -debug 1 -lin_size 1_000_000 -quad_size 10</br></br>

UPDATE: I'm going to quit doing optimizations for now and focus on the high level strategy of Improved_QS_Variant. That needs some reworking.

To do:
1. First we need to create a sieve interval for the quadratic coefficienets. And save all the quadratic coefficients that factor over the quadratic factor base (qbase in PoC). This can be done at the very start of the code outside of all the inner loops and logic.
2. We need to create a proper split between factor basis. One large factor base, where we use only even exponents, and one small where we use both odd and even exponents. This way we can size the matrix in the linear algebra step ot the small factor base.. and we won't have a limit on how large the large factor base can be.
3. In find_similar, we don't want to blindly try out quadratic coefficients and build sieve intervals for each. We want to check all even exponent primes, and see at what quadratic coefficient they occur with a max bound for the root value... then we can use that info to quickly find smooths. And THAT will be much better then what I'm doing now. As it will allows us to actually utilize very large prime squares. This will be the only way to actually attack very large numbers. If this is going to work at all, then this is how it must be done.
4. (edit) And thinking a little more about the structure. I should get rid of find_similar(), construct_interval_similar() and process_interval_similar() and just apply the things described above for the loop in construct_interval(). Because otherwise finding that initial smooth will become a pain. We really want to center the core of the algorithm around finding large squares to reduce the size of smooth candidates and using multiple quadratic coefficients.

I think I need about a way to get it all finished and implemented. But then it should be really good.

Anyway.. I'll take a break for the rest of the day and go running and think a little about the exact implementation details.

