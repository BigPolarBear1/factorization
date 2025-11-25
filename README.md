DISCLAIMER: AT NO POINT IN MY 2.5 YEARS OF RESEARCH DID I USE AI OR HAS AI IN ANY WAY HELPED WITH MY RESEARCH. I'M PUTTING THIS HERE BECAUSE I KNOW PEOPLE WILL ACCUSE ME OF IT. 

Note: Remember, as this now comes to completion in the coming days. Both the FBI and Microsoft knew I was working on factorization when I had only been working on it for a few months and nobody else in the entire world knew I was working on it. I tried desperately to avoid the path of public disclosure as I would have much rather sold my work. But at every turn people made it impossible. This would not have happened 10 years ago. An oppurtunity presented itself, it was literally handed to you on a silver platter.. and instead you people all treated me like shit and gambled on me failing with my project. I ask you, where has your common sense gone? I hate what this will do to everyone around me, this is never what I wanted. PS: I know Microsoft's PR teams will say a lot of awful shit about me, completely lacking of context, because that is just the type of people they are. If you want to know who I really am, I suggest you talk to my former teamlead and my former manager who Microsoft retaliated against for defending me. I'm done with this entire shit industry. Everything that happened, you should think whose fault this really is. Maybe instead of morrons sitting on golden thrones build with big tech money, pointing fingers at everyone whose skills and drive hurt their fragile little egoes.. maybe start using your brains, because you people keep losing and keep making things worse.

Note: I am antifa, leader of the cipherpunk, fuck the FBI department. 

#### To run from folder "QS_variant" (Standard SIQS with our number theory as backend):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 220 -base 20_000 -debug 1 -lin_size 10_000_000 -quad_size 1</br></br>

See below for an improved way of performing what this PoC does.. I'll delete this Proof of Concept once the PoC for Improved_QS_Variant matures a little more<br><br>

Note: With a large enough -base and lin_size this PoC will find smooths for 110 digits. Albeit very slowly, but this is a highly unoptimized cython PoC. However, to push beyond that into novel terroritory for Quadratic Sieve-based algorithms we need to use quadratic coefficients and p-adic lifting, and that is what the PoC below (Improved_QS_Variant) will be for. 

#### To run from folder "Improved_QS_Variant" (Implements more of my number theory and attempts to succeed with fewer smooths by using p-adic lifting):</br></br>

To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 140 -base 10_000 -sbase 5000 -debug 1 -lin_size 1_000_000 -quad_size 100</br></br>

UPDATE: Today I did some testing, to figure out the absolute best way to finish this PoC. 

What we need to do is this:

1. Generate a square modulus with generate_modulus()
2. Have a small factor base and a large factor base. The small factor base we use for trial factorization of smooth candidates and the large one we use to find large squares to reduce the bit size of smooth candidates.
3. This will be the most intensive loop of the algorithm and needs to be hyper optimized, but just keep calculating the root of the modulus*(prime^even_exp) for different quadratic coefficients up to a very large bound (as long as these quadratic coefficients factor over the quadratic factor base.. which is a third factor base outside of the small and large factor base from step 2). If the root becomes small enough we save it.... and repeat for all primes. And we will need to optimize this loop so that going from quadratic coefficient to quadratic coefficient becomes as simple as one multiplication (which I think is possible or close to).
4. Then we feed that into process_interval... and if the smooth candidates is reduced enough in bitsize we do trial factorization using the small factor base.

The only thing that matters is finding very large squares quickly to reduce the bitlength of smooth candidates... this is the ONLY way a quadratic sieve type algorithm can punch past that 110 digit ceiling. I'm going to make it happen. Watch as I perform cyber magic now, hahaha. Give me a few days. Nearly there now.... so close.
