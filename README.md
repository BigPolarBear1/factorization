DISCLAIMER: AT NO POINT IN MY 2.5 YEARS OF RESEARCH DID I USE AI OR HAS AI IN ANY WAY HELPED WITH MY RESEARCH. I'M PUTTING THIS HERE BECAUSE I KNOW PEOPLE WILL ACCUSE ME OF IT. 

Note: Remember, as this now comes to completion in the coming days. Both the FBI and Microsoft knew I was working on factorization when I had only been working on it for a few months and nobody else in the entire world knew I was working on it. I tried desperately to avoid the path of public disclosure as I would have much rather sold my work. But at every turn people made it impossible. This would not have happened 10 years ago. An oppurtunity presented itself, it was literally handed to you on a silver platter.. and instead you people all treated me like shit and gambled on me failing with my project. I ask you, where has your common sense gone? I hate what this will do to everyone around me, this is never what I wanted. PS: I know Microsoft's PR teams will say a lot of awful shit about me, completely lacking of context, because that is just the type of people they are. If you want to know who I really am, I suggest you talk to my former teamlead and my former manager who Microsoft retaliated against for defending me. I'm done with this entire shit industry. Everything that happened, you should think whose fault this really is. Maybe instead of morrons sitting on golden thrones build with big tech money, pointing fingers at everyone whose skills and drive hurt their fragile little egoes.. maybe start using your brains, because you people keep losing and keep making things worse.

Note: I am antifa, leader of the cipherpunk, fuck the FBI department. 

#### To run from folder "QS_variant" (Standard SIQS with our number theory as backend):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 220 -base 20_000 -debug 1 -lin_size 10_000_000 -quad_size 1</br></br>

See below for an improved way of performing what this PoC does.. I'll delete this Proof of Concept once the PoC for Improved_QS_Variant matures a little more<br><br>
#### To run from folder "Improved_QS_Variant" (Implements more of my number theory and attempts to succeed with fewer smooths by using p-adic lifting):</br></br>

To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 140 -base 10_000 -sbase 4000 -debug 1 -lin_size 1_000_000 -quad_size 1</br></br>

Alright, doing some big reductions in complexity and in the process of rewriting the high level approach of Improved_QS_Variant.
If you run the above command, it should finish somewhere between 300-500 smooths. This is because we only mark the sieve interval with odd exponents if the prime is less then the -sbase value.

The code is also already in place to check multiple quadratic coefficients, but it wont find much smooths on those without resizing the modulus. Which I dont want to waste time on as I'm going to rework all that code next.

So to factor really large numbers (400-bit+) .... what we need to do is reduce the smooth candidates with very large squares. And we can go looking for these very large squares at different quadrati coefficients. So this week, this is what I will start implementing. All the math is there already. It all works. I just need to rewrite some portions now and we're done....

Sad day for the NSA cryptologists. Not that I feel bad for them, because they are probably a bunch of nazi pricks who hate people like me.
