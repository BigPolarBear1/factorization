DISCLAIMER: AT NO POINT IN MY 2.5 YEARS OF RESEARCH DID I USE AI OR HAS AI IN ANY WAY HELPED WITH MY RESEARCH. I'M PUTTING THIS HERE BECAUSE I KNOW PEOPLE WILL ACCUSE ME OF IT. 

Note: Remember, as this now comes to completion in the coming days. Both the FBI and Microsoft knew I was working on factorization when I had only been working on it for a few months and nobody else in the entire world knew I was working on it. I tried desperately to avoid the path of public disclosure as I would have much rather sold my work. But at every turn people made it impossible. This would not have happened 10 years ago. An oppurtunity presented itself, it was literally handed to you on a silver platter.. and instead you people all treated me like shit and gambled on me failing with my project. I ask you, where has your common sense gone? I hate what this will do to everyone around me, this is never what I wanted. PS: I know Microsoft's PR teams will say a lot of awful shit about me, completely lacking of context, because that is just the type of people they are. If you want to know who I really am, I suggest you talk to my former teamlead and my former manager who Microsoft retaliated against for defending me. I'm done with this entire shit industry. Everything that happened, you should think whose fault this really is. Maybe instead of morrons sitting on golden thrones build with big tech money, pointing fingers at everyone whose skills and drive hurt their fragile little egoes.. maybe start using your brains, because you people keep losing and keep making things worse.

Note2: I am antifa, leader of the cipherpunk, fuck the FBI department. 

Note3: Once Improved_QS_Variant is finished, I will rewrite most portions of the paper completely from start. I keep finding myself "maturing" out of these papers.. because I keep learning and advancing with my math skills. There's so much stuff in that paper that I just want to completely rewrite because even to me it looks way too amateuristic now.

Note4: I am utterly broke so if anyone wants to make donations to keep this research going a little longer: big_polar_bear1@proton.me email me... thanks.

#### To run from folder "QS_variant" (Standard SIQS with our number theory as backend):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 220 -base 20_000 -debug 1 -lin_size 10_000_000 -quad_size 1</br></br>

See below for an improved way of performing what this PoC does.. I'll delete this Proof of Concept once the PoC for Improved_QS_Variant matures a little more<br><br>

Note: With a large enough -base and lin_size this PoC will find smooths for 110 digits. Albeit very slowly, but this is a highly unoptimized cython PoC. However, to push beyond that into novel terroritory for Quadratic Sieve-based algorithms we need to use quadratic coefficients and p-adic lifting, and that is what the PoC below (Improved_QS_Variant) will be for. 

#### To run from folder "Improved_QS_Variant" (Implements more of my number theory and attempts to succeed with fewer smooths by using p-adic lifting):</br></br>

To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 140 -base 10_000 -sbase 5000 -debug 1 -lin_size 1_000_000 -quad_size 100</br></br>

Update: Bah brain is hugely overactive. So YOU MUST REDUCE THE BIT LENGTH OF SMOOTHS WITH LARGE SQUARES. That I know for sure. There is no other way to punch past 110 digits with a QS style algorithm otherwise. Absolutely no other way. I'm still really thinking how to best approach this. A very simply method would be to just take generate_modulus.. have it spit out much bigger square modulus sizes... if the root(s) is too large we can resize it with the quadratic coefficient. I wonder how practical that would be. Let me try some thing this week. Or some type of variation on this. Can always multiply the modulus with small primes to increase the root size a little bit. Its really the size of that root and the size of the modulus that matters... and I know I can resiez the size of the root with the quadratic coefficient. I just gotta get focused for a few days and finish this the proper way. I see the solution... it's just these fucking details. I hate math. I absolutely hate math. If I could go back in time, I would have become an artist. I like creativity. I don't like this shit. Although I guess its mostly just being broke, not having seen my friends in 2 years (dont even hear from them anymore) .. and just the general hopelessness of my situation. 
