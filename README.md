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

Update: After a little more experimenting, I think I need to approach it differently. I need to be smarter about this then just brute force trying to reduce the bit size of smooth candidates with large squares. We can start lets say zx^2-N=a^2    ... or zx^2=a^2+N  ... hence we can calculate x and z for a given a^2. The goal being to get a small z that factors over a factor base. We can also then just multiply a^2 to try and reduce the z value. Aka... shift the parabola... 

I'll mess around with this a while. I need to have a really focused appraoch. Nothing about my algorithm should be bruteforce. No sieve intervals... nothing of that sort. It needs to be done in a smart and elegant way. I'll explore this angle. Basically we just need a large square to reduce the bit size of smooth candidates. I can start with a large square... and then just add factors... then the goal is simply to generate a small matching z value that factors over a factor base. So if I appraoch it like this... basically building smooths and finding smooths with a small quadratic coefficient (z value) ... I think that is a much stronger approach that makes more sense. 

Anyway, I have to run to court ordered therapy tomorrow for sending angry emails to the FBI (because they are pathetic losers), I'll run 30k to and from... so wont make much progress tomorrow... but lets see... I hope to have something finished this week. I guess in the last few days I did work out quite a bit of the math how to manipulate and recalculate those roots for different quadratic coefficients... so we dont need to build a huge hashmap anymore like we did before. And I am zero-ed in on the only approach that can possibly work.. and that is REDUCING THE BIT SIZE OF SMOOTH CANDIDATES WITH LARGE SQUARES. All that remains now is just implementation details... just got to be smart about it. Thats the problem with algorithms, everything needs to be deliberate and smart... you cant have any bruteforce in it, because thats always going to scale poorly. Every number manipulation needs to bring you closer to a solution.
