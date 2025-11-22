DISCLAIMER: AT NO POINT IN MY 2.5 YEARS OF RESEARCH DID I USE AI OR HAS AI IN ANY WAY HELPED WITH MY RESEARCH. I'M PUTTING THIS HERE BECAUSE I KNOW PEOPLE WILL ACCUSE ME OF IT. 

Note: Remember, as this now comes to completion in the coming days. Both the FBI and Microsoft knew I was working on factorization when I had only been working on it for a few months and nobody else in the entire world knew I was working on it. I tried desperately to avoid the path of public disclosure as I would have much rather sold my work. But at every turn people made it impossible. This would not have happened 10 years ago. An oppurtunity presented itself, it was literally handed to you on a silver platter.. and instead you people all treated me like shit and gambled on me failing with my project. I ask you, where has your common sense gone? I hate what this will do to everyone around me, this is never what I wanted. PS: I know Microsoft's PR teams will say a lot of awful shit about me, completely lacking of context, because that is just the type of people they are. If you want to know who I really am, I suggest you talk to my former teamlead and my former manager who Microsoft retaliated against for defending me. I'm done with this entire shit industry. Everything that happened, you should think whose fault this really is. Maybe instead of morrons sitting on golden thrones build with big tech money, pointing fingers at everyone whose skills and drive hurt their fragile little egoes.. maybe start using your brains, because you people keep losing and keep making things worse.

Note: I am antifa, leader of the cipherpunk, fuck the FBI department. 

#### To run from folder "QS_variant" (Standard SIQS with our number theory as backend):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 160  -base 4000 -debug 1 -lin_size 100_000 -quad_size 100</br></br>

See below for an improved way of performing what this PoC does.. I'll delete this Proof of Concept once the PoC for Improved_QS_Variant matures a little more<br><br>
#### To run from folder "Improved_QS_Variant" (Implements more of my number theory to also take advantage of quadratic coefficients):</br></br>

To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 140  -base 4000 -debug 1 -lin_size 1_000_000 -quad_size 10</br></br>

Update: Alright going to call it a day. Still havn't managed to optimize what I wanted to.
However, in create_hashmap it only has to calculate linear coefficients for quadratic coefficient = 1 now, so that will allow us to generate enourmous factor bases quickly with little memory trade-off. In find_similar() we calculate the linear coefficients for different quadratic coefficient on the fly. And in construct_interval_similar() we now do the same, we precalculate those coefficients mod p in gen_co2 now and then derive the linear coefficiens and roots for different quadratic coefficients from that on the fly. By doing it on the fly... we can make the quad_size variable much bigger... as it wont translate into more memory useage.

Now tomorrow, to begin seeing the really big boost..I need to move lift_root and solve_lin_con out of construct_interval_similar into gen_co2 .. so its only calculated once, we we derive what that should be on the fly. That should be an enormous speed gain. When that is done, all that is left is just code optimizations. I may precalculate some tonelli stuff aswell if that ends up bottlenecking.

So the whole gist of this vs default SIQS is that we can finish much earlier, and also, constructing sieve intervals will be much faster, since sieving in the direction of the quadratic coefficients is very easy as we can derive everything we need with simple calculations from roots and coefficients for quadratic coefficient = 1. So once the optimizations are done. It will sieve faster then normal SIQS plus it can finish faster then normal SIQS. It will be better in every way. But you'll see. It should begin getting there tomorrow as I move more stuff out of construct_interval_similar()

Update: I just did some thinking. What I should do first thing tomorrow is in QS_Variant (not the improved one) ... just have it very quickly build sieve intervals for different quadratic coefficients. Because not having to deal with lifting, these calculations can all be done instantly to derive roots and coefficients for different quadratic coefficients, lifting is kind of a pain I guess...and its going to be easier to work out of I fully implement a non-lifting PoC first. Then once that is done, that should already see a big speed boost, I can continue what I was doing in Improved_QS_Variant.

Update: I'll port these findings over to QS_Variant.. so it checks the quadratic coefficients too. I'll get that done tonight. Then tomorrow I'll do the remainder of the optimizations in QS_Variant first, since it wont have p-adic lifting.. and once that is done I finalize Improved_QS_Variant aswell.

Update: Ok, i've ported my worked to the QS_Variant aswell. It did slightly slow things down... but tomorrow when I do those optimizations it should gain a performance advantage.
