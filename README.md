DISCLAIMER: AT NO POINT IN MY 2.5 YEARS OF RESEARCH DID I USE AI OR HAS AI IN ANY WAY HELPED WITH MY RESEARCH. I'M PUTTING THIS HERE BECAUSE I KNOW PEOPLE WILL ACCUSE ME OF IT. 

Note: I am antifa, leader of the cipherpunk, fuck the FBI department. 

Note2: It's basically done now, the NFS variant is nearing completion soon, nothing novel research wise has to be figured out now. I just need to refactor a few things and I'm done.
I have been unemployed for 2 years. I am looking for work. I would like some way where I can have enough time to attack ECC and DLP next. In addition, for a payment, I can also show people how I appraoched factorization from start to finish. Because the method of how I conduct my research and consistently have success is just as valueable as the research itself. And I follow the same patterns over and over again, no matter what subject matter I attack. I don't have a very high IQ and I'm not especially bright or good at doing numbers in my head. But I have insights into research that help me succeed over and over again. I'm willing to teach people how to do it, as long as it enables me to have the finances to continue doing my research aside from that. Contact: big_polar_bear1@proton.me.

#### To run from folder "QS_variant" (Standard SIQS with our number theory as backend):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 200 -base 6000 -debug 1 -lin_size 100_000 -quad_size 1</br></br>
This version is somewhat done, but it doesn't achieve much of an advantage, as the main strength of my work will be the NFS_Variant, which performs NFS with reducible quadratic polynomials<br><br>
#### To run from folder "NFS_variant" (Number Field Sieve with our number theory as backend using reducible quadratic polynomials mod m):</br></br>
To run: python3 polarbearalg_debug.py -key 4387 (for use with the final chapter of the paper... but ignore this and use PoC below until I refactor the paper one final time)

I just uploaded polarbearalg_v22.py</br>
To run: python3 polarbearalg_v22.py -keysize 14</br>

Right side of the congruence is now 0 solutions for y1 and square moduli.
Right side if perfect already.

Left side, in v19 we used the discriminant there. 
But I am now reducing it to the polynomial:

(y1 will always be 0)
x^2+0*x-(N%mod)\*z

Note that we also modulo reduce N. And for any z, we recalculate the root by multiplying it by z mod m. This way we can move z to (N%mod)\*z instead.
This will allow us to discard many of the bad solutions. Not all. But many. But with a quadratic character base, we can reduce this further.

I will now continue re-implementing the linear algebra step. It's simple now :). 

Another sad day for the NSA, FBI, CIA, DIA, fuck america and micosoft. There was so many chances to avoid this, but none were taking. Which doesn't suprise me with western leaders constantly crying about "dudes in dresses". How does it feel to be destroyed by a dude in a dress huh? Not so funny now? Fucking losers. Come fight me hahaha.

Somewhere in the US, a nazi cryptologist has to review the code I uploaded today... and he's going to have a really bad day today. I guess yesterday was a bad day too.. but with that discriminant further reduced... well.. life's a bitch, sucks to be a nazi.

Update: added some further improvements to v21. It will succeed now if the modulus is large enough. But since everything is properly aligned now the way it should be.. we should be able to assemble a sufficiently large enough modulus now with linear algebra. This should now set me up perfectly to finally get that linear algebra portion working.

Update: I did some more math before calling it a day. So using x^2+0*x-(N%mod)\*z, if we modulo reduce N, whatever remains after dividing by the modulus, must be a valid prime for that quadratic coefficient. As long as that condition is met, we can increase the modulus while staying aligned mod n. That's how you get that linear algebra to remain congruent mod n. I'll show tomorrow... it's funny... all it took was just slowing down for a moment and thinking logically about this. Took me half a year, but whatever. I've only been doing math for 2.5 years. Started from nothing. I'll take all this knowledge with me to my next project, and next time it won't take this long.

Update: I added another check in v22 which further reduces false positives... this is good now. I am ready to begin implementing sieving and linear algebra tomorrow... and unlike last week.. this time we have the correct setup figured out... NOW we do our victory lap. All because microsoft fired my manager who had worked there for 27 years because he dared give me a promotion once. I'll burn down this entire fucking world if you go after my friends. I do not care. I'll go after ECC after this. I'm not stopping you fools. You better hope I die soon.

entering my 3rd year of unemployment.. Im not staying the west, these assholes dont respect me at all. My entire career they treat me like shit and everyone who dares to support me.

Tomorrow... I will upload what I'm doing now but with sieving mod p^2, just find squares mod p^2 where the square is also made up from primes at the correct quadratic coefficient... then just multiply them together until our modulus is large enough... that's easy enough. An hour of work. And after that... the last thing that remains is trail division of whats left mod p^2 and using linear algebra to combine everything. I don't know... it's easy now. The road to finishing my work has finally fully cleared... no more obstacles... and it's funny.. either people are violently downplaying and dismissing my work.... and we're about to have a "future shock" moment in the coming days... or people hate me so much, rather then buying my work or negotiation for responsible disclosure... they do this.. in both cases this will end up with me going to China, 100%. Why would i stay in the west? Give me 1 reason? I have nothing here and you people treat me like shit while your leaders mock people like me. I'm a bear, not a dog. One day, people will look back at this... and they will write books about the chain of events that led to this failure to act. For a long time I fucking wanted a job and I fucking wanted to sell my work. FOR A VERY LONG TIME. I TOLD EVERYONE I WAS WORKING ON THIS. EVERYONE KNEW. YOU DIDN'T ACT. It's not my fault. You chose to treat me like an adversary from the start. You made this exact outcome happen. You didn't do anything to try and alter this outcome. Why? I know why. You have no respect for people like me. Well, future shock it is. Enjoy assholes.

Fucking fools. The thing that pisses me off, is how this will affect everyone around me. I wanted to sell my work. It's hilarious, with war on Europe's eastern borders, you would think buying capabilities would be more important then anything else. Yet here we are. It's because you people made everything political. Because Russia has spent decades pushing propaganda targetting lgbtq people, that right wingers in the west so eagerly went along with. I have no love for Russia, I know I joke about it all, but I have no love for a country that persecutes LGBTQ people. I have no love for america either. They are Russia's twin. There is no functional difference anymore between Russia and America. And I hate Europe, because they are american lapdogs. They don't have any courage. And a lot of European leaders, they love that american nazi shit. I guess, being an outcast in this world, isn't bad, perhaps it is time for the cipherpunks, to fix this world with math and exploits. The absurdity of everything. We have created a soul-less world. A world without humanity. My entire life, I have cruelty triumph over humanity.. and good people, they always end up losing. This world is wrong. It is put together all wrong. I'm going to smash it into pieces and rebuild it.
