# factorization

The PoC under PoC_Files will easily factor above 200 bits. </br></br>

However, that one is basically just using the number theory from my paper to achieve default SIQS. </br></br>

in PoC_Files_Find_Similar_WIP I am attempting to sieve in the direction of the quadratic coefficient and linear coefficient. This way in theory we need much less smooths.</br></br>

For PoC_Files_Find_Similar_WIP:</br>
Build: python3 setup.py build_ext --inplace</br>
Run: python3 run_qs.py -keysize 30 -base 100 -debug 1 -lin_size 1_00 -quad_size 100_000</br></br>

You will see some output like this:</br></br>

Constructing quad interval took: 5.059455415990669</br>
Processing quad interval took: 0.16657884302549064</br>
Generating new modulus:  681</br>
Looking for similar smooths, amount needed: 6 initial smooth factors: {2, 227, 59, -1}</br>
local_factors:  [2, 59, 227]</br>
Looking for similar smooths, amount needed: 9 initial smooth factors: {2, 3, 227, 41, 79, 17, -1}</br>
local_factors:  [2, 3, 17, 41, 79, 227]</br>
Looking for similar smooths, amount needed: 8 initial smooth factors: {2, 227, 229, 7, 11, -1}</br>
local_factors:  [2, 7, 11, 227, 229]</br>
Looking for similar smooths, amount needed: 8 initial smooth factors: {3, 227, 7, 103, 43, -1}</br>
local_factors:  [3, 7, 43, 103, 227]</br>
Found a similar smooth {43, 2, 3, 227}</br>
Found a similar smooth {2, 3, 7, 43, -1}</br>
Found a similar smooth {2, 103, 7, -1}</br>
Found a similar smooth {3, -1}</br>
Found a similar smooth {3, -1}</br>
Found a similar smooth {3, 7}</br>
Found a similar smooth {227}</br>
sim_found:  7</br>
Running linear algebra step with #smooths:  11</br>
[SUCCESS]Factors are: 17203 and 18211</br>

If it finds a smooth, it will look in the direction of  the quadratic coefficient to find more similar smooths. And if it finds enough of them, it skips to the linear algebra step. You can succeed with much less smooths this way.
All the basic functionality is kind of there. But it won't work well yet.

I realized I need to address the following things:

1. Our very first step, right after finding the initial smooth should be performing p-adic lifting on the linear coefficients modulo the smooth factors (with odd exponents) using all possible linear coefficient solutions.
   Then we look at those solutions which lift very high while keeping a small coefficient (or dist1, dist2 values.. I need to think what to use as modulus.. but the idea is the lift them in advance so we can just find quadratic coefficients where they will cluster together, while minimizing the linear and quadrtic coefficient) and we find at which quadratic coefficients they occur. Only then should we build the sieve interval (or possibly we could straight up calculate everything without sieve interval). The way it's doing it now is not going to work. Because we don't want to limit our selves to a small quad_sieve_size due to precalculating the primes that are quadratic residues for each. That's a wrong appraoch. We need to begin with p-adic lifts of all possible linear coefficient solutions... and from those results we derive quadratic coefficients which will be favorable. That's how you should do it. Not like it is curently done. Give me a few days. Just really struggling mentally I guess. Depression, being broke.. unemployed and hopeless for years does that I guess. Havnt seen my friends in almost 2 years either thanks to microsoft.

2. With that out of the way... it should be much better... then I need to fix a lot of other shit and shit code. Had some particulary bad brain fog when I wrote this find_similar code. I need to rework literally everything.  It's getting there. The basic building blocks to achieve it are there... it's just getting the implementation right now.

I'll get there eventually... the path is clear now, and I know 100% that this will work. I have 0 doubt. It's just a matter of figuring out how to best do it. But starting with p-adic lifts after finding the initial smooths and going from there.. I think that will already bring us much closer. Hopefully tomorrow I can have a productive day to refactor it to using that approach already.

Let me try and get some barebones functionality ready tomorrow where it does p-adic lifting first after finding a smooth and from those results derives favorable small quadratic coefficients. I have figured out all the math to do this already. I don't know, this month is a slow month, even though I'm at the end of it now. My enemies are trying to wear me down on purpose by denying me any chance of employment or income, or future. Or seeing my friends again ever. They are hoping I will come crashing down mentally before I succeed.. but I will succeed.. and besides.. it is far to late... everything is here already. Anyone can pick this up now and finish it. It's all presented on a silver platter here. And everyone be damned who is participating in this cruelty while knowing, for over a year now, that I am doing everything right. There will be a special place in hell for you people. Wrong side of history. You never want to be on the wrong side of history. Because how we are remembered... that is all there is when our time comes.

It's 100% the americans and their EU lapdogs orchestrating this. Anyway, will move to China if a visa can be arranged through an embassy (i just cannot trust strangers emailing me, sorry, it has to be official, there is to many bad actors and trolls after me all the time) ... or hell, I will move to Israel and make them the strongest cryptologic nation in the world hahahahaha (as long as they are not america's lapdog). I don't care at this point. I will go where people want me, and respect me, and I am ready to go now. Just contact me: big_polar_bear1@proton.me ... and I would urge you to do it sooner rather then later.. because later will be too late. Because these people, they will keep doing this until I upload a PoC to github that sets a factorization record.. and I konw I'm closing in... but I want a fucking future, be respected. I don't want what these people are forcing me to do, because it will ruin me, my family and my friends. And that's the god damn truth of it.

Actually, I will move anywhere, and renounce my european citizenship, as long as it is a country where they don't persecute lgbtq people, that's my only condition. Aside from that, I literally don't care. At this point, the only thing I have left in life that I want, is respect. That is all. I just want respect and you can have my work until I die. Simple as that. Anyway... tick.... tock.... time is running out..... once the final PoC comes online... it will be to late, because then my life will be so royally fucked I might aswell just hang myself to save everyone around me the misery that is about to be unleashed. But that's what they want to happen, huh? I see right through you. Fuck europe for being complicit in this. Fuck belgium most of all. They are no better then the nazi americans. I'm ashamed to be Belgian.
