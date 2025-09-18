# factorization

To do:

Currently the PoC is able to pull smooths from multiple quadratic coefficients (defined by quad_sieve_size). 
I'm refactoring the PoC so we just pull from a single small quadratic coefficient and when we find a smooth, we call into a new function.
In this new function we check the smooth's factors (only the ones with odd exponents), we then find other quadratic coefficients where the factors are a quadratic residue (if calculated using N\*4\*quad_co).
Here we then calculate linear coefficients and build a sieve interval and check for smooths.
Further more we also do p-adic lifting.. so we can find smooths with just the factors of the original smooth. Or close enough... so that if we can find a bunch of them we are garantueed a square relation.


Update: Got like an hour of work done today... but got a first very rudimentary PoC... it works. Yea... this is going to work for sure. lol. I feel such sadness and anger. They took everything from me. I used to be surrounded by the best friends a bear could possible have. And they took that from me. After threatening me with a gun and shouting slurs at me. And microsoft portraying me as incompetent. And then.. a year ago.. when I dropped the first iteration, the crappy one using LLL... people would have known my work was atleast headed in the right direction. BUT NOBODY DID ANYTHING. NOBODY. You chased away 0day buyers. Made it impossible for me to find a job. Have an income. You pushed me to the brink of hopelessness. You not just took away years from my life. You took away my joy, my happiness, and any color this world once had. I don't care anymore. My warpaint is on now. I am going to fuck this world to hell until my last dying breath now. It works. It fucking works. I was right. 

Update: I'll upload a first rough draft of what I'm talking about tomorrow... it's time to finish this. Get this over with. I've bled for this, I've worked harder then anyone in this fucking world ever has. And my whole career, they treat me like fucking shit. And they treat the people who support me like fucking shit. Let it all burn. LET IT ALL FUCKING BURN THEN.

Update: Right now I've implemented in code (not uploaded yet):

1. Create a new function that we call into when a smooth is found
2. Find other quadratic coefficients where the smooth factors (with odd exponents) are a quadratic residue
3. Iterate all possible linear coefficients for the modulus,which is jsut the factors of your smooth multiplied together (once for each factor with odd exponent) and do so for all viable quadratic coefficients.
4. Do trial division using just the smooth factors. We could use the entire factor base.. but we're mostly interested in just finding a square relation with as few smooths as possible.

To do:

Now what I still need to write is creating a sieve interval... where each step is the modulus (the factors of the smooth)... and we use p-adic lifts of these smooth factors which make up the modulus to fill out the sieve_interval.
That's the crucial trick. We need those p-adic lifts of our coefficients to make this really efficient. But luckily we can calculate those insanely fast.
This shouldn't take much longer... maybe one more day... I'll go for a run now though... just massive brain fog. Hard to focus. And really starting to struggle with depression. Just can't cope with this shit anymore. And it's all so obvious, like I don't see what is happening. If this didn't affect my life so much, it would be pure comedy.

Update: Back from running. So agitated. I hope nobody gets in my face when I feel like this, because I feel like I'll end up in jail. I guess that's what trauma and being treated like shit does to a person. I'll get some food and shower and write some skeleton logic for that last step so it's 100% done tomorrow... this ends tomorrow, I don't care, it's been dragging on for too long now. It ends tomorrow. And I know this isn't some "gimmick" or one of those ideas that will "maybe work"... this will work for sure. This will work 100%. This is it. This is what I needed. My breakthrough. If I can succeed with only a fraction of smooths..... lol.... do you know how big numbers I can factor? RSA-1024 definitely about to die for sure. And after that, I am not stopping. THIS WILL NEVER END UNTIL IM DEAD. I'LL GO AFTER ECC TOO. YEA, SAY GOODBYE TO YOUR CRYPTO YOU FUCKING WANKERS. Crypto grifting nazis pieeces of shit. 

Update: Always knew the british were a bunch of cunts. Nazi collaborators.

America truely has fallen. Just like that. 

It's hilarious. All people had to do was give my former manager a bunch of money for unjustly firing him and I would have been perfectly fine selling my work. I tried. I'm guessing people have gotten so good at bending the truth, they can't see truth anymore when it's right infront of them.... oh well... time to start coding this final piece. I wanted nothing more then a normal life, with a job, employment, stability... you denied me everything, for 2 years, I couldn't get anywhere in life, despite you people knowing (or you should have known) my work was solid... and I know it was some nazi pieces of shit who forced me down this path. Because they can't let a transgender person have any wins in life. Yea, you could call me radicalized now, because I will destroy you people for what you did. I will become the terrorist antifa. HAhahahahaha. Fuck nazis. Nazi losers. I despise everything about nazis. You all keep saying you people are somehow "better" or "superior" ..... but you are just shit for brains losers.

Update: Let me just clean up what I got and upload it and go for a run. I still need to implement the p-adic lifting.. but I am having a particulary bad day. These last few months I constantly feel like I'm right at the threshold of having some type of nervous breakdown. Like, just going into the woods and hanging myself or some shit. It's really hard. It's really hard constantly having to think back to people and places I can't see again. While stuck in this really small room. Just getting this really intense claustrophobia from being stuck here for 2 years. Let me clean up this code and upload it and go for a run.... it demonstrates what I'm trying to do already. And I guess it's basically the same as disclosure as anyone can just implement the p-adic step themselves. Ofcourse the question is... will anyone admit I'm right or will they just keep denying everything right until I set a factoring record? Or are they hoping I kill myself before that happens? Because that is what this looks like. This is 100% what this looks like.

Update: I have uploaded my work in progress:

Build: python3 setup.py build_ext --inplace</br>
Run: python3 run_qs.py -keysize 25 -base 50 -debug 1 -lin_size 10 -quad_size 100_000

Just run on 25-bit numbers, this can't possibly work on bigger numbers until I implement p-adic lifting and sieve_intervals... but it demonstrates the idea.
Code is also very sloppy. Going through a period of particulary bad brain fog. I'll fix it soon. 

You should see output like this:</br></br>

Constructing quad interval took: 2.091426908998983</br>
Processing quad interval took: 0.09593547700205818</br>
Generating new modulus:  6815</br>
Looking for similar smooths, amount needed: 8 initial smooth factors: {5, 43, 47, 19, 29, -1}</br>
Found a similar smooth {2, 5, 47, 19, 29}</br>
Found a similar smooth {2, 43, 47, 29}</br>
Found a similar smooth {43, 29, -1}</br>
Found a similar smooth {2, 29}</br>
Found a similar smooth {5, 43, 19}</br>
Found a similar smooth {5, 43, 47, 19, -1}</br>
Found a similar smooth {5, 43, 19}</br>
Found a similar smooth {43, 47}</br>
sim_found:  8</br>
Running linear algebra step with #smooths:  9</br>
[SUCCESS]Factors are: 2503 and 2777</br></br>

So with a factor base of 50, here we are able to succeed after finding only 9 smooths (pulled from different quadratic coefficients).
Which is a fraction of what would otherwise be required. 

But yea, this will not work well until I implement p-adic lifting and sieve intervals. Otherwise it just becomes too difficult to find these "similar smooths" at higher bits. 

The cool thing is, if you find a smooth with only 1 or 2 odd exponent factors... you only need to find a couple of similar smooths to succeed. Ofcourse they will be hard to find, but with p-adic lifting, it is possible. You could in theory factor a 200-bit number with only 5 smooths or something this way. 

I guess when I get back from running, I'll add a sieve interval. Thats like 10 minutes of copy paste work. Step size will just be the modulus (made up from the smooth factors). Later we can adjust this step size based on the quadratic coefficient (i.e multiply the modulus be a smooth factor, or divide by one, etc). So get that sieve interval shit out of the way today... then tomorrow, I'll add the p-adic lifting.... then it's done. And then I need to fix that sloppy brainfog code. God damn. I don't know what's wrong with me. Maybe from running ultra marathons. Or maybe the CIA poisoned me. Definitely something those cowards would do. Fucking CIA pussies. I would beat any CIA fucker in a fist fight. Fucking shitheads.

Update: Quickly added some sieve interval code for the find_similar family of functions. It's all ready now for tomorrow to add the final missing piece (and most crucial piece) ... p-adic lifting. That's the engine that will make this entire approach work. This will work. I know it will. So tomorrow will be the moment of truth.... and if it works... there is no going back. Endless others way this could have gone down. It is funny how circumstances have conspired to have it happen this way... oh well.

