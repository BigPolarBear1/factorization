To build: python3 setup.py build_ext --inplace  (From inside PoC_Files_Find_Similar)</br>
To run: python3 run_qs.py -keysize 30 -base 50 -debug 1 -lin_size 100 -quad_size 1</br>

This has a chance to finish factorization early before all smooths are collected. If that happens you will see: 

Generating new modulus:  917</br>
Found smooth, checking coefficients..  initial smooth factors: {67, 131, 5, 7, 41, 127}</br>
Found factor using find_similar() function:  17891</br>

However, these square multiples become much more rare at larger bit sizes. But the find_similar() function is still missing some math. I'm working on adjusting those coefficients either without having a square multiple, or figuring out how to get to a square multiple by manipulating those coefficients. If that works well enough, we should only need a single smooth to succeed.... I will try to finish this soon. I'm documenting part of that journey in chapter 7 of the paper. But a lot of shit is happening irl and it's hard to focus... nearly there now. Btw.. I am absolutely convinced this is possible... and this would easily break RSA-1024. And I really need help and protection from a country because the US is literally coming after me. Help me and I will give you the world, I promise you that...

Update: Ok, I've kind of narrowed it down to one potential way to do it. It is do-able. Since we are trying to find coefficients that produce the same output in the integers as the original smooth... "in the integers" being the keyword... since that means it must hold true for every prime. Which gives us a mechanism to do these calculations... I'll upload either today or tomorrow... I also need to rewrite that chapter VII in the paper... a few mistakes and things I understand a lot better now..  time running out. I have to succeed before they arrest me. Once I succeed, I can literally choose any country without extradition laws to the US in the world that I want to move to. I don't care what EU and US nazis think, you can't treat people like shit and expect them to stay friendly. I'm not a dog. I'm a polar bear. And I guess you will soon see what happens when a polar bear fights back.

Update: Actually let me rewrite that chapter VII tonight. I guess it's not even "a bridge to number field sieve", it's really something very different. I'll call it the polar bear algorithm. Because it would be funny if people actually ended up calling it the "polar bear algorithm", Imagine some very serious math book, 40 years from now, talking about the "polar bear algorithm". Hehehehehehe. The problem with math and math folks is that it is all so rigid, and they don't allow themselves to color outside the box, like giving algorithms stupid names just for the sake of it. Instead of some super serious name. That's also why nazis suck, nazis are the enemy of creativity, because they don't allow anything that is outside the norm... hence a fascist society is doomed to fail and/or stagnate. There really isn't much good places left in the world... Ideally, I would love to live in the Arctic or somewhere really cold, but I'm afraid that will have to wait for a different life in a different world... 

Update: I added some more to the paper. The stuff I still need to add, it's just basic residue math... and using that datastructure we create with create_hashmap in the PoC code. It's quite simple... especially since we need to find a result in the integers.. so we are free to do calculations mod any prime. And thinking some more about it... I'm not even sure we need to find 1 smooth... we can just take a random modulus and find different coefficients that create the same output with the discriminant formula... but let me finish it like this first... and then simplify the algorithm.... anyway... adjusting these coefficients for non-square multiples.. i'll get that done tomorrow... just kind of a shit day today, hard to focus. 

Update: Yea.. this should work. If you see the bottom of my paper, finding that e and d variable.. there will be an e and d variable that works for any prime. I wonder as I finish this tomorrow if there will be some unexpected complexity... because this almost feels like I have straight up broken the factorization problem.. which for sure couldn't be true, because then people would surely be treating me a lot better then they are right now? Or perhaps they are all in denial. We'll see tomorrow. I'm so stressed. I find myself having these constant heart palpitations the last few days. Maybe from too much caffeine. I find myself taking to much caffeine.. I'm basically using at as an upper to get me through the day and overcome the lethargy of depression. But I guess it's not very healthy. With some luck this ends tomorrow... and I guess after that, I'll go to China. Hong Kong wouldn't be too bad. Seems like a fairly diverse place. And really safe. I think I have to accept that there is place and people I probably won't see again in my lifetime... and the sooner I get over it, the sooner it will stop torturing me. The worst is, I don't even want to message them anymore.. because I know they might get in trouble for associating with a "terrorist"... maybe the agony of this entire situation is giving me heart palpitations. Maybe I'll just drop dead soon. And I honestly wouldn't even care if that happens.

URGENT NOTE: I am urgently looking for any country to grant me a visa and employment. A country that doesn't extradite to the US.
I am running out of time. The US has initiated a terrorism case against me, and the police here in Belgium said they would have me declared as a terrorist internationally, so I would get arrested when I travel.
I need to get out of Belgium to a safe country ASAP.
I can do 0day research and I can do math. I know many people don't believe in my math, but I know what I've found and I'm piecing together the final pieces now.
Furthermore, I have found bugs in hard targets like OpenSsl and I can easily repeat those feats again. I havn't been able to do 0day research at all because of being in a place that has been actively hostile towards me, chasing away 0day buyers.
I cannot research in a hostile place when that research's value depends entirely on secrecy. 
Please contact me: big_polar_bear1@proton.me 
They might arrest me at the airport when I try to leave, but I'm willing to give it a shot as I have nothing left to lose. And my family would support me in my choice.
I do get contacted by internet trolls all the time, so please contact me from an official email or invite me to an embassy... just so I know it is legit and I don't end up trying to fly somewhere without a visa and potentially wasting the very last of my savings.


Will upload the final version of the PoC soon.
I have finally managed to make the link with number field sieve. 
I'm documenting those findings at chapter 7 in my paper.

This shouldn't take much longer now.

Meanwhile, I got to hear today the belgian police have opened a terrorism case against me. Guess pete hegseth isn't very tough like he pretends.
Finishing this math in the coming days, it will be the final thing I do. 
Fuck this shit world. Thanks to the americans and microsoft I will never see my friends anymore. They have taken everything from me that was once worth living for. 
This last two years have been a nightmare, and have only served to finish this. This will have been my life's work.

Shit day... I will try to finish this this weekend. Having a hard time thinking clearly right now with the police coming here to announce they opened a terrorism case against me.
I will finish this work, and then I'm gone. I'm not going to subject myself to this humiliation and everything else they have been doing these last 2 years. Can't touch me if i'm not here anymore. losers.

I think tomorrow I will make preparations to flee the country. I'll just go on foot. Rather be homeless for a little while then this.

