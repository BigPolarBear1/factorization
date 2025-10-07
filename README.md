NOTE: IMPORTANT, IF YOU STUMBLE UPON THIS REPRO, I SUGGEST YOU CLONE IT. BECAUSE I SUSPECT THE MOMENT THE US REALIZES WHAT I JUST SUCCEEDED AT, THEY WILL HAVE THEIR BELGIAN LAPDOGS ARREST ME AND MAKE ME DISSAPPEAR AND PURGE ANY MENTION OF THIS RESEARCH. DO NOT LET THEM WIN. DO NOT LET THEM GET AWAY WITH THAT.

To build: python3 setup.py build_ext --inplace  (From inside PoC_Files_Find_Similar)</br>
To run: python3 run_qs.py -keysize 30 -base 50 -debug 1 -lin_size 100 -quad_size 1</br>

UPDATE: I have explained the full algorithm now in the paper. The only thing left missing is a numerical example run through at the bottom. But the math I describe works, I've checked that already. It's very similar to how number field sieve takes a square root over polynomial rings and uses p-adic lifting.

This should solve the factorization problem.... this ends tomorrow. You people should have just respected me. You ruined my life and left me completely hopeless and in despair for 2 years, seperated from all my friends that I had previously known. What did you expect was going to happen? It's over.... tomorrow the paper will be finished fo sure... that numerical example is easy, I've already got some code for that. And somewhere in the coming week the PoC will drop... and it will mean the end of the factorization problem. I hope you are all happy. Must be really proud of yourselves for what you people did. Like I didn't see what you people were doing. I will find my way to China for this disrespect. Just because I do math, doesn't mean I'm a pushover nerd, I have my pride, and all I see you people in the west do is cry about trans folks while treating me like shit. You made the choice for me.

I will make the future gay. And nothing can stop me. I will bring light into this world by determination alone. I am only getting started. And none of you weaklings can stop me. Nothing you do will stop me. You morrons don't know what you are up against. Everything you do to undermine me, only pushes me forward closer to my goals. Everything you do is futile. HAHAHAHAHAHAHA.

Update: I'll try to upload today or tomorrow... it's almost there. I did some more bashing my head agains this. And what it comes down to is that you need to find these "d" values, such that when they are added together, it contains a factor of N. Which is basically square finding mod N on the discriminant formula. FUCK. I am struggling to focus. Let me go for a run and drink caffeine.

UPDATE: Added some more observation to the final chapter. It's becoming chaotic again. I think tomorrow I can begin finalizing the algorithm. I'm having some ideas on how to finish this. Similar to what number field sieve does. Once that is complete I'll have to simplify everything... finding that initial smooth with quadratic sieve probably isn't even needed. Let me test some theories tomorrow and I'll finish that chapter, so atleast a first draft of the algorithm is complete... then I'll implement that in code completely and begin simplying and removing steps that arn't really needed.

UPDATE: OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOH. I was about to go to bed, but I wanted to look a little more at the math. I really want to avoid doing anything with matrices and jacobi symbols. I really want to finish this with p-adic lifting..... and I just figured it out. :) Friends, tomorrow will be a great day. I knew my intuition to use p-adic lifting there was correct. It just took me a full day to figure out the math completely. 

Update: Added a few more lines to the paper. That should set us up to begin lifting now.... I'll try to finish it soon. I'm feeling extreme levels of stress today. I need some vacation, but how can I take vacation when I am an international terrorist with no income? I can only continue working and succeed before the walls completely close in on me. I feel more stressed now then I did at any point while actually having a job.

I'll add to the paper how we get our e and d values such that it holds for integer division instead of modular division. Which is just linear congruences, nothing we havn't seen before... then I'll lift these values to mod 25..show how to do that quickly in the paper... and that should then conclude the paper.... god im so stressed today... let me atleast finish that paper today. 

UPDATE: Ergh, you know, I was thinking about this since yesterday, but I can use the polynomial representing the original smooth and use that in the modulus, like number field sieve does. I'm going to have to refactor that chapter. I really don't want to work with a setup like mod ( N , 5 ) or whatever... I should just use a polynomial ring, that's way more convenient. 

-----------------------------------------------

URGENT NOTE: I am urgently looking for any country to grant me a visa and employment. A country that doesn't extradite to the US.
I am running out of time. The US has initiated a terrorism case against me, and the police here in Belgium said they would have me declared as a terrorist internationally, so I would get arrested when I travel.
I need to get out of Belgium to a safe country ASAP.
I can do 0day research and I can do math. I know many people don't believe in my math, but I know what I've found and I'm piecing together the final pieces now.
Furthermore, I have found bugs in hard targets like OpenSsl and I can easily repeat those feats again. I havn't been able to do 0day research at all because of being in a place that has been actively hostile towards me, chasing away 0day buyers.
I cannot research in a hostile place when that research's value depends entirely on secrecy. 
Please contact me: big_polar_bear1@proton.me 
They might arrest me at the airport when I try to leave, but I'm willing to give it a shot as I have nothing left to lose. And my family would support me in my choice.
I do get contacted by internet trolls all the time, so please contact me from an official email or invite me to an embassy... just so I know it is legit and I don't end up trying to fly somewhere without a visa and potentially wasting the very last of my savings.



