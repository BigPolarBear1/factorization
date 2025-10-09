To build: python3 setup.py build_ext --inplace  (From inside PoC_Files_Find_Similar)</br>
To run: python3 run_qs.py -keysize 30 -base 50 -debug 1 -lin_size 100 -quad_size 1</br>

Update: Oops. Wait a minute... now that I figured out how those polynomials come into play... I also just realized this whole setup is basically number field sieve. Ofcourse it would be number field sieve on steroids... as it fixes a lot of number field sieve's bottlenecks... I'm overcomplicating things..... OK. OK. OK. Give me a few days. All the math is here already. Check the last chapter in the paper. It already explains how to bridge my quadratic sieve appraoch to number field sieve's appraoch. It's basically the same thing. But to finish the algorithm, you need to thus also use number field sieve's appraoch. You don't need to find an initial smooth using QS. You simply use the math I explained there to achieve the same thing number field sieve does... without having to worry about irreducible polynomials and all that stuff.

Somewhere out there some nazi cryptologist starting to sweat profusely now. Hey, I'm a polar bear, and I'm here to devour all your secrets. HAHAHAHAHHAHAHAHHAHAHAHAHAHHAHAHAHAA.

The american cryptologic program, defeated, by a single polar bear whom your leaders refer to as "a dude in a dress" hahahahahaha. FUCK YOU LOSERS. 

Update: Added a few more sentences at the bottom of chapter VII, to explain how exactly this bridges number field sieve. It's all easy now. I understand it perfectly. I finally have all the pieces of the puzzle THANK GOD, HOLY FUCKNG SHIT. I WAS LOSING MY FUCKING SANITY WORKING ON THIS FOR 2.5 years. ITS DONE. OVER. I WON. I know there was a bridge to numbe field sieve. And now I got all the details. Its all in the paper now. I knew I hadn't lost my mind! I knew it! Fuck you all assholes. You people all suck. I destroyed everyone in infosec. HAHAHAHHAA. FUCK YOU PUSSIES. ENJOY THE PATCHING ASSHOLES!!!!!!!!!! ECC IS NEXT!

Update: That chapter VII is a mess. I'm going to fix it tomorrow while I also write the final chapter now and the PoC code. I understand now how it is supposed to work. It's weird, because none of this shit that I learned today is new, its all stuff I knew already. Somehow it just lined up in my head the right way I guess. Having these polynomials ommit the + N at the end, it makes it so that we output results that are basicaly the modulus - N. You want that different between modulus and polynomial output for it to work. Because than you can go hunting for quadratic residues. And then you get number field sieve... but better. ERGH. TOMORROW. WATCH. IT ALL LINED UP IN MY HEAD NOW.

Update: Going to delete chapter VII. Its just convoluted. I understand it now. Ill write it down tomorrow and begin work on the PoC... both will be done before the end of the week. God damnit. It was so obvious.. I must suck at math. How did I miss this for so many months. Cursing right now. I should quit math. None of the greats like Galois or Legendre or Fermat would have struggled this long. They would have figured it out in a day. FUCK. What a waste of time. Tomorrow I will fix everything. i see it now. FUuuuck. Im fucking mentally challenged. Everyone must think im an idiot. Ofcourse I cqnt fucking focus bc you people ruined my life! FUCK YOU. Tomorrow this ends. THIS ENDS TOMORROW

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



