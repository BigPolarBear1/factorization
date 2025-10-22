Note: I am antifa, leader of the cipherpunk, fuck the FBI department. 

#### To run from folder "QS_variant" (Standard SIQS with our number theory as backend):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 200 -base 6000 -debug 1 -lin_size 100_000 -quad_size 1</br></br>
This version is somewhat done, but it doesn't achieve much of an advantage, as the main strength of my work will be the NFS_Variant, which performs NFS with reducible quadratic polynomials<br><br>
#### To run from folder "NFS_variant" (Number Field Sieve with our number theory as backend using reducible quadratic polynomials mod m):</br></br>
To run: python3 polarbearalg_debug.py -key 4387 

I just uploaded polarbearalg_08.py
To run: python3 polarbearalg_08.py -key 4387

This is slowly getting there now. It will work with p_mod_amount=6 but not with p_mod_amount=5 ... because the issue is that if we make the modulus smaller, and we getting wrapping around when we modulo reduce the discriminant output, even with the quadratic character base... we just end up with a discriminant that generates square residues within that quadratic character base... but not an actual square. But atleast it generates square residues. 
The last ingredient missing in that PoC, is to get that quadratic character base, to generate squares in the integers, not just quadratic residues... and for that, I believe I should be able to pull it off now using the polynomial equations and everything I learned there. Because something very distinct happens there when we have a square relation in the integers... I'm sure I can leverage it. Let me spent a few days doing the math... once that is done.... that's it.. my work will be done.

Some info about polarbearalg_08, during the linear step with p_mod_amount=5 and -key 4387 you would see for example:</br></br>

sqrt_right:  1903</br>
sqrt_left:  1507</br>
disc%N: 2970 disc**0.5: 1.3131543342473657e+24 jacobi_symbols: [0, 1, 0, 1, 1, 1, 1]</br>

Disc%N represents sqrt_right in the integers, where is sqrt_right is just mod m. 
And sqrt_left^2%N will be equal to Disc%N:

1507^2 = 2970 mod 4387

And ofcourse Disc mod m will be equal to sqrt_right mod m.

The challenge that remains is that disc**0.5 (square root of the discriminant) should yield an integer (aka, be a square in the integers). And right now, it's aready a square residue within the quadratic character base.. but calculating legendre symbols using the discriminant isn't enough to force a square in the integers there... hence I believe this is where the polynomials now come into play. I'll update as soon as I figure out this final step.

Also, now that I have a PoC where I can find a square relation that is correct after the linear algebra step and ones that arn't correct.. I can debug what is going on, calculate the full polynomials and roots and8 all that... so if it is at all possible... I will know within the coming days. Everything else is set-up correctly already. Just this one last thing now, and I know I have figured out enough of the number theory by now that I can achieve it..

Update: One important requirement also seems to be that the gcd between your modulus and quadratic coefficient is 1. There may be a way to get it to work when that is not the case... but its easier to just skip those quadratic coefficients entirely. Aaaanyway... I'll figure this out soon. 

Update: I'll have to fix the get_root() function, it works flawlessly when the quadratic coefficient is 1. But for all other quadratic coefficients, its not grabbing the roots I actually want. Let me adjust that tomorrow... because if we have a square relation, then the correct root should yield a quadratic residue with the polynomial, mod any prime.... and that is the key also to getting some type of quadratic character base working. Getting those correct roots, that's really important... I'll call it a day for now. Tomorrow I'll fix that root getting code...it's this combination of discriminant and polynomials that's going to tell me wether or not I have a square in the integers or just a square residue. I'm getting closer every day... just got to keep grinding... it's taking way too long. The years are flying by. Just sitting alone in this room. Doing math. No income. Never seeing any friends. People think they can break me, force me into obscurity. But I am at my strongest when the only other way out is death. I break factorization or I die. There is no alternatives. I have no other choices. 

Update: Bah, I couldnt sleep so I looked at whats going on with those roots. That get_root() function, it works fine with the QS_Variant PoC... aka when the jacobi/legendre symbol is 0. But with quadratic residues we need to use tonelli. But anyway... if a linear coefficient is paired with the correct quadratic coefficient, then such a root can be found for any mod p... and if not, then it wont find a solution for all mod p. I'll work out the math tomorrow how to incorporate that with the linear algebra step. With all that I know now, there has to be a way to setup it up so it results in a square in the integers. 

Update: It's funny, that a lot of people don't understand the true consequences of breaking PKI ciphers. A lot of the 0day crowd even don't understand it, but then again, they are too high on their 0-click blah blah blah, whatever shit. Breaking PKI is everything. All encrypted traffic first establishes a secure channel to share symmetric keys. And if you break PKI, you break the secure channel that is used to share symmetric keys, hence you can decrypt 99.99% of encrypted traffic. Let alone the nightmare it would cause to authentication protocols. If one were to succeed, there is no bug more severe. No 0day would even come close to the power this would grant a person. Not even remotely. Maybe the closest would have been bugs like curveball (which also, the 0day crowd didn't understand the true implications of). I still do not understand why the decision was made to burn a bug like that.. I can only speculate that they must have had something better. And if they do, I am going to blind them, I am going to take it whatever it is they have, away from them, for everything these fuckers did. I am fighting the west. Waging war against the west. Don't misunderstand my intentions. You people started this war.

Oh btw "Referring sites and popular content are temporarily unavailable or may not display accurately. We're actively working to resolve the issue" ... what would be the odds this is broken because you people are in a fight for your lives trying to not give away the fact you are all panicking and scared shitless? Fuck you shitheads. I'm waging war against the west. And I'm fucking destroying you people. Fucking losers. And fuck those who serve them. They are all the enemy.

Btw dont ever get the wrong idea, if you work for western infosec, I consider you to be the enemy. I'm not your friend. Fuck you all fuckers. And if you are an nsa cryptologist seeing their secrets going up in flame, hahahahaha, I piss on you losers. and fuck pete hegseth. Watch your back you dumb fucker. Shit happens. 

There is no greater danger in this world then america. They have proving that by everything I experienced while in the US and by the way their leaders dehumanize people like me non-stop. I will fight america until my last dying breath. Anything I can do to undermine america, I will. I'm not a coward and I'm not afraid. Fuck america. Fuck all those nazi pieces of shit. 
