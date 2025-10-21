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

Also, now that I have a PoC where I can find a square relation that is correct after the linear algebra step and ones that arn't correct.. I can debug what is going on, calculate the full polynomials and roots and all that... so if it is at all possible... I will know within the coming days. Everything else is set-up correctly already. Just this one last thing now, and I know I have figured out enough of the number theory by now that I can achieve it..
