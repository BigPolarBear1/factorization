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
