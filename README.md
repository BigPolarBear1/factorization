Note: I am antifa, leader of the cipherpunk, fuck the FBI department. 

#### To run from folder "QS_variant" (Standard SIQS with our number theory as backend):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 200 -base 6000 -debug 1 -lin_size 100_000 -quad_size 1</br></br>
This version is somewhat done, but it doesn't achieve much of an advantage, as the main strength of my work will be the NFS_Variant, which performs NFS with reducible quadratic polynomials<br><br>
#### To run from folder "NFS_variant" (Number Field Sieve with our number theory as backend using reducible quadratic polynomials mod m):</br></br>
To run: python3 polarbearalg_debug.py -key 4387 

I just uploaded polarbearalg_07.py
To run: python3 polarbearalg_07.py -key 4387

This one finds smooths (or squares) mod m. This works if the modulus is large enough. 
HOWEVER... we can work with smaller moduli aswell... using all the number theory we have figured out so far... and THAT, I will upload very soon...

And also due to working mod m, we can also take square roots mod m etc.... get that same deal going like number field sieve...


Update: So polarbearalg_07 ... it works for large enough moduli because of having a small quadratic coefficient relative to the modulus. Meaning that if we have one square, and we subtract n\*4\*z, it doesn't end up wrapping around the modulus. But..... and here is the interesting part. If we do have a wrap around... then even after performing the linear algebra step... we'll likely end up with two squares that generate a different value mod N.... but it just becomes a matter of "how often do I need to add the modulus such that they become congruent mod N?". But as long as that number (the one we need to reach by adding the modulus to our square mod m) ... is also a square, which we should be able to garantue using legendre symbols... then we can simply take a square root mod m and we don't need to know what that other huge number will be... as long as both of them are the same square (or just a square residue I suppose?) mod m. I see how it works. I'll have to rewrite that final chapter in my paper once more... but let me finish the code first now... otherwise I keep going round and round. I don't understand why I wasn't able to comprehend this earlier... because now that I see it so clearly, it's kind of simple. I guess the last time I tried to do this, I hadn't yet fully understood quadratic coefficients.... and that greatly helped to understand what is happening I guess.
