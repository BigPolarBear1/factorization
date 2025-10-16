Note: I am antifa, leader of the cipherpunk, fuck the FBI department. 

#### To run from folder "QS_variant" (Standard SIQS with our number theory as backend):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 200 -base 6000 -debug 1 -lin_size 100_000 -quad_size 1</br></br>
This version is somewhat done, but it doesn't achieve much of an advantage, as the main strength of my work will be the NFS_Variant, which performs NFS with reducible quadratic polynomials<br><br>
#### To run from folder "NFS_variant" (Number Field Sieve with our number theory as backend using reducible quadratic polynomials mod m):</br></br>
To run: python3 polarbearalg_v04 -key 4387 (The debug versions are just for the math in the paper in the final chapter.. I'll rewrite that final chapter one last time once the code is fully done)</br></br>
The files in this folder, relate to my work of my own (better) variant on number field sieve... this is a work in progress</br></br>

I have also uploaded polarbearalg_v05. v04 works fully in the integers on both side of the congruence. v05 converts the bigger side to an algebraic one. However, the linear algebra step in v05 is failing now, because it also needs a quadratic character base (i.e you just take the coefficient y0 and calculate jacobi((y0^2-n\*4\*z),prime). That should force square relations with some multiple of N inbetween.... if you don't do that, the linear algebra step will just yield two squares smaller then N.... which we don't want. Additionally, you may also have to take a square root over a finite field then... but I'll have to check how to do that once my quadratic character base is implemented, I have some ideas how that would work though. 

I know this will come to completion very soon now. And when that very final version comes online, I promise you, one day in the future, there will be a reckoning. This is not my fault. This is a disaster of your own making. And either way, at this point, you all deserve this disaster. I couldn't care less anymore. You people have destroyed all that was good in my life. And in the meanwhile, your leaders are openly ridiculing people like me. Yea, you deserve this. I hereby reclaim technology, it belongs to the outcasts again. Because you people don't deserve the promise of technology. You've only brought misery with it into this world.
