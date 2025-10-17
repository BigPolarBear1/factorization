Note: I am antifa, leader of the cipherpunk, fuck the FBI department. 

#### To run from folder "QS_variant" (Standard SIQS with our number theory as backend):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 200 -base 6000 -debug 1 -lin_size 100_000 -quad_size 1</br></br>
This version is somewhat done, but it doesn't achieve much of an advantage, as the main strength of my work will be the NFS_Variant, which performs NFS with reducible quadratic polynomials<br><br>
#### To run from folder "NFS_variant" (Number Field Sieve with our number theory as backend using reducible quadratic polynomials mod m):</br></br>
To run: python3 polarbearalg_v04 -key 4387 (The debug versions are just for the math in the paper in the final chapter.. I'll rewrite that final chapter one last time once the code is fully done)</br></br>
The files in this folder, relate to my work of my own (better) variant on number field sieve... this is a work in progress</br></br>

I have also uploaded polarbearalg_v05. v04 works fully in the integers on both side of the congruence. v05 converts the bigger side to an algebraic one. However, the linear algebra step in v05 is failing now, because it also needs a quadratic character base (i.e you just take the coefficient y0 and calculate jacobi((y0^2-n\*4\*z),prime). That should force square relations with some multiple of N inbetween.... if you don't do that, the linear algebra step will just yield two squares smaller then N.... which we don't want. Additionally, you may also have to take a square root over a finite field then... but I'll have to check how to do that once my quadratic character base is implemented, I have some ideas how that would work though. 

UPDATE: EUREKA!!!! I ENDED UP LOOKING SOME MORE AT THE MATH TODAY. I GOT IT. YOU NEED TO USE THOSE 0 SOLUTIONS (like the QS_Variant does), NOT ALL QUADRATIC RESIDUES. AND IF YOU USE THOSE 0 SOLUTIONS.... I NOTICED SOMETHING INTERESTING :) :) :) :) :) :). FRIENDS, TOMORROW, WE'RE GOING TO PERFORM MAGIC. TOMORROW, WILL BE A GREAT DAY! 
