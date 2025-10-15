Note: I am antifa, leader of the cipherpunk, fuck the FBI department. 

#### To run from folder "QS_variant" (Standard SIQS with our number theory as backend):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 200 -base 6000 -debug 1 -lin_size 100_000 -quad_size 1</br></br>
This version is somewhat done, but it doesn't achieve much of an advantage, as the main strength of my work will be the NFS_Variant, which performs NFS with reducible quadratic polynomials<br><br>
#### To run from folder "NFS_variant" (Number Field Sieve with our number theory as backend using reducible quadratic polynomials mod m):</br></br>
To run: python3 polarbearalg_v04 -key 4387 (The debug versions are just for the math in the paper in the final chapter.. I'll rewrite that final chapter one last time once the code is fully done)</br></br>
The files in this folder, relate to my work of my own (better) variant on number field sieve... this is a work in progress</br></br>

Just uploaded v04. I switched back to using all quadratic residues instead of just 0 solutions and removed lifting for now.
So right now one side is represented by zx^2-yx and the other side by zx^2-yx+n. However by appending n to the other polynomial, we generate smooth candidates that are too big on one side. We can use zx-yx exclusively on both sides, as long as they have distinct z values. I'll get that to work tomorrow... in theory that should work... and if I get that to work... then people have a big problem tomorrow lol. Anyway, I'm going to try and relax and get a few hours of sleep soon first... if that doesn't work (which I doubt), I'll attempt number field sieve's approach.. we'll see.

God damnit. Agitated. So damn agitated. If this works.. if it is as simple as iterating quadratic coefficients for the same linear coefficient.. then tomorrow I am breaking factorization so hard, it wont even be funny. You start to cope with long term isolation and trauma in really strange ways. Sometimes I want to believe that this has some cosmic importance, that the entire universe is nudging me on. Because the events that have driven me to dedicate my life to this are so surreal. It does feel like someone has been interfering with my life. But that train of taught isn't very healthy. But what if... and what does that mean to the future I am helping create? It better be the really gay and queer future or I'm going to be so mad. Anyway.. potentially a big day tomorrow... we'll see... if it doesn't work.. plenty of other approaches left to try. My understanding keeps growing, and for now that is all that matters.
