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

My brain does some weird things sometimes.</br> Actually doing this:
if y0= 148, z=1 and N=4387</br>

left side: 1\*74^2-148\*74 = -(74^2)</br>
right side: 1\*74^2-148\*74+4387 = -(33^2)</br>

but the right side can also be expressed as: 1\*33^2-66\*33</br>

However, if we have a linear coefficient, for example 5 mod 13, then calculating what that coefficient will be at +4387, is simply taking the square root of (5**2-4*4387) mod 13.
So thats how you generate smaller values on that other side of the congruence instead of adding +4387 to the polynomial.
Ofcourse THEN you need to take a square root over a finite field somehow after the linear algebra step. AND THAT I WILL FIGURE OUT NOW IN THE COMING DAYS. 
