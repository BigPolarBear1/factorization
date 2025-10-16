Note: I am antifa, leader of the cipherpunk, fuck the FBI department. 

#### To run from folder "QS_variant" (Standard SIQS with our number theory as backend):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 200 -base 6000 -debug 1 -lin_size 100_000 -quad_size 1</br></br>
This version is somewhat done, but it doesn't achieve much of an advantage, as the main strength of my work will be the NFS_Variant, which performs NFS with reducible quadratic polynomials<br><br>
#### To run from folder "NFS_variant" (Number Field Sieve with our number theory as backend using reducible quadratic polynomials mod m):</br></br>
To run: python3 polarbearalg_v04 -key 4387 (The debug versions are just for the math in the paper in the final chapter.. I'll rewrite that final chapter one last time once the code is fully done)</br></br>
The files in this folder, relate to my work of my own (better) variant on number field sieve... this is a work in progress</br></br>

Just uploaded v04. I switched back to using all quadratic residues instead of just 0 solutions and removed lifting for now.
So right now one side is represented by zx^2-yx and the other side by zx^2-yx+n. However by appending n to the other polynomial, we generate smooth candidates that are too big on one side. We can use zx-yx exclusively on both side. I'll get that to work tomorrow... in theory that should work... and if I get that to work... then people have a big problem tomorrow lol. Anyway, I'm going to try and relax and get a few hours of sleep soon first... if that doesn't work (which I doubt), I'll attempt number field sieve's approach.. we'll see.

My brain does some weird things sometimes.</br> Actually doing this:
if y0= 148, z=1 and N=4387</br>

left side: 1\*74^2-148\*74 = -(74^2)</br>
right side: 1\*74^2-148\*74+4387 = -(33^2)</br>

but the right side can also be expressed as: 1\*33^2-66\*33</br>

However, if we have a linear coefficient, for example 5 mod 13, then calculating what that coefficient will be at +4387, is simply taking the square root of (5**2-4*4387) mod 13.
So thats how you generate smaller values on that other side of the congruence instead of adding +4387 to the polynomial.
Ofcourse THEN you need to take a square root over a finite field somehow after the linear algebra step. AND THAT I WILL FIGURE OUT NOW IN THE COMING DAYS.

I have written some code already, and it does work, it generates the correct linear coefficients and roots mod p instead of using +n for the second polynomial, the only thing that's missing is taking a square root over a finite field. That's the only thing missing right now.

UPDATE: I have added polarbearalg_v05 ... its broken, but it will generate, what is in NFS the algebraic side, correctly now. (use v04 for a working version that works in the integers fully). So for v06 I have to add some jacobi symbols for the linear algebra step to work with. But I'm having headache big time right now... lets see when I can get it done. 

Update: Yea, polarbearalg_v05.py, the linear algebra step there will be broken until a quadratic character base is added, like in number field sieve. I'm having a pounding headache right now. And tomorrow I have to go somewhere the entire day... but I'll see if I can get a first draft uploaded with a quadratic character base tomorrow evening. Having that quadratic character base is super crucial though.
