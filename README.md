Note: I am antifa, leader of the cipherpunk, fuck the FBI department. 

#### To run from folder "QS_variant" (Standard SIQS with our number theory as backend):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 200 -base 6000 -debug 1 -lin_size 100_000 -quad_size 1</br></br>
This version is somewhat done, but it doesn't achieve much of an advantage, as the main strength of my work will be the NFS_Variant, which performs NFS with reducible quadratic polynomials<br><br>
#### To run from folder "NFS_variant" (Number Field Sieve with our number theory as backend using reducible quadratic polynomials mod m):</br></br>
python3 polarbearalg_debug.py -key 4387 -z 1 (this code just demonstrates the math, use the paper to understand.. full algorithm will be released in the coming days)</br></br>
The files in this folder, relate to my work of my own (better) variant on number field sieve.</br></br>
polarbearalg_debug.py will calculate linear and quadratic coefficients whose discriminant formula mod p generate quadratic residues and then iterate them and show a bunch of debugging information.
I will finish this algorithm in the coming day. But it showcases the bridge to Number Field Sieve already. In addition see the paper. 

Polarbearalg_debug.py also has highly minimized code compared to some of the earlier code I shared. After 2.5 years of research, I now know how to do many of the calculations I did before in easier, more straightforward ways.

Update: I finished a first complete draft of the paper. The entire full algorithm is described there now. I also uploaded two pocs to NFS_Variant, polarbearalg_debug and polarbearalg_debug2, read the last chapter of the paper to understand them.
I will now work on finishing the full algorithm in code, which will be uploaded shortly...

Note: This absolutely destroys the factorization problem. And either a massive deception is ongoing and people know already or you are all about to find out.... either way, good luck I guess.

Update: I'm suddenly super anxious. I see how to do it now... it's in the paper too. People would have known. Someone math educated with experience in factorization would have known the moment they saw my work. Then why? There has to be some deception going on. I'm not sure what it means. They may be moving against me to try and deny me the credits of my work or something. Write me out of history like they always try to do. I will begin writing the final algorithm tomorrow. I know how to do it now. That shouldn't take more then a week to completely finish. 

For tomorrow I want to atleast the very least do this:

Do what polarbearalg_debug2 is doing. But we should lift all those primes to even exponents and only focus on those coefficients that generate 0 mod p for y1 (like the QS_Variant does) ... because we need those solutions to build up what is in NFS the rational side. And then perform trial factorization on what is generated at the algebraic side. Then I can upload that already... and then focus on getting the linear algebra working in the days after, aswell as using jacobi symbols in that step to get better odds of getting a square in the integers for the algebraic side.
