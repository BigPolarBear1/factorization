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

Update: So zx^2-yx will generate one side of the congruence, zx^2-yx+kN will generate another. Where kN is a multiple of N. One of them needs to be a square in the integers... then the other one you just take the square root mod m. Ok ok. Let me fix everything tomorow. Holy fuck shit, I am legit mentally challenged. Fucking taking 2.5 years to make sense of fucking quadratics. Wtf is wrong with me. This is why I dropped out of highschool.. I'm just not cut out for this shit. Just don't have that mental sharpness required to do math... Oh well... still broke factorization I guess. 

Update: Updated both the PoC and Paper to demonstrate this... now all we do is complete this by using the smallest polynomial to find a square in the integer, and the other one with +kN to find one in mod m.... and down goes the factorization problem... good luck assholes.

To do:

-Show an example where we have a polynomial with a square in the integers (zx^2-yx) and another polynomial (zx^2-yx+kN) with a quadratic residue mod m, and take the square root over a finite field. I should demonstrate that in code and in the paper

-Demonstrate how to combine polynomials to generate that square relation we just talked about.. just like NFS does... use some linear algebra... easy enough. I'll write this in code first and release it.. then demonstrate a full numerical example in the paper also. 

Update: Oops, that PoC that I uploaded for the NFS_Variant last night was broken. I fixed it now.
