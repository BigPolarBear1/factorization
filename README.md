Note: I am antifa, leader of the cipherpunk, fuck the FBI department. 

#### To run from folder "QS_variant" (Standard SIQS with our number theory as backend):</br></br>
To build: python3 setup.py build_ext --inplace  (From inside PoC_Files_Find_Similar)</br>
To run: python3 run_qs.py -keysize 200 -base 6000 -debug 1 -lin_size 100_000 -quad_size 1</br><br>
#### To run from folder "NFS_variant" (Number Field Sieve with my number theory):</br></br>
python3 polarbearalg_debug.py -key 4387 -z 1 </br></br>
The files in this folder, relate to my work of my own (better) variant on number field sieve.</br></br>
polarbearalg_debug.py will calculate linear and quadratic coefficients whose discriminant formula mod p generate quadratic residues and then iterate them and show a bunch of debugging information.
I will finish this algorithm in the coming day. But it showcases the bridge to Number Field Sieve already. In addition see the paper. 

Polarbearalg_debug.py also has highly minimized code compared to some of the earlier code I shared. After 2.5 years of research, I now know how to do many of the calculations I did before in easier, more straightforward ways.

Update: Oops, I get it now. I see the entire picture now lol. OMG. It's not even difficult. 

So zx^2-yx will generate one side of the congruence, zx^2-yx+kN will generate another. Where kN is a multiple of N. One of them needs to be a square in the integers... then the other one you just take the square root mod m. Ok ok. Let me fix everything tomorow. Holy fuck shit, I am legit mentally challenged. Fucking taking 2.5 years to make sense of fucking quadratics. Wtf is wrong with me. This is why I dropped out of highschool.. I'm just not cut out for this shit. Just don't have that mental sharpness required to do math... Oh well... still broke factorization I guess. 
