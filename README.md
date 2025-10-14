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

Update: I have begun finalizing the algorithm under Polarbearalg_v01.py

This will do number field sieve with my own number theory up until the sieving portion.

What still needs to be done is:

If an algebraic smooth is found, we should construct a matrix row.
Once enough such rows are found, we should then do the linear algebra step.
And once a solution is found, we just take a square root over a finite field.

Should all be finished somewhere this week... 

Update: I'll get the linear algebra working tomorrow... straightforward enough... just copy paste from my QS PoC... then I can verify all the math. Then that's either going to work and I can probably finish the PoC completely very quickly.. or there will still be some hidden snakes in the math that I'll need to figure out (but I'm sure I can manage as I now understand things far better then before)... we'll see tomorrow. Today was a slow day.. just depression. I need to finish this soon. I think come January I'm just going to spent a year in the wilderness by myself... these last two years changed me... I feel deeply disgusted by humanity.
