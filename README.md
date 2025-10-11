Note: I am antifa, leader of the cipherpunk, fuck the FBI department. 

#### To run from folder "QS_variant" (Standard SIQS with our number theory as backend):</br></br>
To build: python3 setup.py build_ext --inplace  (From inside PoC_Files_Find_Similar)</br>
To run: python3 run_qs.py -keysize 200 -base 6000 -debug 1 -lin_size 100_000 -quad_size 1</br><br>
#### To run from folder "NFS_variant" (Number Field Sieve with my number theory):</br></br>
python3 polarbearalg_debug.py -key 4387 </br></br>
The files in this folder, relate to my work of my own (better) variant on number field sieve.</br></br>
polarbearalg_debug.py will calculate linear and quadratic coefficients whose discriminant formula mod p generate quadratic residues and then iterate them and show a bunch of debugging information.
I will finish this algorithm in the coming day. But it showcases the bridge to Number Field Sieve already. In addition see the paper. 

Polarbearalg_debug.py also has highly minimized code compared to some of the earlier code I shared. After 2.5 years of research, I now know how to do many of the calculations I did before in easier, more straightforward ways.

To do:

1. Add to the paper how to find the correct linear and quadratic coefficient if we find a square mod m with the quadratic polynomial. Since we must adjust it so it generates that square with the polynomial in the integers instead. This the analogue for NFS's taking square roots over finite fields step.

2. Add that to the debug PoC

3. Add to the paper how to find squares mod m, which we can then finish by adjusting the coefficients and taking the gcd.

4. Write the non-debug PoC which combines all these things.

5. Taunt the FBI and everyone who has wronged me in my life. Because I won. I know I won. And I told all my Chinese friends already and they will be destroying your networks soon. Because go to hell for how you people treat me. Fucking losers.
