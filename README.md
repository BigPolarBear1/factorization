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

Update: Actually what I was doing before was better. Using the polynomials with the +N and -N ommitted. You probably want to use those instead. Since they generate smaller values in the integers. It's the same thing like number field sieve... where that entire integer value needs to factor over your algebraic factor base. Because else you can't really take your square root over a finite field. It becomes complicated then. Let me show tomorrow how to do it... slow day today.. my way of doing math research, is very different from what most people do probably... just a lot of python code and print functions, like with polarbearalg_debug.py ... and then I look at the output and make sense of the patterns. Thats how I have been doing math for 2.5 years now. Alright, I need a break for today now.
