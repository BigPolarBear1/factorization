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

Update: Uploaded polarbearalg_v02:

Use: python3 polarbearalg_v02 -key 4387

Compared to v01, I have further expanded the sieving logic and it also makes sure both sides of the congruence are smooth now (the same way NFS does).
Next I will implement the linear algebra step... if that doesn't work well enough, I'll have to investigate adding a quadratic character base (legendre symbols) as well. 

My initial goal is just to get it to work, no matter how convoluted or bad.. once it works, I can begin optimizing and simplifying everything. I may have to change some things on how I represent both sides of the congruence to get the linear algebra step to work... we'll see. Let me implement the linear algebra step first, and then I can debug the math if it doesn't work straight out of the box.

Update: Thinking some more about it... I'll have to change somethings eventually. What is in NFS the algebraic side, we definitely want that to be just zx^2+yx (without -N)  (or the other way around? i'll figure it out soon)... just to make sure we get small values there. However, what is in NFS the rational side, that I will need to represent differently then I'm doing right now... but for now, let me get this working first... and then I'll simplify and optimize. Get the easiest path working first, and then figure out the rest. 

Like, it doesn't matter if the first iteration will be some ugly NFS/QS hybrid. All that matters is that I keep inching closer to that final solution. And you can only arrive at that final solution by following the path of least resistance. Anyway... I'll take a break and go for a run. Tomorrow I'll upload the linear algebra step if I can get it to work. And after that we optimize and simplify.
Even if I just start by implementing some ugly QS/NFS hybrid that finds a square relation in the roots instead of the coefficients.. .that is again one step closer to full NFS.. and I'll be able to use the code to further debug and figure out the math.
