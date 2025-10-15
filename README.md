Note: I am antifa, leader of the cipherpunk, fuck the FBI department. 

#### To run from folder "QS_variant" (Standard SIQS with our number theory as backend):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 200 -base 6000 -debug 1 -lin_size 100_000 -quad_size 1</br></br>
This version is somewhat done, but it doesn't achieve much of an advantage, as the main strength of my work will be the NFS_Variant, which performs NFS with reducible quadratic polynomials<br><br>
#### To run from folder "NFS_variant" (Number Field Sieve with our number theory as backend using reducible quadratic polynomials mod m):</br></br>
To run: python3 polarbearalg_v03 -key 4387 (The debug versions are just for the math in the paper in the final chapter.. I'll rewrite that final chapter one last time once the code is fully done)</br></br>
The files in this folder, relate to my work of my own (better) variant on number field sieve... this is a work in progress</br></br>

LINEAR ALGEBRA IS WORKING!!!!!!!!! WE DID IT!!!!! However, this is now a NFS-QS hybrid. We now do square finding in the integers on both sides of the congruence. However, one side we should reduce that to mod m and take a square root over a finite field. Let me figure out this final step now :) Almost!

Ok.... the fact that the linear algebra works, that's a big morale boost. Even though it isnt yet true NFS... but it is only one final step away from true NFS (well, true NFS, but better and faster).

Update: Ok I got it. I had some clarity and did some deep thinking. So trial factorization on zx^2+yx mod p is easy, since they will be small values. But zx^2+yx-N is not easy since it will be big values. HOWEVER. I know from what I was doing earlier that -N (or -kN) is interchangeable with the quadratic coefficient. So what I need to do, is enumerate other possible quadratic coefficients mod p, and make that my second polynomial instead of appending -N. BOOM. DONE. HAH. Holy fucking shit, I was so close to solving it early summer... but then I made a detour to SIQS instead. I guess I did learn a lot so it wasn't a waste of time. I'll go for a run soon... I'll probably upload tomorrow, I'll have to write a bit of code again. Lets see how long it takes.
