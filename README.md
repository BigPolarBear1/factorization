Note: I am antifa, leader of the cipherpunk, fuck the FBI department. 

#### To run from folder "QS_variant" (Standard SIQS with our number theory as backend):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 200 -base 6000 -debug 1 -lin_size 100_000 -quad_size 1</br></br>
This version is somewhat done, but it doesn't achieve much of an advantage, as the main strength of my work will be the NFS_Variant, which performs NFS with reducible quadratic polynomials<br><br>
#### To run from folder "NFS_variant" (Number Field Sieve with our number theory as backend using reducible quadratic polynomials mod m):</br></br>
To run: python3 polarbearalg_debug.py -key 4387 

Ignore the other PoCs in that folder. 
I learned some things over the weekend. 
What I was doing last month with square multiple of a smooth. 
I managed to finally understand that much deeper and generalizing it to just having to find a square mod m. 
so if we have y0^2-N\*4\*z instead of finding a square in the integers, we find one mod m. 
I'll try to upload a PoC tomorrow.

