Note: I am antifa, leader of the cipherpunk, fuck the FBI department. 

#### To run from folder "QS_variant" (Standard SIQS with our number theory as backend):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 200 -base 6000 -debug 1 -lin_size 100_000 -quad_size 1</br></br>
This version is somewhat done, but it doesn't achieve much of an advantage, as the main strength of my work will be the NFS_Variant, which performs NFS with reducible quadratic polynomials<br><br>
#### To run from folder "NFS_variant" (Number Field Sieve with our number theory as backend using reducible quadratic polynomials mod m):</br></br>
To run: python3 polarbearalg_debug.py -key 4387 

I just uploaded polarbearalg_09.py
To run: python3 polarbearalg_09.py -key 4387  (or -keysize 14 to generate a random 14 bit modulus)

I have fixed the quadratic character base! Atleast with p_mod_amount=6 (which is still a very large modulus compared to the N when N=4387) .... but I have seen cases where we have a "wrap around" of the modulus, and still get a square in the integers now. I'm able to generate plenty such examples now (test with -keysize 14). Now I need to reduce the modulus.. but with atleast the quadratic character base working, I assume the next step is getting some type of square root over finite field thing working now. That should be the very last thing. And it should definitely be do-able, even without a polynomial ring like NFS. I'm fairly confident I have figured out enough of the number theory now to get it working. When that is done... that's it. My vengence will have been completed. I will have destroyed the american cryptologic advantage. For everything they have done. 
