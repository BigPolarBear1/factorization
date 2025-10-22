Note: I am antifa, leader of the cipherpunk, fuck the FBI department. 

#### To run from folder "QS_variant" (Standard SIQS with our number theory as backend):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 200 -base 6000 -debug 1 -lin_size 100_000 -quad_size 1</br></br>
This version is somewhat done, but it doesn't achieve much of an advantage, as the main strength of my work will be the NFS_Variant, which performs NFS with reducible quadratic polynomials<br><br>
#### To run from folder "NFS_variant" (Number Field Sieve with our number theory as backend using reducible quadratic polynomials mod m):</br></br>
To run: python3 polarbearalg_debug.py -key 4387 

I just uploaded polarbearalg_09.py</br>
To run: python3 polarbearalg_09.py -key 4387  (or -keysize 14 to generate a random 14 bit modulus)</br>

I have fixed the quadratic character base! Atleast with p_mod_amount=6 (which is still a very large modulus compared to the N when N=4387) .... but I have seen cases where we have a "wrap around" of the modulus, and still get a square in the integers now. I'm able to generate plenty such examples now (test with -keysize 14). Now I need to reduce the modulus.. but with atleast the quadratic character base working, I assume the next step is getting some type of square root over finite field thing working now. That should be the very last thing. And it should definitely be do-able, even without a polynomial ring like NFS. I'm fairly confident I have figured out enough of the number theory now to get it working. When that is done... that's it. My vengence will have been completed. I will have destroyed the american cryptologic advantage. For everything they have done. 

Update: Oh yea, I just realized.. if a wrap around happens, the PoC is still using the square it found mod m.. but you should use the square it finds from the full discriminant for taking the GCD.. anyway... I'm going to sleep for now... I'll upload v10 tomorrow and make some progress toward reducing that modulus and using NFS's approach of taking a square root. My hope is that the quadratic character base now forces it to the correct quadratic coefficient atleast when using very small moduli.. else I will need to figure out how to build up a big modulus with the linear algebra step... which is what NFS does... so there is two different approaches I'll need to experiment with. 

You know, I was just thinking about western infosec and the whole tech-bro exploit dev scene. It's really a disgusting industry with disgusting people. They really hate people like me. I ever meet a western exploit dev, I'm going to punch them. Losers. I hate spooks too. I ever see a spook irl, I will do worse things then just punch. I'm not afraid of people like that. Fuck them. I'm about to take on the entire fucking world. I'm well aware of what's happening. And I'm ready. I've nothing left to lose. I'm ready to fight the entire fucking world. Bring it on fucking losers.
