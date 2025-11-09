DISCLAIMER: AT NO POINT IN MY 2.5 YEARS OF RESEARCH DID I USE AI OR HAS AI IN ANY WAY HELPED WITH MY RESEARCH. I'M PUTTING THIS HERE BECAUSE I KNOW PEOPLE WILL ACCUSE ME OF IT. 

Note: I am antifa, leader of the cipherpunk, fuck the FBI department. 

#### To run from folder "QS_variant" (Standard SIQS with our number theory as backend):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 200 -base 6000 -debug 1 -lin_size 100_000 -quad_size 1</br></br>

See below for an improved way of performing what this PoC does.. I'll delete this Proof of Concept once the PoC for Improved_QS_Variant matures a little more<br><br>
#### To run from folder "NFS_Variant" (Sort of a QS/NFS hybrid that improves on Quadratic Sieve using the number theory I developped):</br></br>

I just uploaded polarbearalg_v42.py</br>
To run: python3 polarbearalg_v42.py -keysize 50</br>

UPDATE: I actually managed to replace the trial factorization of zx^2 completely with jacobi symbols. This works, but as the bit size goes up, the amount of required jacobi symbols also goes up. I have some ideas to keep it smaller by just performing trail factorization mod m of zx^2 combined with jacobi symbols.. I will upload soon. 

UPDATE: v41 can do 50 bit using a massive amount of jacobi symbols. On the bright side, this does allow us to generate incredibly small smooth candidates on the right side of the congruence. I will reduce the amount of required jacobi symbols tomorrow... I have an idea how to do it. If we manage to reduce that suffeciently... then we have a superior algorithm which is better then both QS and NFS.

UPDATE: v42 halves the required jacobi symbols amount by also adding as condition that zx^2 mod zx^2-n factors over the factor base. Which is still infinitely better then doing trial factorization over zx^2 in the integers. I will further improve it tomorrow. I have a fairly good idea how to do it. But I need to refactor large portions of the PoC first. Ofcourse everything else about the PoC will need to be rewored eventually also... it being slow is as expected for now.
