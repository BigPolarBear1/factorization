DISCLAIMER: AT NO POINT IN MY 2.5 YEARS OF RESEARCH DID I USE AI OR HAS AI IN ANY WAY HELPED WITH MY RESEARCH. I'M PUTTING THIS HERE BECAUSE I KNOW PEOPLE WILL ACCUSE ME OF IT. 

Note: I am antifa, leader of the cipherpunk, fuck the FBI department. 

#### To run from folder "QS_variant" (Standard SIQS with our number theory as backend):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 200 -base 6000 -debug 1 -lin_size 100_000 -quad_size 1</br></br>

See below for an improved way of performing what this PoC does.. I'll delete this Proof of Concept once the PoC for Improved_QS_Variant matures a little more<br><br>
#### To run from folder "NFS_Variant" (Sort of a QS/NFS hybrid that improves on Quadratic Sieve using the number theory I developped):</br></br>

I just uploaded polarbearalg_v44.py</br>
To run: python3 polarbearalg_v44.py -keysize 50</br>

UPDATE: Before I attempt anything complicated. Let me finish a PoC that does the following. 

If we do regular SIQS and we start with generating x^2-n divisible by the modulus. If the smooth is too large. In theory by reducing x and adding some small z, we should be able to reduce the size of the smooth using zx^2-n. 

Like, you know how the sieve interval in standard SIQS, creates a parabole around 0. But using z we can shift where exactly it wraps around 0. We can shift the parabole left or right. That's the secret power of this appraoch. I have to quickly work out the details.. fucking headache all the time. grr. 

Update: Ok I'm going to begin implementing this. First I will take my SIQS variant. Change the discriminant of y0^2-N*4*z to x^2-N ... then after that, we can shift the parabole the sieve interval creates with quadratic coefficients. Yea, that will be FUCKING beautiful. People are so fucked. AND NOBODY BELIEVED IN ME. GO TO HELL. 
