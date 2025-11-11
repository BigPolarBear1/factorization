DISCLAIMER: AT NO POINT IN MY 2.5 YEARS OF RESEARCH DID I USE AI OR HAS AI IN ANY WAY HELPED WITH MY RESEARCH. I'M PUTTING THIS HERE BECAUSE I KNOW PEOPLE WILL ACCUSE ME OF IT. 

Note: I am antifa, leader of the cipherpunk, fuck the FBI department. 

#### To run from folder "QS_variant" (Standard SIQS with our number theory as backend):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 200 -base 6000 -debug 1 -lin_size 100_000 -quad_size 1</br></br>

See below for an improved way of performing what this PoC does.. I'll delete this Proof of Concept once the PoC for Improved_QS_Variant matures a little more<br><br>
#### To run from folder "Improved_QS_Variant" (Sort of a QS/NFS hybrid that improves on Quadratic Sieve using the number theory I developped):</br></br>

I just uploaded polarbearalg_v44.py</br>
To run: python3 polarbearalg_v44.py -keysize 50</br>

Or to run my work in progress from this folder which doens't yet demonstrate the new number theory but it swaps coefficients for roots:

To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 200 -base 6000 -debug 1 -lin_size 100_000 -quad_size 1</br></br>

UPDATE: Before I attempt anything complicated. Let me finish a PoC that does the following. 

If we do regular SIQS and we start with generating x^2-n divisible by the modulus. If the smooth is too large. In theory by reducing x and adding some small z, we should be able to reduce the size of the smooth using zx^2-n. 

Like, you know how the sieve interval in standard SIQS, creates a parabola around 0. But using z we can shift where exactly it wraps around 0. We can shift the parabola left or right. That's the secret power of this appraoch. I have to quickly work out the details.. fucking headache all the time. grr. 

Update: Ok I'm going to begin implementing this. First I will take my SIQS variant. Change the discriminant of y0^2-N*4*z to x^2-N ... then after that, we can shift the parabole the sieve interval creates with quadratic coefficients. Yea, that will be FUCKING beautiful. People are so fucked. AND NOBODY BELIEVED IN ME. GO TO HELL. 

Update: First step complete... it now uses roots instead of linear coefficients. I uploaded that PoC already. Next... we're going to shift parabola with quadratic coefficients :)... this may take a few days.. because it requires quite a bit of code and refactoring..

Update: Ok I'm done uploading for now. It's obvious now. If 1*x^2-N produces a parabola when iterating x, then zx^2-N produces a shifted parabola. And polarbearalg_v44 shows how that can be made to work. Dude... people are so fucking screwed. I know exactly what this implies. I don't need NFS when I have this lol.

Update: Alright another day. Lets get some work done. Goal for today is to implement into my SIQS variant the use of quadratic coefficients in the shape of zx^2-N like polarbearalg_v44.py does

