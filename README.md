Note: I am antifa, leader of the cipherpunk, fuck the FBI department. 

#### To run from folder "QS_variant" (Standard SIQS with our number theory as backend):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 200 -base 6000 -debug 1 -lin_size 100_000 -quad_size 1</br></br>
This version is somewhat done, but it doesn't achieve much of an advantage, as the main strength of my work will be the NFS_Variant, which performs NFS with reducible quadratic polynomials<br><br>
#### To run from folder "NFS_variant" (Number Field Sieve with our number theory as backend using reducible quadratic polynomials mod m):</br></br>
To run: python3 polarbearalg_v04 -key 4387 (The debug versions are just for the math in the paper in the final chapter.. I'll rewrite that final chapter one last time once the code is fully done)</br></br>
The files in this folder, relate to my work of my own (better) variant on number field sieve... this is a work in progress</br></br>

I have also uploaded polarbearalg_v05. v04 works fully in the integers on both side of the congruence. v05 converts the bigger side to an algebraic one. However, the linear algebra step in v05 is failing now, because it also needs a quadratic character base (i.e you just take the coefficient y0 and calculate jacobi((y0^2-n\*4\*z),prime). That should force square relations with some multiple of N inbetween.... if you don't do that, the linear algebra step will just yield two squares smaller then N.... which we don't want. Additionally, you may also have to take a square root over a finite field then... but I'll have to check how to do that once my quadratic character base is implemented, I have some ideas how that would work though. 

UPDATE: I AM A MORRON. OVERCOMPLICATING EVERYTHING OMG. GOT IT NOW. ILL UPLOAD THE CORRECT WAY TO DO IT TODAY OR TOMORROW.

UPDATE2: That earlier work I did, with the find_similar stuff, last month. Finding a square multiple of a smooth and then adjusting the coefficients.
That's what I needed to do. But I need to combine that approach with what I'm doing now, using linear algebra. Damn. Let me try to get a working PoC asap.
 
UPDATE3: GOD DAMNIT. I'M TO SLOW IN MY HEAD. FUCK. THIS SHOULDN"T HAVE TAKEN AN ENTIRE MONTH TO CONNECT THE DOTS HERE. MY HEROES GALOIS. FERMAT, LEGENDE AND EULER ARE LAUGHING AT ME FROM THE AFTERLIFE. FUCK. Give me a few days to completely refactor everything and get the correct fucking appraoch working this time (the one where I combine what I was doing a month ago with what I'm doing now)

UPDATE: Almost now. I see now how it can be done. I may be able to finish my work finally either tomorrow or early next week. When this is finally over... it better birth into this world the gay future, or I am going to be seriously angry. Going to be a waste of 2.5 years otherwise. 

UPDATE: O M G. 
I found it.
I found it.
I found it. 
I found it. 
ERGH. I'm such an idiot. FUCK.
Somewhere, some cryptologist, who knows exactly what I am trying to do, is laughing at me because I couldn't connect these final dots for THIS FUCKING LONG.
God damnit. 

Example:

If we calculate linear coeficients mod 15015 where N=4387.

One linear coefficient and one quadratic coefficient will be: y0 = 13502 and z = 1.

If we calculate the discriminant formula mod 15015: 13502^2-4387\*4\*1 = 4,356 mod 15015
Boom. So now, instead of finding a square mod N, we can find a square mod M. 

And further more, calculating the polynomials will also reveal squares (74^2 and 33^2).
SO BASICALLY... WE NEVER NEED TO DO ANYTHING MOD N. WE HAVE CONSTRUCTED MATH SUCH THAT ALL OUR CALCULATIONS MOD N CAN BE DONE MOD M INSTEAD. HELLO GAYS.HHHHHHHHHHHHHHHELLLLOOOOOOO????!!?!!!?!!!!!?. HELLO NSA CRYPTOLOGISTS. It's easy to finish this now. Give me a few days. Fucking hell. I AM CURSING MYSELF RIGHT NOW. I AM A FUCKING IDIOT. WHY DID IT TAKE ME 2.5 YEARS? FUCK. Good thing Fermat is long gone and all my math heroes are dead. Because if any of them would see this embarassment, I would fucking kill myself.

