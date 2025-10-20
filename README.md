Note: I am antifa, leader of the cipherpunk, fuck the FBI department. 

#### To run from folder "QS_variant" (Standard SIQS with our number theory as backend):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 200 -base 6000 -debug 1 -lin_size 100_000 -quad_size 1</br></br>
This version is somewhat done, but it doesn't achieve much of an advantage, as the main strength of my work will be the NFS_Variant, which performs NFS with reducible quadratic polynomials<br><br>
#### To run from folder "NFS_variant" (Number Field Sieve with our number theory as backend using reducible quadratic polynomials mod m):</br></br>
To run: python3 polarbearalg_debug.py -key 4387 

I just uploaded polarbearalg_07.py
To run: python3 polarbearalg_07.py -key 4387

This one finds smooths (or squares) mod m. This works if the modulus is large enough. 
HOWEVER... we can work with smaller moduli aswell... using all the number theory we have figured out so far... and THAT, I will upload very soon...

And also due to working mod m, we can also take square roots mod m etc.... get that same deal going like number field sieve...


Update: I added a quadratic character base to polarbearalg_07 and it appears to be working now... I'll check it some more tomorrow just to be absolutely sure... because there is always a chance my math is producing false positives due to testing on small integers (for now). But if this was all I had to do.... easy. Damnit. I hope the actual solution will be a little bit more complicated... because if this is it... god damnit... I remember experimenting with this appraoch back in april, but I couldn't get it to work back then.. ofcourse I know a little bit more then I knew back then... but still.... god damn. The problem with math research is, I managed to get very far into this, with no math education, but then fully going from smooth finding in the integers to mod m... that's when things become much more abstract and difficult. And if you don't do everything exactly right, it straight up won't work, and it's incredibly hard to debug or figure out why it isn't working. I should have seen all of this half a year ago... I'm giving people too much time to react to my work. I should just break factorization, be merciless, and get revenge for what microsoft did to my manager and for seperating me from all my friends. The best friends I ever had in my life. I never had many friends, but my former managers and all the people I got to know in Vancouver... the best years of my life. And I can't go back to it. All there is left is sadness and anger now. I have to succeed. And then I must attack diophantine equations and ECC. Maybe I wasted a lot of time figuring out basic shit while attacking factorization, but it's all part of the learning curve. Next project will be faster, because I will know how to appraoch it now.
