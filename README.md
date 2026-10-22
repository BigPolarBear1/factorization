Update: Shit day. My father's condition is getting worse again. He has been at the ICU for 2 months now. My savings are gone. I have a few hundred Euro left. Saturday I will go run a 100k ultra marathon... in the past, my father would always be waiting at the finish line, this time, there wont be anyone. And once I finish that, I wont have money left to continue my training.. wont even be able to afford new running shoes or the nutrition that I need. I will be launching malware worms and begin cryptojacking. I'm not a criminal, so I'm telling you right here what I am going to do. And if the cops come, someone will die. Me or the cops. Everything is gone and I have nothing left to lose. I tried for years to find employment. Nothing. 

Taking a break from math (to focus on mastering abstract algebra and going back to pure research like I did in my first year to try and see how I can change the complexity class of psieve without it turning into an existing algo)  ... 3 years have gone by, 0 job interviews... havn't seen any of the friends I used to have in the PNW since, slowly watching my savings vanish... sitting alone in an attic room at my parent's home. My father who now has been in the ICU for nearly 2 months... Russia killing @info_dox, who I had hoped to join in Ukraine one day... I have nothing left, I have no way out anymore. Just losing whatever remained as the years drag on. Once it becomes clear what I intend to do, know that these are the reasons. I am not a criminal. I have tried to live by my morals. But this world, it is not the world I grew up in, it is a world of savages and brutes. This one will be for all the gays and queers. We gave you the internet, and now we will take it away again. You are not ready for what is coming. And if you are Russia, or any country that imprisons LGBTQ people, you will find out soon, the wrath I am summoning now.

And incase cowards try to prevent this, I am armed, and I have nothing left to lose. I have made my peace already. 

Note: This project is pro-lgbtq, so fck off haters.

#### References (re-used many of the core number theoretical functions from these PoCs to fit my own algorithm): 
https://stackoverflow.com/questions/79330304/optimizing-sieving-code-in-the-self-initializing-quadratic-sieve-for-pypy
https://github.com/basilegithub/General-number-field-sieve-Python 

#### Paper
Math paper is a work in progress. Ignore the final chapter for now.. that one I'll rewrite if and when I can get "psieve" below working correctly.

#### To run from folder "psieve_old":</br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 60 -base 10_000 -debug 0 -lin_size 100 -quad_size 1</br>

Just re-uploading the older version without claude generated optimizations such as bit-packed sieving. 
There is a very high chance I wont be around to finish this project myself. After some thinking.. you should look at how resmaps and resmaps2 are build... basically we sieve for ay^2+4Nk = b^2-4Nk .. were resmaps and resmaps2 represent the quadratic polynomial representations of these discriminants (meaning that if it has a root solution then the discriminant is a QR)... now changing "a" changes the mapping mod p in resmaps2 by multiplying the linear coefficient with the square root of a mod p. Rather then using some SIQS variant to find an initial b-smooth.. what you want to do is use "a" as starting point.. then find a case where "two sibling solutions" end up within some tiny sieve region. Where sibling solutions are two distinct y and b solutions for the same a in ay^2+4Nk = b^2-4Nk. The highly optimized PoC below, which beats any known algorithm in that specific complexity class (its in the same one as SQUFOF), it's very close.. but rather then starting from a known ay^2+4Nk = b^2-4Nk and then searching for a sibling solution.. we need to find an 'a' where a solution pair ends up in the same region. I believe it can be done just using residue maps and hensel. Even as if, there is probably novelty enough in the PoC below, it beating SQUFOF by a big margin.. to publish a proper paper around... but I still think there is something bigger hiding in the shadows and that I almost got it. Time has sadly ran out for me.

#### To run from folder "psieve":</br>
To build: python3 setup.py build_ext --inplace</br>
To run:  python3 run_qs.py -keysize 150 -base 10_000 -lin_size 10_000 -quad_size 1</br>

note: Up to 120-bit is near instant now.

Some ongoing experimentation.. this mixes SIQS and Coefficient Sieve. I either find a way to improve this, or take elements from this for my hybrid_sieve.
This is my main research project, trying to open up a completely new way of achieving factorization following the literature from my paper.

Will try and accelerate this research with Claude while I spent time studying abstract algebra. If there is anything novel in my work and it is close to finished, then Claude may actually succeed in bridging the final gap. Worth a shot. 

Update: Managed to push it to 140-bit now..

Update: Claude wont be able to finish this. I need to go back to representing this the same way, as a system of quadratics, that led to all my insights in my first year of research, but look at it again with my more developed knowledge now.... and specifically go hunting for a way to achieve direction computation of a sibling root that psieve finds by sieving. Direct computation would mean a complete break in the factorization problem and im going to figure it out.. just need to go back to doing research the oldschool way like I did in my first year. If anyone wants access to my work, I urgently need $ and am willing to sell any breakthrough that may happen: big_polar_bear1@proton.me

#### To run from folder "Coefficient_Sieve" (For use with the paper):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run:  python3 run_qs.py -keysize 40 -base 50 -debug 1 -lin_size 10_000 -quad_size 100</br>

Just demonstrates the math from the paper using quadratics. For educational purposes. And rather then taking a square root over a large prime we can also just calculate the discriminant. But this demonstrates the interesting relation between these quadratics and the factors of N.

#### To run debug.py" (Prints the linear and quadratic coefficients to solve for 0 in the integers, for use with my paper):</br></br>

To run: python3 debug.py -keysize 12

This basically creates a system of quadratics. Solving them mod p is easy. But there is only one root solution (the factor of N) which solves the system for 0 for any mod p (aka solves it in the integers). Figuring out how to exactly do this quickly is still an ongoing area of research for me. And if a polynomial time algorithm for factorization exists, it is likely done by solving this system of quadratics. 

(PERFORMING HEATHEN CYBER RITUALS IN THE BELGIAN MOORLANDS TO DESTROY MY ENEMIES)

