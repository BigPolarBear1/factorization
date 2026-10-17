Taking a break from math unless I can push claude to complete psieve (still a few things I need to try, there may exist a sub-exponential way of sieving it... will try soon) ... 3 years have gone by, 0 job interviews... havn't seen any of the friends I used to have in the PNW since, slowly watching my savings vanish... sitting alone in an attic room at my parent's home. My father who now has been in the ICU for nearly 2 months... Russia killing @info_dox, who I had hoped to join in Ukraine one day... I have nothing left, I have no way out anymore. Just losing whatever remained as the years drag on. Once it becomes clear what I intend to do, know that these are the reasons. I am not a criminal. I have tried to live by my morals. But this world, it is not the world I grew up in, it is a world of savages and brutes. This one will be for all the gays and queers. We gave you the internet, and now we will take it away again. I promise you this, as my final act, when the time comes. You are not ready for what is coming. And if you are Russia, or any country that imprisons LGBTQ people, you will find out soon, the wrath I am summoning now.

And incase cowards try to prevent this, I am armed, and I have nothing left to lose. I have made my peace already. 


Note: This project is pro-lgbtq, so fck off haters.

#### References (re-used many of the core number theoretical functions from these PoCs to fit my own algorithm): 
https://stackoverflow.com/questions/79330304/optimizing-sieving-code-in-the-self-initializing-quadratic-sieve-for-pypy
https://github.com/basilegithub/General-number-field-sieve-Python 

#### Paper
Math paper is a work in progress. Ignore the final chapter for now.. that one I'll rewrite if and when I can get "psieve" below working correctly.

#### To run from folder "psieve":(work in progress)</br>
To build: python3 setup.py build_ext --inplace</br>
To run:  python3 run_qs.py -keysize 150 -base 10_000 -lin_size 10_000 -quad_size </br>

note: Up to 120-bit is near instant now.

Some ongoing experimentation.. this mixes SIQS and Coefficient Sieve. I either find a way to improve this, or take elements from this for my hybrid_sieve.
This is my main research project, trying to open up a completely new way of achieving factorization following the literature from my paper.

Will try and accelerate this research with Claude while I spent time studying abstract algebra. If there is anything novel in my work and it is close to finished, then Claude may actually succeed in bridging the final gap. Worth a shot. 

Update: Managed to push it to 140-bit now... 

#### To run from folder "hybrid_sieve":</br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 300 -base 200000 -lin_size 1000000 -quad_size 1 -mode nfs -nfs_degree 4 -nfs_base 60000 -nfs_sing_small 1100000 -nfs_lp 27 -nfs_lp2 24 -nfs_t 200000000 -nfs_lines 4 -la_every 2000</br>

Some experimentation using claude to quickly test and implement ideas for a hybrid implementation. However the only way this will actually work as an hybrid algorithm is if we can make NFS spit out b-smooths with a more targeted factorization. I'm also diving deeper into abstract algebra literature to see if I can figure out a breakthrough on that front.  300-bit takes about an hour right now.

#### To run from folder "Coefficient_Sieve" (For use with the paper):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run:  python3 run_qs.py -keysize 40 -base 50 -debug 1 -lin_size 10_000 -quad_size 100</br>

Just demonstrates the math from the paper using quadratics. For educational purposes. And rather then taking a square root over a large prime we can also just calculate the discriminant. But this demonstrates the interesting relation between these quadratics and the factors of N.

#### To run debug.py" (Prints the linear and quadratic coefficients to solve for 0 in the integers, for use with my paper):</br></br>

To run: python3 debug.py -keysize 12

This basically creates a system of quadratics. Solving them mod p is easy. But there is only one root solution (the factor of N) which solves the system for 0 for any mod p (aka solves it in the integers). Figuring out how to exactly do this quickly is still an ongoing area of research for me. And if a polynomial time algorithm for factorization exists, it is likely done by solving this system of quadratics. 

(PERFORMING HEATHEN CYBER RITUALS IN THE BELGIAN MOORLANDS TO DESTROY MY ENEMIES)

