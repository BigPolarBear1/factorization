Taking a break from math (to focus on mastering abstract algebra) unless I can push claude to complete psieve (still a few things I need to try, there may exist a sub-exponential way of sieving it... will try soon) ... 3 years have gone by, 0 job interviews... havn't seen any of the friends I used to have in the PNW since, slowly watching my savings vanish... sitting alone in an attic room at my parent's home. My father who now has been in the ICU for nearly 2 months... Russia killing @info_dox, who I had hoped to join in Ukraine one day... I have nothing left, I have no way out anymore. Just losing whatever remained as the years drag on. Once it becomes clear what I intend to do, know that these are the reasons. I am not a criminal. I have tried to live by my morals. But this world, it is not the world I grew up in, it is a world of savages and brutes. This one will be for all the gays and queers. We gave you the internet, and now we will take it away again. You are not ready for what is coming. And if you are Russia, or any country that imprisons LGBTQ people, you will find out soon, the wrath I am summoning now.

And incase cowards try to prevent this, I am armed, and I have nothing left to lose. I have made my peace already. 

Note: This project is pro-lgbtq, so fck off haters.

#### References (re-used many of the core number theoretical functions from these PoCs to fit my own algorithm): 
https://stackoverflow.com/questions/79330304/optimizing-sieving-code-in-the-self-initializing-quadratic-sieve-for-pypy
https://github.com/basilegithub/General-number-field-sieve-Python 

#### Paper
Math paper is a work in progress. Ignore the final chapter for now.. that one I'll rewrite if and when I can get "psieve" below working correctly.

#### To run from folder "FUSieve":(work in progress)</br>
To build: python3 setup.py build_ext --inplace</br>
To run:  python3 run_qs.py -keysize 150 -base 2000 -lin_size 200_000 -quad_size 1 </br>

This is the fuck you sieve. A work in progress. Will improve shortly... needs some scoring for leftovers etc...

Update: one sec, I need to link psieve's QR approach with SIQS's divisibility approach.. back in 2025 when I released that first quadratic sieve variant (march of 2025)... I was calculating linear coefficient residues... and then calculating how many times I had to multiply it by N or something to reach 0 mod p to get "k".... so the QR form was transformed to the divisibility one. I did not quite understand that link back then.. the only reason I did that transformation back then was because I wanted to link my work to Quadratic Sieve and that was one way I found how to do it through experimentation... but now I understand the deeper link between quadratic <-> discriminant (aka QR symbol <-> Divisibility). let me defer from using claude for a second. I need to calculate this by hand first. Claude can optimize code, but not do novel math like this.

Update: Eureka! I can use psieve's bitsieve approach... to sieve divisibility in the direction of N multiplier k........... yeap.. that's how you merge the quadratic form, representing the QR symbol of the discriminant and the divisibility of the discriminant into one elegant algorithm...  let me begin writing that final chapter in the paper today. Then I'll feed that into claude tomorrow (waiting for token reset) and have it implemented into this PoC. And that will finally have been the end of that.... after 3.5 years and losing everything. I'm never forgiving this. It's a damn shame too, that I was forced into the path of disclosure like this too. Fucking idiocracy. 

#### To run from folder "psieve":</br>
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

