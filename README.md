Note: This project is pro-lgbtq, so fck off haters.

Disclaimer: This is a research project I have worked on for 3.5 years. I have not used AI for any of the research in the paper. Trying to merge an SIQS algorithm with NFS is something I have tried for a long time and my own findings indicated that it should be possible. After feeding in code of a nearly finished project and aggressively prompting Claude with specific directions on how to finish it, I did manage to finally achieve this. But none of the research math was done by Claude, I layed down those foundations myself.

#### References (re-used many of the core number theoretical functions from these PoCs to fit my own algorithm): 
https://stackoverflow.com/questions/79330304/optimizing-sieving-code-in-the-self-initializing-quadratic-sieve-for-pypy
https://github.com/basilegithub/General-number-field-sieve-Python 

#### About the paper
Math paper is a work in progress. Ignore the final chapter for now.. that one I'll rewrite if and when I can get "psieve" below working correctly.

#### To run from folder "psieve":</br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 100 -base 2000 -lin_size 100_000 -quad_size 1 -mode nfs</br></br>

This merges SIQS and NFS into an hybrid algorithm. SIQS finds b-smooths with a square.. the larger the square the better NFS will perform. Then NFS runs and feeds b-smooths back to the SIQS algorithm. And we keep repeating this process... my paper also indicates that these number fields should work on quartics. But this has yet to be implemented. 

I have worked on this research project for 3.5 years, without AI, but after the recent OpenAI math drop, I decided to aggresively push claude to try and finish my project with very specific prompting (prompting it to replace the psieve() function with an nfs implementation and giving advice on how to do it.. such as re-using the square part of an SIQS generated b-smooth). Uploaded version proves this works. Next support for quartics... 

Update: I realized that I messed around with exactly this NFS setup in the past... I'll make some modifications tomorrow. I know what to do now :)
The thing about nfs_launch_sq() is that it needs to spit back out b-smooths that reduce the matrix rank significantly.. similar to what psieve() did... claude doesnt see it, but I experimented with this setup in the past (actually spent weeks messing around with it)... so I can do the edits myself. Will do it tomorrow. 

Update: I've come up with the following: SIQS style sieve with a large factor base. Then run  NFS (of arbitrary degree) with smaller factor base to achieve a rank reduction of the SIQS matrix. And just keep repeating that. So don't see NFS as a factoring algorithm but rather as a b-smooth finding algorithm... 

#### To run from folder "Coefficient_Sieve" (For use with the paper):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run:  python3 run_qs.py -keysize 40 -base 50 -debug 1 -lin_size 10_000 -quad_size 100</br>

Just demonstrates the math from the paper using quadratics. For educational purposes. And rather then taking a square root over a large prime we can also just calculate the discriminant. But this demonstrates the interesting relation between these quadratics and the factors of N.

#### To run debug.py" (Prints the linear and quadratic coefficients to solve for 0 in the integers, for use with my paper):</br></br>

To run: python3 debug.py -keysize 12

This basically creates a system of quadratics. Solving them mod p is easy. But there is only one root solution (the factor of N) which solves the system for 0 for any mod p (aka solves it in the integers). Figuring out how to exactly do this quickly is still an ongoing area of research for me. And if a polynomial time algorithm for factorization exists, it is likely done by solving this system of quadratics. 

(PERFORMING HEATHEN CYBER RITUALS IN THE BELGIAN MOORLANDS TO DESTROY MY ENEMIES)

