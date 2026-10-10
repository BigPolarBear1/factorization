Note: This project is pro-lgbtq, so fck off haters.

Disclaimer: This is a research project I have worked on for 3.5 years. I have not used AI for any of the research in the paper. Trying to merge an SIQS algorithm with NFS is something I have tried for a long time and my own findings indicated that it should be possible. After feeding in code of a nearly finished project and aggressively prompting Claude with specific directions on how to finish it, I did manage to finally achieve this. But none of the research math was done by Claude, I layed down those foundations myself.

#### References (re-used many of the core number theoretical functions from these PoCs to fit my own algorithm): 
https://stackoverflow.com/questions/79330304/optimizing-sieving-code-in-the-self-initializing-quadratic-sieve-for-pypy
https://github.com/basilegithub/General-number-field-sieve-Python 

#### About the paper
Math paper is a work in progress. Ignore the final chapter for now.. that one I'll rewrite if and when I can get "psieve" below working correctly.

#### To run from folder "psieve":</br>
To build: python3 setup.py build_ext --inplace</br>
To run:python3 run_qs.py -keysize 200 -base 12000 -lin_size 100_000 -quad_size 1 -mode nfs -nfs_degree 4 -nfs_base 6000 -nfs_sing_small 80_000 -nfs_lines 100 -nfs_want 1000 -nfs_lp 26</br></br>

Requirements: You should probably pip install gmpy2, not required but will help big num performance

This is a hybrid SIQS / NFS algorithm. Both feed the same matrix but serve a different purpose. SIQS finds b-smooths using a very large factor base while NFS tries to optimize the matrix. This idea is an ongoing area of research I have been working on for years now. I recently discovered that Claude has finally matured enough to rapidly prototype these ideas for me... so that's what I'm doing as I can now do in a day what would otherwise take me weeks of manual coding and labour. I am still skeptical about AI as a research tool, since it lacks creativity.. but for implementing documented things, even complex math code.. it has definitely matured enough now and I'm becoming a convert.

I'll do some more research myself now. I'll also fix the -mode psieve again.. I need to study that to figure out how to optimize the NFS part.... those two modes are in idea somewhat related.. what still has to be done now is making the NFS algorithm better at improving the overall matrix shared between both algorithms. 

#### To run from folder "Coefficient_Sieve" (For use with the paper):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run:  python3 run_qs.py -keysize 40 -base 50 -debug 1 -lin_size 10_000 -quad_size 100</br>

Just demonstrates the math from the paper using quadratics. For educational purposes. And rather then taking a square root over a large prime we can also just calculate the discriminant. But this demonstrates the interesting relation between these quadratics and the factors of N.

#### To run debug.py" (Prints the linear and quadratic coefficients to solve for 0 in the integers, for use with my paper):</br></br>

To run: python3 debug.py -keysize 12

This basically creates a system of quadratics. Solving them mod p is easy. But there is only one root solution (the factor of N) which solves the system for 0 for any mod p (aka solves it in the integers). Figuring out how to exactly do this quickly is still an ongoing area of research for me. And if a polynomial time algorithm for factorization exists, it is likely done by solving this system of quadratics. 

(PERFORMING HEATHEN CYBER RITUALS IN THE BELGIAN MOORLANDS TO DESTROY MY ENEMIES)

