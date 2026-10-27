Update: Trying some stuff.. and added a few sentences to the final chapter of the paper. Have an idea... let me see first if it works then ill finish fhe paper and code...

Update: Bah close... very close... let me look at the roots of the quadratic representations of these discriminants.. there should be a shared root or something for these sibling solutions since they can also be re-ordered as a difference of squares.... that might help me drive an algorithm to pull these cases from a residue map.... fck it man, let me try for one more week.

Note: This project is pro-lgbtq, so fck off haters.

#### References (re-used many of the core number theoretical functions from these PoCs to fit my own algorithm): 
https://stackoverflow.com/questions/79330304/optimizing-sieving-code-in-the-self-initializing-quadratic-sieve-for-pypy
https://github.com/basilegithub/General-number-field-sieve-Python 

#### Paper
Math paper is a work in progress.
Note: This paper was written at different stages in my  research project, and as such I will need to go back to refactor some of the older chapter as I now have a better understanding of the entire structure.

#### To run from folder "fu_algo":</br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 25 -base 5000 -debug 0 -lin_size 100 -quad_size 1

Pronounced the fuck you algorithm

This will be my final iteration.. reworking psieve so it works purely with residues and direct computation of solutions inside a modulus.
I've added some code to the bottom of process_interval2d(), which replaces the psieve() function. Just showcases the math I recently added to the paper. I need to re-org that code a little bit now though. 

So there is these two conditions that I talk about in the paper. First is two sibling roots both yielding the same value for x^2-4Nk = s ... If we find a b-smooth with SIQS.. on the right hand side we can factor out any squares to get the form x^2-4Nk = s, then we need to find a second sibling root. We can calculate residues inside a modulus. Next, knowing that both roots are also a difference of square itself (the one fermat's factorization algorithm looks for, lets start calling these fermat squares) and thus there must exist some "k" where they have a shared root.. we can then use those residues to figure out how often we need to add the modulus.. this last part I still need to implement but all the math is already in the code.. let me see if I can get it done soon.. very tired after running a 100k ultra marathon last weekend. I know I'm right. I know I almost got it. I'm not going to give up on this until I get it to work, even if it takes the rest of my life. And I know people know this, but they are all just a bunch of fcking cowards and treating me this way is very deliberate. And if I ever meet one of the people behind this, I will become extremely violent, so you better fucking start hiding already you fcking cowards.

I've often skimmed very close to this solution in the last year... guess it just took a little more mathematical growth.. only have been doing math for 3.5 years starting from nothing as a highschool dropout.

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

