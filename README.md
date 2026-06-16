Disclaimer: At no point did AI contribute anything to this research project. Copilot can't even calculate basic congruences without making mistakes. People who think AI can do novel math research are delusional. I did not use AI for the paper, and I did not use AI for any of the code. I attacked this problem, from the very basics 3 years ago, slowly burrowing deeper and deeper. I did experiment from time to time with AI, but it always failed to be useful at math. It is useful for writing simple and short python functions for simple things.. but that is the extend of its usefulness. Saying this here because I am closing in fast now and I know people will claim I just used AI or something bc of my non-math background and being transgender. And I'm honestly tired of the trolls feeding my work into AI and then attacking my work using the garbarge analysis AI makes of it.</br>

Disclaimer2: I never finished highschool. I started working on this 3 years ago, not even knowing basic algebra (I can explain exactly how I can to all the conclusions I came to if someone ever cares to hear about it). This is 3 years. I'm only going to get better. I'm still far away from plateauing skill wise. 

The day I break factorization will be day 0 of the gay future. A better world will be born. A world of human creativity. A world of art, science and spirituality. A world without suffering. I will see that world become reality, no matter the cost. 

Been unemployed for years without income, still looking for work: big_polar_bear1@proton.me ... able to relocate to anywhere in the world except the US.

#### References (re-used many of the core number theoretical functions from these PoCs to fit my own algorithm): 
https://stackoverflow.com/questions/79330304/optimizing-sieving-code-in-the-self-initializing-quadratic-sieve-for-pypy
https://github.com/basilegithub/General-number-field-sieve-Python 

#### About the paper
Math paper is a work in progress. The final chapters are a bit rushed and building an algorithm around p-adic lifting isnt as straight forward as I had assumed. I do think there is an angle there I can exploit, but I'll do some further experimentation first and get a working PoC before I make edits to the paper again.

#### To run from folder "polysieve" WORK IN PROGRES...extremely early version:</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run:  python3 run_qs.py -keysize 70 -base 500 -debug 0 -lin_size 10_000 -quad_size 1

EDIT: Important note to self so I dont forget. Going to degree higher then 2 we get two linear polynomials and some leftover part... basically splitting it in 3 parts rather then 2. I should at a minimum make note of that even if I dont end up using it in an algorithm without number fields (its likely a lot more interesting to leverage with number fields in the picture).

This basically starts with a SIQS variant using my own residue sieving that I figured out (so not encoding the modulus into the leading coefficient) then once a B-smooth is found, it calls into find_same.. where I'm using a setup that begins to blur the lines between QS and NFS. But I'm hitting a wall there as I really need to be able to move up in degree, but it's not practical without implementing number fields. And I need a few weeks to study up on that subject. Once I understand that, I will continue. Will optimize everything and get my paper up to date then take a break to study. If anyone wants to buy future work: big_polar_bear1@proton.me .. just uploading shit to github was never my first choice.

Update: Uploaded some work in progress. Moved all the recent findings into find_same() ... now I need to polish that up and in fact, we can presieve the values of x+offset and gval and abuse some mechanism that changes how to get paired up. Then I should try to get that same thing working for higher degree.. but instead of two linear polynomials.. we have more smallers ones that we can presieve and whose individual parts we can re-arrange. 

I will optimize (uploaded code is crap need to fix a lot of things) what is here and document everything in my paper. My intuition of wanting to move up to higher degree is correct.. but I need to learn about number fields.. I just cant punch past these limitations without number fields. So I will documented what I got so far in my paper.. learn number fields and in a year from now upload the next iteration of my paper and hopefully by then I am able to do what I'm actually trying to do. Ofcourse I would be going into the direction of NFS, but hopefully all the intuition and understanding build with this is going to help to push the envelope there.

Update: mainly been focused on making the paper more presentable this week. Also just now quickly fixed that offset calculation so we end up with a constant that divides by the modulus. Tomorrow I'll start with the pre-sieving of the two linear polynomials... then write all that logic... because that's by large the biggest refactor and the sooner I get it over with the better... 

Update: It is now presieving one of the linear polynomials already, need to implement support for the other one aswell. Then we can gut all these loops in find_same and just pair up the indexes of when both the linear polynomials factorize (the idea being that we can pair up arbitrary indexes, not just when both just happen to be smooth.. important distinction). So tonight, I'll add presieving for the other linear poly and then before I gut all these loops I'll also optimize that fval.. because that one needs to be kept as small as possible after dividing out all the known factors. Idea is to shift all the burden to the linear polynomials bc those we can presieve and pair up... this pairing up does have an offset cost depending on how far the indices for the two linear polynnomials are that gets added to fval, so the smaller fval starts with after dividing out the known factors, the better.

Update: Alright, both linear polynomials now get pre-sieved. nothing too complicated about that. Next step now is getting fval as small as possible after dividing out the known factors because we want to shift all of the factorization burden to the two linear polynomials. And after that I am hoping I can use that to generate near identical b-smooths... because sieving this way gives us quite a bit more freedom.

Update: thinking a bit more about this... so fval doesnt need to be super small... just reasonably small after dividing out the modulus factors. And then the factors of the linear polynomials either need to be small or square. And that's how you do the whole "find_same" strategy.. because I have proven in the past that this reduces the required number of B-smooths. Now finding these cases.. that's like a parameter search. Something I am sure can be optimized. Let me write down the core ideas in the paper today so atleast that is out of the way...

Update:  Today was mostly paper day.. figured out some more of the math... so we can change linear poly pairing either by adding the modulus to our offset value or by adding the modulus to the coefficient... and in addition we can also add a skew with that b parameter... similar to what NFS does.. So I should be good to start implementing some proper code now...

#### To run from folder "Coefficient_Sieve" (For use with the paper):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run:  python3 run_qs.py -keysize 40 -base 50 -debug 1 -lin_size 10_000 -quad_size 100</br>

Just demonstrates the math from the paper using quadratics. For educational purposes. And rather then taking a square root over a large prime we can also just calculate the discriminant. But this demonstrates the interesting relation between these quadratics and the factors of N.

#### To run from folder "CUDA_QS_variant" (Failed Experiment):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run: To run:  python3 run_qs.py -keysize 240 -base 100_000 -debug 1 -lin_size 100_000_000 -quad_size 100</br></br>
 
Prerequisites: </br>
-Python (tested on 3.13)</br>
-Numpy (tested on 1.26.2)</br>
-Sympy</br>
-cupy-cuda13x</br>
-cython</br>
-setuptools</br>
-h5py</br>
(please open an issues here if something doesn't work)</br></br>

Additionally cuda support must be enabled. I did this on wsl2 (easy to setup), since it gets a lot harder to access the GPU on a virtual machine.

This was an attempt at finding smooths with similar factorization using an SIQS variant. By using quadratic coefficients. But it didnt end up working as I had hoped so I abondoned this approach, but perhaps someone will get some use out of it.

#### To run debug.py" (Prints the linear and quadratic coefficients to solve for 0 in the integers, for use with my paper):</br></br>

To run: python3 debug.py -keysize 12

This basically creates a system of quadratics. Solving them mod p is easy. But there is only one root solution (the factor of N) which solves the system for 0 for any mod p (aka solves it in the integers). Figuring out how to exactly do this quickly is still an ongoing area of research for me. And if a polynomial time algorithm for factorization exists, it is likely done by solving this system of quadratics. Finding a polyomial time algorithm is my ultimate goal, as this would make progress toward solving p = np as well. 


