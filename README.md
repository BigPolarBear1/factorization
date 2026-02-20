Disclaimer: At no point did AI contribute anything to this research project. Copilot can't even calculate basic congruences without making mistakes. People who think AI can do novel math research are delusional. I did not use AI for the paper, and I did not use AI for any of the code. I attacked this problem, from the very basics 3 years ago, slowly burrowing deeper and deeper. 

#### To run from folder "NFS_Variant" (Will be adding more NFS related code in the coming days/weeks):</br></br>

Note: Experimental WORK IN PROGRESS.</br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 14 -base 30 -debug 1 -lin_size 100 -quad_size 10</br></br>

NFS related code is borrowed from: https://github.com/basilegithub/General-number-field-sieve-Python (note: Very impressively written, helped me big time, thanks)

This code is a work in progress. I'm trying to merge some of my findings from https://github.com/BigPolarBear1/factorization/tree/main/NFS_Variant_simple into this one.

To run from NFS_Variant_simple, which is an intermediate step between QS and NFS: python3 run_qs.py -keysize 50 -base 300 -debug 1 -lin_size 100_000 -quad_size 10

Update: Just run NFS_Variant using python3 run_qs.py -keysize 14 -base 30 -debug 1 -lin_size 100 -quad_size 10, it may or may not succeed, if it fails generate another number. I basically just removed the linear algebra step for now. 
So we could just use the hashmap and lift all the solutions to even powers to find these cases. But the easier idea would just be to add the linear algebra step again, but rather then restricting ourselves to sieving only multiples of a single polynomial.. I want to be able to sieve much more then that.. now to do that when we multiply polynomial values and zx+y together... we must reconstruct a polynomoial for the product and use that to take a square root over a finite field.. it shouldn't be too hard, because I have worked out all that number theory. Plus it works already with NFS_Variant_Simple, so I dont see why I wouldn't be able to pull it off with this setup as well. Let me begin tomorrow, hopefully this wont take more then a few days. 

Update: Crippling depression. Until my code is done, there will always be that doubt  "what if I'm wrong" .. I struggle with dark thoughts every day, and lately they have been becoming very intense. If I cant innovate on the NFS algorithm.. its over, my life will literally be over. And eveything I experienced. The things I experienced in the US, getting harassed, threatened with a gun outside my apartment, they way microsoft retaliated against my manager for defending me while portraying me as incompetent.. it deeply changed me as a person. I dont think I'll ever be the same again. Seen so much ugliness in this shit world. Its really becoming insanely difficult to keep pushing forward. Anyway.. I started properly investigating what I found earlier. There is cases where we can combine different quadratics, which are not multiples of eachother with linear algebra. Cases where the symbols in the quadratic character base dont get affected after multiplying. I'm still trying to figure out the exact details, I think I somewhat got it though. If I can pull that off, thats going to be great. Just really depressed man. I enjoy the days I'm angry and pissed, because it sure beats the days when I'm feeling like this.

I havnt seen my friends in years bc of what microsoft did. I'm going to kill myself if this doesnt work now. I'm tried man. Getting nowhere in life. Being broke. Just nothing. Its been years since i've last seen any friends. Nobody even reaches out anymore. Just wasting away the years. Just no respect for anything I have done in my career. Just pissed all over everything I have ever done. Microsoft even want back to downgrade some of my bugs (even though I had black on white PoCs proving otherwise) to I would assume retroactively justify them portraying me as incompetent. This is the world we live in. What a nightmare.

Update: I spent all day looking at this math. Let me begin working on what I now I can achieve for sure. 

To do:

1. Index a hashmap by zx+y, zx+y must be square, return by the hashmap will be the possible z,x,y values.
2. All solutions in the hashmap should be lifted to even exponents. This way we can quicly calculate if the polynomial value will also be square
3. Step 1 and 2 allows us to really quickly pull relations from the hashmap where zx+y is square and the polynomial value is square. Then the last thing will be the using the rest of the factor base as a quadratic character base.

Now this I can code in one day easily. Let me just start with that. 

Update: Depression. It is so many years ago now. I think, rollerskating with my teamlead around the seawall in Vancouver, or drinking beers on the beach. And all the other people I got to know there.. that was the peak of my life. All feels so distant now. I can barely even remember the faces of all these people. My life is over. I'm finished. I tried man. I dont care what they will say. What those people at microsoft will say. late 2023 was the worst time in my life. I was threatened with a gun, I got shouted at, called a "p*do" while outside running just minding my own business. In addition to stress of having to deal with the fact I had to leave all my friends behind and come back to a country where I have felt stuck and alone my entire life. I dont think a lot of people would have handled that with a lot of grace. Anyway, its over, bye.

#### To run from folder "CUDA_QS_variant":</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run:  python3 run_qs.py -keysize 240 -base 100_000 -debug 1 -lin_size 100_000_000 -quad_size 100</br></br>
 
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

