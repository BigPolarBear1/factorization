Disclaimer: At no point did AI contribute anything to this research project. Copilot can't even calculate basic congruences without making mistakes. People who think AI can do novel math research are delusional. I did not use AI for the paper, and I did not use AI for any of the code. I attacked this problem, from the very basics 3 years ago, slowly burrowing deeper and deeper. 

#### To run from folder "NFS_WIP" and "NFS_Variant_Simple" (Will be adding more NFS related code in the coming days/weeks):</br></br>

Note: Experimental WORK IN PROGRESS.</br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 50 -base 500 -debug 1 -lin_size 100_000 -quad_size 1</br></br>

You can run the above command from either two. NFS_Variant_Simple is an intermediate step between QS and NFS and NFS_WIP is merging all my findings into a proper NFS algorithm...

!!!NFS related code is borrowed from: https://github.com/basilegithub/General-number-field-sieve-Python (note: Very impressively written, helped me big time, thanks)

Update: So I'm porting my findings from NFS_Variant_Simple into NFS_WIP. I figured out the math today (february 27th). And began writing some code.. I hope a first rough draft is ready this weekend.
So the insight that I got from NFS_Variant_Simple is that if we perform trial factorization on something like this a-N = poly_val ... as long as poly_val factorizes we can use it. Because a can contain many different roots and coefficients of the shape zx^2+yx. And at the very minimum 1^2+(a-1)\*1 will be a valid solution. Any other roots can be found by factorizing 'a'. We control the size and factorization of a .. hence we can sieve for many different coefficients at the same time, while still being deliberate about the range of coefficients we are sieving for. Functionally this would be the same as running many instances of NFS using different polynomials in parallel for the same effort of only one instance.  This is definately a breakthrough.  I did it... I actually finally did it. I just got to bring it home now.

Update: Finished a first initial proof of concept.. but now I will need to spent time improving it. I'll upload somewhere next week for sure. 

Update: And as an extension of the above.. I am fairly confident it is even possible to multiply sieving results to a common coefficient... such that we are not restricted by the coefficients that can be generated from 'a' (as in a-N = poly_val). That would be ideal to achieve.. as that would mean I could probably outperform existing NFS algorithms using quadratics. 

Update: hehehehe... yea ofcourse you can multiply everything to a common coefficient easily. Its very straightforward. Alright... this is going to set the internet on fire very soon... 

Update: I hope to finish my work next week. So right now I already have a working PoC (not yet uploaded) that just sieves a-N = poly_val and from 'a' deduces all possible coefficients and roots and saves them and then tries to run NFS on these results. Ofcourse this works only if we end up with a bunch of root solutions for the same coefficients. Ultimately what I need to do is multiply all these results so they all share a common coefficient. I did the math and I do believe it can be done... so tomorrow... I will do it.. once that works I should be ready to begin publishing my work for real. And if that doesn't work we can always just select 'a' (as in a-N = poly_val) such that we sieve for the same grouping of coefficients, at the very least it should allow us to sieve for multiply coefficients in one go.. that I know is possible for sure and is already novelty enough.. but I want to destroy factorization.

Update: Ofcourse microsoft deleted my linkedin account. Too much of a coward to face what they have done. Just like those cowards who threatened me with a gun. Guess cowardice is a thing about people living in Redmond.

Update: Yea... can just change coefficients by multiiplying poly_val with a square   ... this also changes 'a' in a-N = poly_val and we get a different set of possible coefficients and roots. Then its just a matter of multiplying a bunch of different smooths until they all share the same coefficients and plug that into NFS... and done. And that's how you make history by sieving with arbitrary polynomials for the number field sieve algorithm. Ofcourse everyone is too much of a coward to admit I've made a breakthrough. I guess the goal for the coming days is to now write a function to multiply smooth values we found by sieving a-N = poly_val with squares, such that the new 'a' value for all smooths can be decomposed into the same coefficients. And that will conclude my work.. because sieving a-N = poly_val is infinitely easier then sieving a^2 - N = poly_val. Literally turned a quadratic problem into a linear problem. The irony is... if people know.. and this is how I am being treated, history isnt going to be kind on you cowards. Either people know, and they are telling people to be quiet for "national security reasons" or people are being willfully ignorant of what is happening. And I garantuee you, had it been any other (non-transgender) person who had achieved this.. they atleast would have had an income and oppurtunities.. they wouldn't have been treated this way. Good luck explaining your way out of that. There was absolutely no reasons for all the despair of the last few years. I'm sure as hell never going to forgive any of you people. Never.

Update: Its a really bad day today. I always used to have a couple of irl friends and  handful of online friends who would regularly reach out. I cant see my irl friends anymore bc they live on the other side of the world and I'm broke. And this last year.. nobody really talk to me anymore. Just slowly the messages stopped coming as months went by. I'm thinking a lot about suicide lately. A lot. It is becoming very difficult to have any hope left. I'm trying to finish this math, but what if it doesn't end up working the way I'm thinking it will work? What then? What will be left in my life? Then the last glimmer of hope will be gone. 

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

