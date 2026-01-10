One final time before I do what I must, looking for work, big_polar_bear1@proton.me. Cant work for big tech bc of what happened in 2023. Cant travel to the US anymore bc I would get arrested (hence cant work for US companies either). Cant say I didnt try when all the infosec shitheads start pointing fingers again.

Disclaimer: At no point did AI contribute anything to this research project. Copilot can't even calculate basic congruences without making mistakes. People who think AI can do novel math research are delusional.

#### To run debug.py" (Prints the linear and quadratic coefficients to solve for 0 in the integers, for use with my paper):</br></br>

To run: python3 debug.py -keysize 12

This basically creates a system of quadratics. Solving them mod p is easy. But there is only one root solution (the factor of N) which solves the system for 0 for any mod p (aka solves it in the integers). Figuring out how to exactly do this quickly is still an ongoing area of research for me. And if a polynomial time algorithm for factorization exists, it is likely done by solving this system of quadratics. Finding a polyomial time algorithm is my ultimate goal, as this would make progress toward solving p = np as well. 

#### To run from folder "CUDA_QS_variant" (WIP):</br></br>
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

Update: While making this PoC, I suddenly had an idea, linking back to some things I tried earlier last spring. See below (I think this approach is interesting, but its hard to get an advantage.. however using these different multiples of N.. there is a much easier way to exploit this, hence see NFS_variant):

#### To run from folder "nfs_variant" (WIP):</br></br>

Note: Experimental WORK IN PROGRESS.</br>
To build: python3 setup.py build_ext --inplace</br>
To run:  python3 run_qs.py -keysize 26 -base 50 -debug 1 -lin_size 1000 -quad_size 1</br></br>

This uses the number theory from the paper and full quadratics. It will succeed if two distinct roots for the same coefficient produce the same results mod N. 
I'm still trying to figure out how to sieve this somehow.

Update: Nearly there now. When we check the hashmap, we should get that to work for higher exponents as well. Anyway, the more important part is quickly finding that other root which produces the same polynomial value.. and now that I have a way to know what the residues are garantueed to be for that other root, that k variable (zx^2+yx-nk) and the discriminant.. I actually may be able to use p-adic lifting now to quickly find it... I should try to write a function to do this tomorrow.. I'll begin that first thing tomorrow.. bc if I can get that to work, then thats a simple and elegant solution to very quickly find that other root once we have a smooth.. hmm.

note to self: I should also check that instead of finding two roots that produce the same polynomial value... can I do the same with 2 different roots, whose polynomial values multiplied together are equal to a third root's polynomial value... and so on.. because if it still holds true in that case... then we can just use linear algebra easily.. I should check this first thing tomorrow... before I check anything else... that surely would be the path of least resistance toward finishing this.

</br></br></br></br></br></br></br></br></br></br></br></br></br></br></br></br></br></br></br></br></br></br>
-------------------------------------------------------------------------
#### Rants
NOTE: Starting 2026, none of my research will be published. Only people who treat me with respect will be allowed access to my work. And NATO countries/big tech are very low on that list after harassing me for years and treating me like shit. And I garantuee you, I will succeed at finding a polynomial time algorithm. There is no one else alive in this fucking world more determined then me to succeed at this. Fucking losers.

ps: If I one day succeed at finding a polynomial time algorithm, I will make sure it is used to defeat america. Fuck you all. Nazi losers. Little man pete hegseth. Pathetic piece of shit. And fuck europe too, all these european leaders love that american nazi shit. Fuck them all. Spineless cowards. All of this could have had a different outcome, yet you people decide to be fucking nazis and harass people like me. Pay the price shitheads.

pps: 2 years ago. I desperately wanted to sell my work. People who could have made that happen knew I was working on factorization. Yet I was treated like shit. And after I released the first iteration of my work, you could have still negotiated with me. Yet you people chose the adversarial route. And this is 100% due to politics. Make no mistake. People will tell plenty of lies to justify what happened. These western nazis just really hate transgender people, so much it destroys their common sense. They will bend the truth. This could have so easily gone another way. In the beginning, I wouldn't even have asked for much. Its far too late. The amount of suffering this last year I had to endure, despite knowing I was right about my math.. the silence.. being in this dead end, broke, unemployed, hopeless... not seeing my friends anymore because of what microsoft did, etc. I hate the west. This wasn't always the case, and in the past when I said it, it was more of a pendantic joke. But now I truely hate the west, I so deeply hate the west now, nothing will ever change that again. You people have shown your disgusting nazi colors. I will do what I can to make sure the west loses. You people chose this route by treating me like shit. You could have just given me employment 2 years ago, and back then, I would have happily shared everything and kept stuff private. But I guess I'm trans, so I must be DEI huh? Fuck you dumbasses. Fuck you all. YOU PEOPLE MADE THIS HAPPEN. NO MATTER WHAT LIES YOU WILL COME UP WITH. THIS IS THE TRUTH.
