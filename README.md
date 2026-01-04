One final time before I do what I must, looking for work, big_polar_bear1@proton.me. Cant work for big tech bc of what happened in 2023. Cant travel to the US anymore bc I would get arrested (hence cant work for US companies either). Cant say I didnt try when all the infosec shitheads start pointing fingers again.

Disclaimer: At no point did AI contribute anything to this research project. Copilot can't even calculate basic congruences without making mistakes. People who think AI can do novel math research are delusional.

#### To run debug.py" (Prints the linear and quadratic coefficients to solve for 0 in the integers, for use with my paper):</br></br>

To run: python3 debug.py -keysize 12

This basically creates a system of quadratics. Solving them mod p is easy. But there is only one root solution (the factor of N) which solves the system for 0 for any mod p (aka solves it in the integers). Figuring out how to exactly do this quickly is still an ongoing area of research for me. And if a polynomial time algorithm for factorization exists, it is likely done by solving this system of quadratics. Finding a polyomial time algorithm is my ultimate goal, as this would make progress toward solving p = np as well. 

#### To run from folder "CUDA_QS_variant" (WIP):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run:  python3 run_qs.py -keysize 180 -base 4000 -debug 1 -lin_size 10_000_000 -quad_size 50</br></br>

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

Update: Bah, woke up in the middle of the night. So did some more work. I re-added p-adic lifting. Now when it tries to find similar smooth candidates, it uses p-adic lifting and only sieves even exponents if prime is larger then dupe_max_prime. Additionally thresvar2 sets the threshold tighter then when we do our normal sieving. This has the effect that any smooths we find wont have any new large factors. Hence they have a great chance of forcing the linear algebra step to succeed early. Anyway.. I need some more sleep. I'll fine tune the parameters tomorrow and optimize everything. Feeling very optimistic... I know I tried this before.. but somehow, it works much better this time around.. I think its because Im just worrying about the large factors now... which is what really matters anyway.

You know... seeing it finally behave like it should.. aside from finetuning parameters, optimizing and finishing the paper.. Im basically done. After that I'm taking my research private.. and start the hunt for a polynomial time algorithm.... anyway... glad atleast this chapter is about to come to an end now...

Update: Didn't do much work due to insomnia today. I'll continue tomorrow. So what we need to do now is, finetune those parameters... and also determine in advance if a quadratic coefficient is even worth sieving when we try to find smooth candidates with similar factorization... preferably we pass the test with legendre symbols for many small primes.. that would dramatically increase the odds of finding what we need.

Update: I'll try to finish this when I can.. it's basically finished already.. needs finetuning and optimization work now.. all the math is there though. Just struggling with really severe depression and dark thoughts. This world has lost its color a long time ago.. its all so god damn ugly. And like I also dont see what is going on. I'm not stupid. I know I'm right with my math, and I know people must know. And the only reason they woudln't tell me is that they must be freaked out and that it opens up attacks on RSA-1024. Which is funny... because I had been telling microsoft, the fbi, everyone I was working on this before I published anything at all.. and now that I did publish out of desperation and not seeing another way forward anymore, this happened. A year of silence, unemployment, hopelessness, isolation... while western leaders are mocking people like me. I just want to be alone in the arctic mountains right now. I want to witness the beauty of the natural world, the empty frozen landscapes... because atleast there the world still has some color left.

-------------------------------------------------------------------------
NOTE: Starting 2026, none of my research will be published. Only people who treat me with respect will be allowed access to my work. And NATO countries/big tech are very low on that list after harassing me for years and treating me like shit. And I garantuee you, I will succeed at finding a polynomial time algorithm. There is no one else alive in this fucking world more determined then me to succeed at this. Fucking losers.

ps: If I one day succeed at finding a polynomial time algorithm, I will make sure it is used to defeat america. Fuck you all. Nazi losers. Little man pete hegseth. Pathetic piece of shit. And fuck europe too, all these european leaders love that american nazi shit. Fuck them all. Spineless cowards. All of this could have had a different outcome, yet you people decide to be fucking nazis and harass people like me. Pay the price shitheads.

pps: 2 years ago. I desperately wanted to sell my work. People who could have made that happen knew I was working on factorization. Yet I was treated like shit. And after I released the first iteration of my work, you could have still negotiated with me. Yet you people chose the adversarial route. And this is 100% due to politics. Make no mistake. People will tell plenty of lies to justify what happened. These western nazis just really hate transgender people, so much it destroys their common sense. They will bend the truth. This could have so easily gone another way. In the beginning, I wouldn't even have asked for much. Its far too late. The amount of suffering this last year I had to endure, despite knowing I was right about my math.. the silence.. being in this dead end, broke, unemployed, hopeless... not seeing my friends anymore because of what microsoft did, etc. I hate the west. This wasn't always the case, and in the past when I said it, it was more of a pendantic joke. But now I truely hate the west, I so deeply hate the west now, nothing will ever change that again. You people have shown your disgusting nazi colors. I will do what I can to make sure the west loses. You people chose this route by treating me like shit. You could have just given me employment 2 years ago, and back then, I would have happily shared everything and kept stuff private. But I guess I'm trans, so I must be DEI huh? Fuck you dumbasses. Fuck you all. YOU PEOPLE MADE THIS HAPPEN. NO MATTER WHAT LIES YOU WILL COME UP WITH. THIS IS THE TRUTH.
