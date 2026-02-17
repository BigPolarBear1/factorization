Disclaimer: At no point did AI contribute anything to this research project. Copilot can't even calculate basic congruences without making mistakes. People who think AI can do novel math research are delusional. I did not use AI for the paper, and I did not use AI for any of the code. I attacked this problem, from the very basics 3 years ago, slowly burrowing deeper and deeper. 

#### To run from folder "NFS_Variant" (Will be adding more NFS related code in the coming days/weeks):</br></br>

Note: Experimental WORK IN PROGRESS.</br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 14 -base 100 -debug 1 -lin_size 100 -quad_size 1</br></br>

NFS related code is borrowed from: https://github.com/basilegithub/General-number-field-sieve-Python (note: Very impressively written, helped me big time, thanks)

This code is a work in progress. I'm trying to merge some of my findings from https://github.com/BigPolarBear1/factorization/tree/main/NFS_Variant_simple into this one.

To run from NFS_Variant_simple, which is an intermediate step between QS and NFS: python3 run_qs.py -keysize 50 -base 300 -debug 1 -lin_size 100_000 -quad_size 10

Update: Just run NFS_Variant using python3 run_qs.py -keysize 14 -base 100 -debug 1 -lin_size 100 -quad_size 1, it may or may not succeed, if it fails generate another number. I basically just removed the linear algebra step for now. 
So we could just use the hashmap and lift all the solutions to even powers to find these cases. But the easier idea would just be to add the linear algebra step again, but rather then restricting ourselves to sieving only multiples of a single polynomial.. I want to be able to sieve much more then that.. now to do that when we multiply polynomial values and zx+y together... we must reconstruct a polynomoial for the product and use that to take a square root over a finite field.. it shouldn't be too hard, because I have worked out all that number theory. Plus it works already with NFS_Variant_Simple, so I dont see why I wouldn't be able to pull it off with this setup as well. Let me begin tomorrow, hopefully this wont take more then a few days. 

I guess what I'll do tomorrow is this: If we find a smooth, just change the linear coefficient and try to find the same smooth... so that the product is square.. then modify create_solution() so it works. Once that works, its just a matter of readding the linear algebra step etc. So right now, its already working by succeeding with a single smooth.. next step we make it work with 2 smooths... I'm reducing it to the simplest cases because that's the easiest way to make sure all the number theory gets implemented correctly. Now if I can get it to work with two smooths having distinct coefficients already tomorrow.. then things are going to start escalating really quickly.. because that alone would proof everything. 

I'm struggling big time with depression. Haunted by the past. I wish I could go back in time, to that day when I first arrived in Vancouver Canada, and relive those few good years, with all the friends I had there. The contrast between those years and what my life has become.. some days that alone makes me want to k*ll myself. And also what microsoft did to my former manager, I think that is the hardest thing to cope with. The only person in this industry to believe in me and support me. He even came all the way to Belgium in 2019 to come be my friend while this entire shit industry was harassing me for dropping 0days. Honestly, if I ever come face to face with the people who retaliated against him because he defended me, there would be extreme violence. I'm not a violent person, but I honestly just want to kill those fuckers for what they did. Just treating people like shit. Good people getting screwed over again and again by corporate snakes and their stupid politics. I cant ever have a corporate job again after what happened.. bc I would just lose my shit over all the hypocritical corporate bullshit. MSFT used to keep rambling about integrity and trust and they even had thise shitty tv show everyone was supposed to watch, preaching their corporate values.. its so bullshit, its the most rotten company in the world. Microsoft is a company of fucking snakes. Pathetic snakes.

Btw, I know stalkers would be reporting this repo over and over again, everytime I say shit like "fuck america' or write violent shit when I'm fucking trauma dumping about the past and in a shit mood.. becuase that is what these stalkers had been doing for years. Yet, this repo has never been taken offline. Is it because people know I would just stop sharing my work? huh? Maybe I should have just kept it to myself, sell it to the Chinese when finished.. fuck this entire situation man, i hate everything. I hate belgium most of all. Fucking piece of shit country. worst country in the world. I never felt as isolated and hopeless then the years i've lived in belgium. I dont know what I'm going to do. I might just legit kill myself soon. Fuck all of this. 

Anyway, incase all my suspicions are correct, then I know you will see this: Hello nsa fuckers. Fuck you shitheads. Quit your jobs. And the FBI too. Punch you in your stupid faces shitheads.

What does it matter. I know this will work, sieving with many different quadratics. I dont care what people think. I know what I got. I dont care that so many people laughed at me in the last few years for doing this. I dont care about the silence for years. The desperation. The inability to find employment. Because I'm at the threshold now. I'm going to make history and everything can go to fucking hell. What a shit world this is anyway. Btw, I piss on western intelligence. You could have avoided this so easily. Yet all you incompetent shitheads losers let it escalate right up to this point. Dont ever cross my path. You will die. Dumbass shitheads. Quit your jobs. Go do something that requires less intelligence, bc you are all clearly lacking in that department. HAHAHA> Fucking losers.

Aaaanyway. Maybe I know everything. Maybe I see everything. Maybe I am an opponent you cannot hope to defeat. But I got to play stupid, because knowing things in ways people cannot comprehend, I think people get scared. And then you just end up even more isolated then before. I miss my former manager a lot. There havnt been a lot of humans in this world who didnt make me feel like an outcast. Who for atleast a short time, made me feel like I belonged somewhere. Until everything came crashing down again. Now there is nothing left. I started learning math in the summer of 2023. I didn't know any math and had dropped out of highschool. In a couple more years, I dont think this will be funny anymore. And I'll be even more miserable and even more isolated. Or dead. Oh well. I would really like to do that months long ski expedition in Scandinavia though. Just sleeping in a tent, as I pull a sled/pulk into the arctic in the middle of winter. I can be free of everything alone in the cold. Nature doesn't notice me that I'm there. Nothing can hurt me in nature.

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

