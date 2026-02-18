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

Update: I'll see what I get done today. I'm going really dark places again mentally. I know I've won. I dont even need to modify create_solution. I just need to write some code that will use the hashmap to reconstruct a quadratic when we multiply results together, such that we have a quadratic producing that product, and then we just plug that into create_solution() like I'm doing now. Anyway, tomorrow will be a wasted day again, because the losers of the belgian judicial system are forcing me to do therapy, so I have to run to therapy for 30km bc I have legs that can run and public transit is for losers. I want to leave this shit country so badly. There is just nothing for me here. My entire life in belgium has been hopeless and socially isolated. Belgium is the worst country in the entire world. I think if I could, I would probably accept a job offer in China if it came my way..  better then rotting away in Belgium, wasting my life like this.

Update: I did some testing. It seems to be that we can really just multiply any polynomials together as long as zx+y is the same.. which relates to the factorization of the constant in the discriminant. Which makes sense and is in-line with what I suspected earlier. That gives us still an awful big range of polynomials to sieve with. I will begin finishing my work now and prepare to formally publish toward the end of the month. Please follow me on linkedin where I will announce when I do, thankms (seems fitting to drop it on a platform owned by microsoft after the misery they caused).

ps: I will go dark now until I publish. Fuck the NSA, you dumbass losers. Fuck the FBI, morrons. Fuck microsoft, shitheads. I was harassed and threatened with a gun in the US, fired, havnt seen any of my friends since because you fuckers left me with no other options but to return back to this shity country where I never had any oppurtunities, just endless harassment by the shit judicial system. The only person to ever support me and believe in me in this shit industry was fired after working for 27 years at microsoft, why? Because he defended me. They portrayed me as incompetent and him as someone who defend someone who was incompetent. FUck you all. I have every right to be angry. If I see any of you shitheads, I will become violent. I'm going to Asia after I publish my work. Fuck this. How can I ever stay in the west after these last few years? I'm not a dog, I have my pride. And I know what I got. I know what I have sunk my teeth into. Yet you people dared to gaslight me, treated me like I was psychotic or some shit. Yea, you nearly drove me to madness these last few years, anyone would, fucking shitheads. Dont ever come near me, because bad shit will happen. I am not human, I am polar bear dressed like a human. Fuck you.

update: Got back from running 30k. Too tired now to work but tomorrow I begin finalizing everything, everything just be ready to drop here on github around the end of the month, I know exactly how to implement my number field sieve variant now. I just want to say this last thing. I know people will portray me like a raging angry lunatic. Bc that is this industry's default. Everyone knew I was working on this. I told microsoft and the FBI long before (literally a year before) I even shared with the public for the first time that I was working on factorization. Nobody acted. NOTHING WAS DONE. Do you think I want to be in the middle of this shit storm? I would have much prefered fucking financial security for a while and some fucking peace and quiet. Like, every fucking time this shit industry forces this outcome where I am forced to publish my work. Thanks for destroying my life and thanks for what this shitstorm will do to everyone around me. Go to hell. I am angry bc I understand the consequences. And that is why, I ever see any of you fuckers who led it come this far, I will be violent. Morrons. They should rebrand the intelligence sector to the fucking mentally challenged shithead sector because everyone working in western intelligence is mentally challenged (and I dont mean that as disrespect to actually mentally challenged people, bc atleast those people arn't careless on purpose like this). Fuck you all. Pay my former manager 10 billion for the misery you caused or fuck off and I'll be accepting the first job offer coming out of Asia. I'm only going to get better at this. I'm never stopping. Should atleast have tried to kill me if you wernt going to make a deal. Morrons. Losers. Literal big time fucking losers. Had it coming I guess. Shitheads. Drop dead fuckers.

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

