One final time before I do what I must, looking for work, big_polar_bear1@proton.me. Cant work for big tech bc of what happened in 2023. Cant travel to the US anymore bc I would get arrested (hence cant work for US companies either). Cant say I didnt try when all the infosec shitheads start pointing fingers again.

Disclaimer: At no point did AI contribute anything to this research project. Copilot can't even calculate basic congruences without making mistakes. People who think AI can do novel math research are delusional.

NOTE: Starting 2026, none of my research will be published. Only people who treat me with respect will be allowed access to my work. And NATO countries/big tech are very low on that list after harassing me for years and treating me like shit. And I garantuee you, I will succeed at finding a polynomial time algorithm. There is no one else alive in this fucking world more determined then me to succeed at this. Fucking losers.

Update: Getting pretty close now. Just so you people know, had I been treated with respect, I would have been open to many things. You people made this an adversarial situation by treating me like shit. And yea, when I get threatened with a gun, kicked out of a part of the world I had lived for 4 years (3 years in Vancouver and 1 year in Redmond), forced to leave all my friends. Just to come back to this European shithole where I couldnt find employment my entire life.. yea.. ofcourse I'm going to entertain other options. Yes, had I just been treated with respect instead of this.. things wouldn't have had to happen this way. Its all far too late now. These last 2.5 years I wont ever forgive. Time is all we got on this planet, and there's so little of it. And you've wasted my time. So should I succeed soon.. then good luck.. I really dont understand the logic behind any of this. I'm guessing transphobia and mainstream politics has something to do with it. Plus by treating me this way you are also gambling that I'll never have another breakthrough in my life again. Because after these last 2 years... I want nothing to do with the west. I have very violent feelings towards all you people. 

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

Update: Alright changed a few more things..  will first sieve using a square modulus and then sieve at different quadratic coefficients with the modulus set as the odd exponent factors of any smooths found. The hope is that this setup will allow us to succeed earlier when factorizing very large numbers with huge factor bases.. Anyway.. suddenly got really bad headache. Gonna take a break.

Update: I just had an insight:

We can also add smooth candidates generated with zx^2-N to the linear algebra step... its still congruent mod N. It really doesn't matter.
I should try sieving both zx^2+N and zx^2-N tomorrow...  just in case there is something there we can leverage to finish faster. 

Aha.

So if N=4387

33^2+4387 = 5476  lets say we find this by sieving... now if we know the factorization of 33, we can use that as modulus and try to sieve for zx^2-4387 = 0 mod 33 where we construct zx from the smooth candidate we found. 
And we might find 74^2-4387 = 1089

Actually, that intuitively makes more sense. Let me try it tomorrow. Thats one hour of coding max. I will need to create a split again for the linear algebra step since we'll need to find a square for the left and right side. 
It just makes the most sense... hmm. I just got to try it.. it might just yield squares faster.. I wont know for sure until I try it. Maybe it will work like shit... but I cant risk not trying it and potentially missing something huge. It does make intuitive sense...

Like... it makes so much intuitive sense... I am internally cursing myself. This is literally what Ive been talking about in my paper for 3 years. Why just sieve one side when I can set it up like this....  ergh. What if... what if....  fuck it, time to get to work tomorrow. I can definitely check this in one day tomorrow. Easy. 





