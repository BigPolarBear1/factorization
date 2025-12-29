Disclaimer: At no point did AI contribute anything to this research project. Copilot can't even calculate basic congruences without making mistakes. People who think AI can do novel math research are delusional.

Disclaimer extended version: Fuck AI. Tech bros (especially those at microsoft) are idiots who dont understand algorithms and computing. Literally, if I see some of the shit people write about AI, glorifying that trash, it just makes me mad. Literally, the average tech worker is an idiot (especially those at microsoft). There is nothing worse then someone in an overpaid job, with average intelligence who thinks he's better then everyone else (the archetype to which elon musk belongs, and most people at microsoft).. people like that are going to destroy this world. None of you know shit about anything in this world, so stop pretending like you do. Talk is easy. Doing is hard. Maybe do more and talk less. Fucking AI tech bro losers. I piss on you all.

NOTE: Starting 2026, none of my research will be published. Only people who treat me with respect will be allowed access to my work. And NATO countries/big tech are very low on that list after harassing me for years and treating me like shit. And I garantuee you, I will succeed at finding a polynomial time algorithm. There is no one else alive in this fucking world more determined then me to succeed at this. Fucking losers.


#### To run debug.py" (Prints the linear and quadratic coefficients to solve for 0 in the integers, for use with my paper):</br></br>

To run: python3 debug.py -keysize 12

This basically creates a system of quadratics. Solving them mod p is easy. But there is only one root solution (the factor of N) which solves the system for 0 for any mod p (aka solves it in the integers). Figuring out how to exactly do this quickly is still an ongoing area of research for me. And if a polynomial time algorithm for factorization exists, it is likely done by solving this system of quadratics. Finding a polyomial time algorithm is my ultimate goal, as this would make progress toward solving p = np as well. 

#### To run from folder "CUDA_QS_variant" (WIP):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run:  python3 run_qs.py -keysize 120 -base 2000 -debug 1 -lin_size 1_000_000 -quad_size 1_00</br></br>

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

UPDATE: Bah, moved some shit around again and put the roots in arrays on file to be memory efficient. It's a bit slow atm, but I rearranged it so I can finally do what I've been trying to do from the start.... which I will implement tomorrow. 

Ok, so tomorrow here's what I'll do:

-quad_size is going to become huge. We want to use enourmous values here. As much as our hard drive allows us. 
In build_database2interval(), this loop:

while j < len(quadlist):

We will change that, and simply iterate 1 to prime. Because all quadratic coefficients mod p will have the same root. So we only need to do 1 to prime at most. And we just bulk build all the quadratic coefficients mod p with the same root. Just construct a massive amount of sieve intervals in one go. Going to save us from repeating calculations over and over again.
And because while i < len(primeslist): is now the outer loop, we can preconstruct per prime, two large sieve intervals for both roots, just slice as needed and add them to our actual sieve intervals.. so we make use of vector addition in the GPU.

That should be easy to implement tomorrow.. its not that much work. Going to take a break now bc I have been working all day and I regret not having gone running today. I just feel agitated if I skip running. 

These last few days have been difficult mentally. This will either work or it wont. There is moments where I feel trapped, like I cant breath, and very dark thoughts to end it all flash through my head. If this works, I'm never going to forgive people. I'll be 36 in 3 months.. my life has just flashed by. I had a few happy years in Vancouver for a little while, with friends.. but that has been taken away again by microsoft. Aside from that, I'm always working day and night, I'm nearly 36 and aside from those 3 years, I never had time to live life. Just always working and grinding. And for what? I dont have a career. The only people who respect me are from countries the west is hostile towards. Life is a shitshow. I just want to do the impossible, and for once in the history of mankind, shout loudly "a polar bear was here". It's the only way any of this would have been worth it.
