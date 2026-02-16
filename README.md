Disclaimer: At no point did AI contribute anything to this research project. Copilot can't even calculate basic congruences without making mistakes. People who think AI can do novel math research are delusional. I did not use AI for the paper, and I did not use AI for any of the code. I attacked this problem, from the very basics 3 years ago, slowly burrowing deeper and deeper. 

#### To run from folder "NFS_Variant" (Will be adding more NFS related code in the coming days/weeks):</br></br>

Note: Experimental WORK IN PROGRESS.</br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 50 -base 1000 -debug 1 -lin_size 100_000 -quad_size 1</br></br>

NFS related code is borrowed from: https://github.com/basilegithub/General-number-field-sieve-Python (note: Very impressively written, helped me big time, thanks)

This code is a work in progress. I'm trying to merge some of my findings from https://github.com/BigPolarBear1/factorization/tree/main/NFS_Variant_simple into this one.

To run from NFS_Variant_simple, which is an intermediate step between QS and NFS: python3 run_qs.py -keysize 50 -base 300 -debug 1 -lin_size 100_000 -quad_size 10

I will continue working on this in the weekends. But I dont know how long it is going to take. I want to be able to sieve using different quadratic coefficients in the NFS implementation like I'm doing in my _simple version and CUDA_QS_Variant. Because together with the precalculated hashmap, that would allow for very strong sieving. But I dont know how long that could take to implement, if at all possible. Could be days, or could be weeks. I am broke, I havnt seen my friends in years because of what Microsoft did. I have nothing here in belgium. And every day, the dark thoughts become stronger and stronger and stronger. And if I dont do a hail mary now and hope to find a good 0day exploit that I can sell.. so I can maybe go hike in the arctic, find some peace again... I'm just going to end it, because there is absolutely nothing worth living in my life anymore. Everything was destroyed. Thats just how it is, and thats just the kind of world we live in. It doesn't matter all the work I have done in my career. People were jealous from the start. From my very first bugs.. I was an intruder into an elitist industry.. and they managed to remove me again, I never belonged in this industry.

I started this project a little less then 3 years ago. I didn't know basic algebra. Never finished highschool. I've learned a lot. I dont care what people think. Most people wouldn't have the courage to do what I did, like some mediocre dipship with a university degree and 6 figure salary with less CVEs then me, pointing fingers and laughing. Thats the type of shit world we live in. You cant get anywhere with determination alone. Everything is politics. They even went after the only person to ever support me in this industry and ruined his career as well after he worked at microsoft for 27 years while portraying me as absolutely incompetent. And there is not a single person in this entire world who spoke up, because everyone wanted to see me fail from the start, like those loser stalkers who had stalked and harassed me for 7 years. And thats just the truth of it. Thats the type of world this is. Just an endless meat grinder, destroying hopes and dreams and devouring people alive.

UPDATE: OH SHIT. Wait a minute. I couldn't fall asleep. And I suddenly got struck by lightning. So when using NFS, if the factorization of the root doesn't matter. Then whats stopping me now from just lifting all my solutions in the hashmap to even powers and just grabbing cases where the polynomial value and zx+y is square and doing a square root over a finite field. RIGHT? WHATS TO STOP ME? NOTHING IS!! HOLY FUCK. FUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUCK! I GOT IT! I need to sleep first.

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

