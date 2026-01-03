One final time before I do what I must, looking for work, big_polar_bear1@proton.me. Cant work for big tech bc of what happened in 2023. Cant travel to the US anymore bc I would get arrested (hence cant work for US companies either). Cant say I didnt try when all the infosec shitheads start pointing fingers again.

Disclaimer: At no point did AI contribute anything to this research project. Copilot can't even calculate basic congruences without making mistakes. People who think AI can do novel math research are delusional.

NOTE: Starting 2026, none of my research will be published. Only people who treat me with respect will be allowed access to my work. And NATO countries/big tech are very low on that list after harassing me for years and treating me like shit. And I garantuee you, I will succeed at finding a polynomial time algorithm. There is no one else alive in this fucking world more determined then me to succeed at this. Fucking losers.

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

UPDATE: FUCK AMERICA. FUCK ANY COUNTRY ALIGNED WITH AMERICA. Pete hegseth is a fucking pussy. Coward. Its funny how america only ever invades much poorer countries (and still loses). Biggest cowards on this entire world. I piss on america.
