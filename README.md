One final time before I do what I must, looking for work, big_polar_bear1@proton.me. Cant work for big tech bc of what happened in 2023. Cant travel to the US anymore bc I would get arrested (hence cant work for US companies either). Cant say I didnt try when all the infosec shitheads start pointing fingers again.

Disclaimer: At no point did AI contribute anything to this research project. Copilot can't even calculate basic congruences without making mistakes. People who think AI can do novel math research are delusional.

NOTE: Starting 2026, none of my research will be published. Only people who treat me with respect will be allowed access to my work. And NATO countries/big tech are very low on that list after harassing me for years and treating me like shit. And I garantuee you, I will succeed at finding a polynomial time algorithm. There is no one else alive in this fucking world more determined then me to succeed at this. Fucking losers.


#### To run debug.py" (Prints the linear and quadratic coefficients to solve for 0 in the integers, for use with my paper):</br></br>

To run: python3 debug.py -keysize 12

This basically creates a system of quadratics. Solving them mod p is easy. But there is only one root solution (the factor of N) which solves the system for 0 for any mod p (aka solves it in the integers). Figuring out how to exactly do this quickly is still an ongoing area of research for me. And if a polynomial time algorithm for factorization exists, it is likely done by solving this system of quadratics. Finding a polyomial time algorithm is my ultimate goal, as this would make progress toward solving p = np as well. 

#### To run from folder "CUDA_QS_variant" (WIP):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run:  python3 run_qs.py -keysize 160 -base 2000 -debug 1 -lin_size 10_000_000 -quad_size 100</br></br>

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

UPDATE: ERGH. I woke up in the middle of the night feeling very restless. So I started reading about SIQS. They use zx^2+yx+c where y^2-n = zc   .... I'm basically sieving y^2-n = zc... but we should be sieving zx^2+yx+c. ................lol. I am a MORRON! I AM STUPID! I WILL QUIT MATH FOREVER AFTER THIS! All my heroes, Galois, Fermat, Legendre and Euler would be laughing at me if they were still alive. I will try to sleep for a few more hours.. then work non-stop to correct this. I dont get it. Because I know I read about SIQS earlier.. but somehow, reading about it tonight, it suddenly clicked... I guess the previous time I read about it, I just hadn't progressed enough with my math skills to correlate what I was reading with what I was doing myself.

GOD DAMNIT! I should have made the connection months ago! GOD DAMNIT. It was so obvious. I gave my enemies too much time to react. They must have really thought me a fool who wouldn't figure this out... and shamefully, they were almost correct. I figured it out.. but too late. I need one day to correct this. Prepare yourselves!
