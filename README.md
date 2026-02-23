Disclaimer: At no point did AI contribute anything to this research project. Copilot can't even calculate basic congruences without making mistakes. People who think AI can do novel math research are delusional. I did not use AI for the paper, and I did not use AI for any of the code. I attacked this problem, from the very basics 3 years ago, slowly burrowing deeper and deeper. 

#### To run from folder "NFS_Variant_Simple" (Will be adding more NFS related code in the coming days/weeks):</br></br>

Note: Experimental WORK IN PROGRESS.</br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 50 -base 300 -debug 1 -lin_size 100_000 -quad_size 10</br></br>

NFS related code is borrowed from (Mainly sieve(), just to demonstrate we can sieve the same way): https://github.com/basilegithub/General-number-field-sieve-Python (note: Very impressively written, helped me big time, thanks)

I was trying to port some of my insights to the number field sieve algorithm. NFS_Variant_Simple is an intermediate step between QS and NFS. I suspect it may be possible to use the number field sieve approach while sieving many different polynomials rather then just a single polynomial and its multiples. Or at the very least different quadratic coefficients since it relates to multiples of the constant (see CUDA_QS_Variant which tries to leverage different Quadratic Coefficients).

Update: Uploaded NFS_WIP, run with: python3 run_qs.py -keysize 14 -base 50 -debug 1 -lin_size 100 -quad_size 10  
This is work in progress. I'm still trying to figure out the math to get it to sieve many different polynomials. Mainly by changing the quadratic coefficient, because we already know how we can translate that to multiples of N. The current PoC currently just sieves different quadratics where the constant is multiples of N until it finds one where x+y is square, the polynomial value is square and everything in the quadratic character base is square... then it takes a square root over a finite field. I'm now trying to figure out how to get some linear algebra working again while sieving all these different polynomials. I know it can be done and I know I have a good foothold already into a breakthrough. Almost there now.

Update: I'm being an idiot.... its actually very simple now to finish my work. Just trying to make things more difficult then they actually are inside my head lol. *sigh* Let me begin constructing my final version now. Thats the problem with my brain.. it creates all these wildly abstract and complicated thoughts.. and sometimes I just need to do the easy and straightforward thing thats already right infront of me.

Update: You know, some days, it almost felt like someone was trying to reach out to me across space and time. And I know its just in my head due to trauma and the nearly complete social isolation in the aftermath of that trauma. Plus the only reason why anyone would try to reach out across space and time would be if one day I broke factorization. Which is unlikely, even though I'm going to give it my best shot. It would be funny if I was actually breaking factorization and the people reading this are probably spooked as hell now. Anyway, on the offchance that time travel is real, make it so that all my friends and family will have good fortune and no more bullshit like microsoft did to my former manager, thanks, oh and make everyone involved in firing my former manager absolutely miserable... like, destroy their lifes, I dont care. 

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

