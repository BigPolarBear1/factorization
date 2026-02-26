Disclaimer: At no point did AI contribute anything to this research project. Copilot can't even calculate basic congruences without making mistakes. People who think AI can do novel math research are delusional. I did not use AI for the paper, and I did not use AI for any of the code. I attacked this problem, from the very basics 3 years ago, slowly burrowing deeper and deeper. 

#### To run from folder "NFS_WIP" and "NFS_Variant_Simple" (Will be adding more NFS related code in the coming days/weeks):</br></br>

Note: Experimental WORK IN PROGRESS.</br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 50 -base 500 -debug 1 -lin_size 100_000 -quad_size 1</br></br>

You can run the above command from either two. NFS_Variant_Simple is an intermediate step between QS and NFS and NFS_WIP is merging all my findings into a proper NFS algorithm...

!!!NFS related code is borrowed from: https://github.com/basilegithub/General-number-field-sieve-Python (note: Very impressively written, helped me big time, thanks)

To do for NFS_WIP:

UPDATE: I GOT IT! So I need to merge my findings from NFS_Variant_Simple with NFS_WIP. So right now NFS_WIP is bottlenecked by the fact that we must sieve with the same coefficients else we cannot combine smooths during the linear algebra step. HOWEVER... using my findings from NFS_Variant.. we see how we can transform smooths sieved with arbitrary coefficients to something of the shape x^2+N\*k. HOWEVER... we can use those exact same methods to transform smooths sieved with arbitary coefficients to something of the shape x^2+y\*x-N\*k ... where y is a coefficient of our choosing. The exact same principles apply.. its easy actually.

Update: Oh and unlike NFS_Variant_Simple its really just the polynomial value which must be square, the factorization of the constant of the discriminant (+ an offset determined by the linear coefficient we are sieving for) must just be smooth over the factor base. It doesn't need to be square. Just stressed today. I'll finish it tomorrow. But anyway.. so what we can do now probably is have 2 linear algebra steps.. where the first linear algebra step is to combine sieving results from different coefficients into a shape where we can use it to  plug into the final linear algebra step, which is the one for NFS using a quadratic character base and all that. I knew it could be done.. none of this makes sense. People would have known I was right.. this literally doesn't make sense.

Update: Wait, do I even need a non-zero linear coefficient to sucessfully take a square root over a finite field? If not.. .then its just straight up copy past from NFS_Variant_Simple.. that would be by far the easiest and most elegant solution... let me check... first however.. I'm going for a run. 


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

