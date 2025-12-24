Disclaimer: At no point did AI contribute anything to this research project. Copilot can't even calculate basic congruences without making mistakes. People who think AI can do novel math research are delusional.

NOTE: Starting 2026, none of my research will be published. Only people who treat me with respect will be allowed access to my work. And NATO countries/big tech are very low on that list after harassing me for years and treating me like shit. And I garantuee you, I will succeed at finding a polynomial time algorithm. There is no one else alive in this fucking world more determined then me to succeed at this. Fucking losers.


#### To run debug.py" (Prints the linear and quadratic coefficients to solve for 0 in the integers, for use with my paper):</br></br>

To run: python3 debug.py -keysize 12

This basically creates a system of quadratics. Solving them mod p is easy. But there is only one root solution (the factor of N) which solves the system for 0 for any mod p (aka solves it in the integers). Figuring out how to exactly do this quickly is still an ongoing area of research for me. And if a polynomial time algorithm for factorization exists, it is likely done by solving this system of quadratics. Finding a polyomial time algorithm is my ultimate goal, as this would make progress toward solving p = np as well. 

#### To run from folder "QS_variant" (Standard SIQS with our number theory as backend):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 220 -base 20_000 -debug 1 -lin_size 10_000_000 -quad_size 1</br></br>

See below for an improved way of performing what this PoC does.. I'll delete this Proof of Concept once the PoC for Improved_QS_Variant matures a little more<br><br>

Note: With a large enough -base and lin_size this PoC will find smooths for 110 digits. Albeit very slowly, but this is a highly unoptimized cython PoC. However, to push beyond that into novel terroritory for Quadratic Sieve-based algorithms we need to use quadratic coefficients and p-adic lifting, and that is what the PoC below (Improved_QS_Variant) will be for. 

#### To run from folder "Improved_QS_Variant" (Building smooths from the ground up... see paper on clues to finish the code):</br></br>

To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 14 -base 20 -debug 1 -lin_size 100_000 -quad_size 1</br></br>

Note: Broken for now... this is some work in progress research

#### To run from folder "CUDA_QS_Variant"  (QS variant using mainly vector operations and CUDA support):</br></br>

To run: python3 QScu.py -keysize 100 -base 500 -lin_size 100_000 -quad_size 10_000</br></br>

Note: You need to run from the host. I'm using wsl2 and ubuntu with CUDA support. I don't believe this will work from hyper-v.

Update: I've quickly added some more stuff. This needs a lot more work though. First things first though... tomorrow I will set it up so we create a sieve interval for atleast 4 quadratic coefficients (encode it in a 64-bit number) in one go. And I just start optimizing everything. And also make sure we keep making copies of intervals and all that to a minimum.. I have some ideas for that. Producing huge intervals should be lightning fast. 

Update: Aha! Eureka! So the way this should be setup eventually is you have a worker process/thread just having the GPU create sieve intervals and do sum()... and write the result of the sum to disk. Then we have multiple worker processes just reading this and performing trial factorization. If all is well, one worker thread using the GPU should create more then enough workload for a bunch of workers doing trial factorization. Plus the more we can optimize creating sieve intervals.. the higher we can put the threshold variable.. making it less about trial factorization and more about fast indexing.

ps: There is also some bug in the code because i'm using zx^2+n rather then zx^2-n this time around.. so I need to check some things so our smooth candidates arn't too big. I may just swap it around again.. but it shouldn't matter too much in theory.

