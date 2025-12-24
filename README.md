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

To run: python3 QScu.py -keysize 50 -base 500 -lin_size 10_000 -quad_size 10_000</br></br>

Note: You need to run from the host. I'm using wsl2 and ubuntu with CUDA support. I don't believe this will work from hyper-v.

UPDATE: I did some brain storming today. I need to use array slicing.... it will be by far the fastest method to shift/roll rows as we change quadratic coefficients. Then the only computation heavy thing that has to be offloaded to the GPU will be taking the sum() of that 2d interval. And in addition we can precompute by how much we have to shift rows when we change quadratic coefficients.. because if we do it in the hot loop we end up with a lot of duplicate computations, especially for small primes. 











You know, I really hate microsoft. They destroyed everything that was good in my life. They also went after the only person in this industry who actually believed in me and supported me. And nobody at microsoft did anything. Fuck they probably patted the guy who threatened me with a gun on the back, wouldn't surprise me if they knew exactly who that guy was. All these people who work in security at microsoft, they all suck at their job, they are all incompetent. I dont even know why i'm mad at the fbi, bc its really microsoft that i'm mad at. But then again, the fbi are also losers.

I realize after what happened in 2023, I dont want to have social attachments anymore. I will become a recluse, and when the time comes that everyone has long forgotten about me, I'm headed into the Arctic mountains, never to be seen again. Thats going to be the rest of my life and nothing will change that. My life was great before microsoft tore it apart. I had friends and people who supported me. Never again do I want to experience that pain. Better to just be alone. It doesn't even bother me, like when I'm alone in the mountains, I dont really feel alone, there's something about the natural world that brings me peace. 
