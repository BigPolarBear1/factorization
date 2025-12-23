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

To do: It is still seriously bad. I'll fix stuff tomorrow. Right now its initializing the 2d interval with zeroes... then making a copy for each quadratic coefficient. Which isn't good. What I would really like is that as we go up in quadratic coefficient, we call a function like roll() on the entire 2d interval... let me investigate tomorrow if this would be possible. Because then we basically just roll() and sum() like there is no tomorrow... after that I need to use moduli for the interval step size and make sure we arnt bottlenecking with slow python indexing. I really want to reduce this to just roll() and sum(), I know its possible.. if I can pull that off I win, I just know it.

Update: Eureka! I believe I have a way to get roll() working. Nothing stops us from having multiple smaller "2d intervals"... so we can set the length to a modulus.. that way its all nicely contained. Something like that... let me give it a try tomorrow. roll() or death.

You know, I really hate microsoft. They destroyed everything that was good in my life. They also went after the only person in this industry who actually believed in me and supported me. And nobody at microsoft did anything. Fuck they probably patted the guy who threatened me with a gun on the back, wouldn't surprise me if they knew exactly who that guy was. All these people who work in security at microsoft, they all suck at their job, they are all incompetent. I dont even know why i'm mad at the fbi, bc its really microsoft that i'm mad at. But then again, the fbi are also losers.

