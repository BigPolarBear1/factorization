Disclaimer: At no point did AI contribute anything to this research project. Copilot can't even calculate basic congruences without making mistakes. People who think AI can do novel math research are delusional.

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

Update: Just uploaded a quick PoC ... it just shows the condition that must be met to take the discrimimant when sieving the root and linear coefficient.
The pieces of the puzzle are starting to fall into place. Let me fix chapter 8 in the paper and start working toward an actual algorithm now... you really just want to combine stuff so that the condition is met where you are garantueed a smooth when taking the discriminant... which also means we dont need to worry about that damned factorization of the quadratic coefficient.

Update: Just got rid of what I previously wrote in chapter 8 in the paper and restarting that. I think the way I'm doing it now is way more promising. So basically we will create an algorithm now that does sieving with the quadratic polynomial (with non-zero linear coefficients) .. and then assembles them together such that if we call get_root() (see PoC) we get the same root mod N and modulo the smooth candidate. When this happens we know for sure the discriminant formula also yields a smooth... and then we dont need to factorize quadratic coefficients.... while still being able to do all the heavy lifting with the quadratic polynomials rather then the discriminant formula.

Update: I'll go for a run now. I'll work out some algorithm to combine quadratic polynomials such that the discriminant is also smooth in the coming days. Its much easier to generate extremely small values using the quadratic polynomials then the discriminant.. so if I can get it to work, that would be a big improvement over traditional Quadratic Sieve based algorithms. Plus I also don't need to worry about going to jail because apparently they decided to drop the case against me, given that I keep doing therapy. Maybe I can think clearly again now, because all this stress is throwing my head into chaos, it's not great for doing math.

Update: Anyway, tomorrow I have to run 30k to and from therapy again (or maybe I'm going to dissappear and run to Russia hahaha), so wont get much work done tomorrow. Threat intel losers stalking everything I upload online probably think my math is psychotic ramblings (because they suck at math). But I know what I've got. And as time passes, I get closer and closer. And I have a feeling that sieving with the linear coefficient is going enable me to finally finish my work for real now.. 
