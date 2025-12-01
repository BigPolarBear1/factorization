DISCLAIMER: AT NO POINT IN MY 2.5 YEARS OF RESEARCH DID I USE AI OR HAS AI IN ANY WAY HELPED WITH MY RESEARCH. I'M PUTTING THIS HERE BECAUSE I KNOW PEOPLE WILL ACCUSE ME OF IT. 

Note: I am antifa, leader of the cipherpunk, fuck the FBI department. 

Note2: Once Improved_QS_Variant is finished, I will rewrite most portions of the paper completely from start. I keep finding myself "maturing" out of these papers.. because I keep learning and advancing with my math skills. There's so much stuff in that paper that I just want to completely rewrite because even to me it looks way too amateuristic now.

Note3: I am utterly broke so if anyone wants to make donations to keep this research going a little longer: big_polar_bear1@proton.me email me... thanks.

Note4: If what I'm trying to do currently doesn't work, I will completely abandon the quadratic sieve type of approach. Because using the number theory I worked out, another approach would be is to directly find a factor of N. Since using the factors of N solves the quadratic for 0 for any mod m when the correct linear and quadratic coefficients are used. So this will be the next avenue I will explore. And thinking a little more about this... actually let me check something today also really quick... I woke up way too early today.. but I suddenly had this insight which I must absolutely verify today just incase I overlooked a very obvious solution using all the number theory I worked out...

#### To run debug.py" (Prints the linear and quadratic coefficients to solve for 0 in the integers, for use with my paper):</br></br>

To run: python3 debug.py -keysize 12

This basically creates a system of quadratics. Solving them mod p is easy. But there is only one root solution (the factor of N) which solves the system for 0 for any mod p (aka solves it in the integers). Figuring out how to exactly do this quickly is still an ongoing area of research for me. And if a polynomial time algorithm for factorization exists, it is likely done by solving this system of quadratics. Finding a polyomial time algorithm is my ultimate goal, as this would make progress toward solving p = np as well. 

#### To run from folder "QS_variant" (Standard SIQS with our number theory as backend):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 220 -base 20_000 -debug 1 -lin_size 10_000_000 -quad_size 1</br></br>

See below for an improved way of performing what this PoC does.. I'll delete this Proof of Concept once the PoC for Improved_QS_Variant matures a little more<br><br>

Note: With a large enough -base and lin_size this PoC will find smooths for 110 digits. Albeit very slowly, but this is a highly unoptimized cython PoC. However, to push beyond that into novel terroritory for Quadratic Sieve-based algorithms we need to use quadratic coefficients and p-adic lifting, and that is what the PoC below (Improved_QS_Variant) will be for. 

#### To run from folder "Improved_QS_Variant" (Implements more of my number theory and attempts to succeed with fewer smooths by using p-adic lifting):</br></br>

To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 80 -base 10_000 -sbase 1_000 -debug 1 -quad_size 10_000</br></br>

Update: Calling it a day. I added a sieve interval for the quadratic coefficient. In the code you will see a list called many_primes .. we use this to mark large squares in the sieve interval for the quadratic coefficient.. this way we can divide the quadratic coefficient by any squares it find and reduce the bit length. The advantage of this sieve interval setup is that we dont need root calculations. I had cases where both the quadratic coefficient and smooth candidates (poly_val) ended up being of very small bitlength. It still happens too rarely... but its those cases, the cases where we manage the reduce the size of the quadratic coefficient with large squares... those cases we need to hunt down to attack very large numbers. One idea I'm playing with is to  just keep generating primes for a while and keep marking its squares in the sieve interval.. and using even larger sieve intervals. That should help to find those cases. Anyway, I'll try to get that to work for the rest of the week. If I cant get it to work, I'm going to sell some 0day exploits to make some money because I'm so broke, I cant even buy things anymore. I'll have to move my math to my spare time... and I'll probably go back to eploring some of my initial ideas of attacking it by solving a system of quadratics instead, because if this doesn't work then I feel like I've pretty much exhausted the quadratic sieve approach.

So for tomorrow what I'll do is: 

1. Much bigger sieve intervals. Atleast a gig worth of elements. (and use array.array)
2. Many many many many more primes... we dont need to hold them in memory... can just chunk it... just got to keep improving the quality of that sieve interval by marking it with more and more squares. We need those cases where we can reduce the size of the quadratic coefficient by as many bits as we possibly can... because factoring the quadratic coefficient is now the bottleneck... the size of poly_val is now simply decided by the mod_mul variable (smaller = smaller poly_val).. so the only thing that matters now is that quadratic coefficient.

Then the ultimate test will be, given a 400-bit number... we can set mod_mul to a low enough value so that it generates smooth candidates (poly_val) below 40 bits. And then hopefully we can reduce the quadratic ccoefficient enough in size so we can factor it over a factor base without requiring a ridiciously sized factor base (the many_primes list isn't part of the factor base since we only use it with even exponent.. hence they dont increase the amount of required smooths). If I can achieve this one thing... my work is done and ready to publish. Lets see tomorrow.. I'm so fucking broke. I really need that breakthrough right now. Also I taunted the FBI so much this last year, I better succeed now with my math, otherwise I will be the loser. 

Update: Oh yea, one thing I should also explore is if I cant just build smooths in the integers from the ground up with everything I know now. It might actually be possible now. 
