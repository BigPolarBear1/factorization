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
To run: python3 run_qs.py -keysize 80 -base 5_000 -sbase 1000 -debug 1 -quad_size 10_000</br></br>

Update: Uploaded some work in progress. So right now, using mod_mul variable we can change the ratio of the bitlength between the smooth candidate and the size of the quadratic coefficient. The lower it is, the larger our quadratic coefficient is and the smaller our smooth candidates produces by the polynomial. I.e if we generate 400 bit numbers and set mod_mul to 0.05, then we would get smooth candidates of about 30-bit and quadratic coefficient of 300+ bit ... since the quadratic coefficient is linear, we can add or subtract the modulus to it aloooooooooooooooooot of times without really affecting the size of the smooth candidate. Now the only step that remains... if we have these 300+ bit quadratic coefficients .... we need to add or subtract a modulus x amount of times  and reduce the size of the quadratic coefficient with large squares... because then we can move those over to the root instead. I dont think sieve intervals would be the best way to do this, since if we're only adding a 20 bit modulus to a 300 bit quadratic coefficient.. then the range to find large squares without affecting the bit length of the smooth candidates produced by the entire polynomial is thus enourmous.. and I should take advantage of that large search space. I just got to think. There has to be some good heuristic for this.

I guess this is similar to what I was doing recently, trying to reduce the length of the smooth candidate with large squares by adding the modulus to the root x amount of times (i.e the sieve interval in QS_Variant)... but now instead we are changing that strategy to reducing the length of the quadratic coefficient with large squares, since that is in theory more approachable since we dont have to worry about quadratic growth.

Update: Doing some testing. So I can keep the smooth candidate small. Then it just becomes a matter of finding large squares such that New_z = z + a*mod ... where new_z is divisible by a large square and where 'a' must be bounded by some value. Then we can divide New_z and move it over to the root. 'a' needs to be bounded so we dont end up increasing the size of the smooth candidate. But since we are dealing with linear growth, that bound can be very large. It does work... I just need to write proper code and a proper strategy for maximizing the size of those large squares.
