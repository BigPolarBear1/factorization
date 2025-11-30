DISCLAIMER: AT NO POINT IN MY 2.5 YEARS OF RESEARCH DID I USE AI OR HAS AI IN ANY WAY HELPED WITH MY RESEARCH. I'M PUTTING THIS HERE BECAUSE I KNOW PEOPLE WILL ACCUSE ME OF IT. 

Note: I am antifa, leader of the cipherpunk, fuck the FBI department. 

Note2: Once Improved_QS_Variant is finished, I will rewrite most portions of the paper completely from start. I keep finding myself "maturing" out of these papers.. because I keep learning and advancing with my math skills. There's so much stuff in that paper that I just want to completely rewrite because even to me it looks way too amateuristic now.

Note3: I am utterly broke so if anyone wants to make donations to keep this research going a little longer: big_polar_bear1@proton.me email me... thanks.

Note4: If what I'm trying to do currently doesn't work, I will completely abandon the quadratic sieve type of approach. Because using the number theory I worked out, another approach would be is to directly find a factor of N. Since using the factors of N solves the quadratic for 0 for any mod m when the correct linear and quadratic coefficients are used. So this will be the next avenue I will explore. And thinking a little more about this... actually let me check something today also really quick... I woke up way too early today.. but I suddenly had this insight which I must absolutely verify today just incase I overlooked a very obvious solution using all the number theory I worked out...

#### To run from folder "QS_variant" (Standard SIQS with our number theory as backend):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 220 -base 20_000 -debug 1 -lin_size 10_000_000 -quad_size 1</br></br>

See below for an improved way of performing what this PoC does.. I'll delete this Proof of Concept once the PoC for Improved_QS_Variant matures a little more<br><br>

Note: With a large enough -base and lin_size this PoC will find smooths for 110 digits. Albeit very slowly, but this is a highly unoptimized cython PoC. However, to push beyond that into novel terroritory for Quadratic Sieve-based algorithms we need to use quadratic coefficients and p-adic lifting, and that is what the PoC below (Improved_QS_Variant) will be for. 

#### To run from folder "Improved_QS_Variant" (Implements more of my number theory and attempts to succeed with fewer smooths by using p-adic lifting):</br></br>

To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 80 -base 10_000 -sbase 1_000 -debug 1 -quad_size 100_000</br></br>

Update: Before I finish this PoC, let me do something different today ... let me write a script that is going to show the number theory of the paper. I suddenly have this nagging feeling I may have overlooked a much easier approach to solve this rather then trying to ram my work into a quadratic sieve type of algorithm. 
