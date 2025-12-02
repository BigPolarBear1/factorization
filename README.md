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


Update: Actually, the uploaded PoC isnt a good approach. I need to build smooths. I.e if N = 4387  just start like this 4388\*1^2-4387 = 1 ... and then we use all our number theory to reduce the size of that quadratic coefficient until it factors of the quadratic factor base by adding factors on the right side. This approach will be the only elegant solution. Nothing else is going to work. I build smooths, I succeed...if I cant build smooths... I cant punch past 110 digits. Until the end of the year I will attempt to build smooths. I know all the math now. I know all the mechanims. All the calculations. I know everything there is to know now. Now I make this work or I fail. One more month to succeed. I know I can do this, watch me perform real magic fuckers. Watch in awe as I destroy this fucking world! hahahahahaha.

Update: Feeling weak today. Ran 15k today, but felt like I had no energy in my legs at all. Also been struggling with severe bloating for nearly a year now (something I never struggled with in my entire life up until now). Maybe CIA fuckers tried to poison me. I guess I probably have some vitamin deficiency .. I should really see a doctor.. I just cant be bothered. Tomorrow.. let me work out the math to build smooths. Its much easier to just add smooth factors and resize the quadratic coefficient that way...that way we dont need to perform trial factorization on the smooths and just need to worry about the quadratic coefficient. Shouldn't be too difficult to work out... probably should have done it this way some time ago instead of wasting time with whatever the fuck I was doing these last few weeks.

Woudlnt even surprise me if these coward americans tried to poison me. For the same reason I was threatened with a gun in america while I was unarmed. All spineless fuckers. It doesn't even matter. My body is temporary anyway. Time to finish this math. Try harder to kill me fuckers. Losers all of you.

Trying to come up with a procedure to construct smooths:

So let us say N = 4387.

We start the algorithm with this: 4388\*1^2-4387 = 1
Next we perform trail factorization onf 4388. We notice it has 2^2 as factors.
Hence we can change it to: 1097\*2^2-4387 = 1
1097 is a prime number. So we can't further divide it. Now we need to increase or decrease the quadratic coefficient until we find a square again.
The next square we find needs to be as close as possible to 1097 as with each increasement or decrement the polynomial value/smooth candidate changes by 4.
We notice that 1098 has 9 as factor hence we can divide and multiply the root:
122\*6^2-4387 = 4

Just keep repeating this procedure until we meet the conditions to add it as a valid smooth (if both polynomial value and quadratic coefficient factor over the factor base).
Then I'll also need some way to prevent the algorithm from just hitting the same smooth over and over again. 

Anyway... it would be very similar to the uploaded PoC I guess. So I guess it wasn't too far off. Just need to improve some things.

God damnit man. This bloating, its like for months now. And I feel like its getting worse. Maybe I got like heavy metal poisoning or something. Maybe some fucking transphobe put mercury in my food or some shit. There's a lot of sick people in this world who would do that type of stuff. 

Whatever, even if people try to attack my health, because obviously ruining the rest of my life hasn't stopped me, I'll find a way. I will succeed through determination alone. I am a polar bear. The hunt never stops for as long as my heart beats. You'll have to kill me to stop me. But ofcourse, you're all weaklings and you stand no chance against me.
