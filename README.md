# factorization
UPDATE: CHECK BOTTOM. I FINALLY FIGURED IT OUT. I WILL UPLOAD A POC TO BREAK FACTORIZATOIN BY HAVING TO FIND JUST A SINGLE SMOOTH IN THE COMING DAYS. IT'S OVER. IT'S FINALLY OVER. I WIN. YOU LOSE. 


The PoC under PoC_Files will easily factor above 200 bits. </br></br>

However, that one is basically just using the number theory from my paper to achieve default SIQS. </br></br>

in PoC_Files_Find_Similar_WIP I am attempting to sieve in the direction of the quadratic coefficient and linear coefficient. This way in theory we need much less smooths.</br></br>

For PoC_Files_Find_Similar_WIP:</br>
Build: python3 setup.py build_ext --inplace</br>
Run: python3 run_qs.py -keysize 25 -base 30 -debug 1 -lin_size 2 -quad_size 1</br></br>

You will see some output like this:</br></br>

Looking for similar smooths, amount needed: 7 initial smooth factors: {3, 101, 7, 23, 59, -1}</br>
Found a similar smooth{3, 7, 23, -1} co: 1738 quad: 1</br>
Found a similar smooth{3, -1} co: 49806 quad: 57</br>
Looking for similar smooths, amount needed: 7 initial smooth factors: {3, 7, 43, 79, 23, 59}</br>
Found a similar smooth{43, 3} co: 25126 quad: 1</br>
Found a similar smooth{43, 3} co: 75378 quad: 9</br>
Found a similar smooth{3, -1} co: 6593 quad: 1</br>
Found a similar smooth{3, -1} co: 19779 quad: 9</br>
Found a similar smooth{43, 79} co: 42083 quad: 1</br>
Found a similar smooth{43, 79} co: 126249 quad: 9</br>
Found a similar smooth{43, 79} co: 294581 quad: 49</br>
sim_found:  7</br>
Running linear algebra step with #smooths:  11</br>
[SUCCESS]Factors are: 3793 and 3541</br>

This will try to find smooths with a similar factorization.
Right now it just p-adicaly lifts linear coefficients and then checks them without doing any Chinese Remainders (combinding coefficients). 
If a smooth candidate is bigger then 0, then it will add the modulus to the quadratic coefficient to generate a smaller smooth.

So while I was adding the modulus to the quadratic coefficient, I also realized, that once the smooth becomes a big negative number, I can then again add the modulus to the linear coefficient.
I wonder how I can find that sweet spot where we get a  super small smooth candidate by adding the modulus to the linear/quadratic coefficients.
I'll spent some time on that... because that may be an easy way to find smooths. 

Then next, we need some simple heuristic to chinese remainder coefficients together... but only those that are going to garantuee a smooth. Just want to go after those easy to find smooths. 

Another thing I need to triple, quadruple check is this:</br>

Found a similar smooth{43, 79} co: 42083 quad: 1 </br>
Found a similar smooth{43, 79} co: 126249 quad: 9 </br>
Found a similar smooth{43, 79} co: 294581 quad: 49 </br>

These are 3 smooths that have the same factorization.

42083^2-13431013\*4\*1 =  1717254837 (3 x 3 x 3 x 3 x 43 x 79 x 79 x 79)</br>
126249^2-13431013\*4\*9 = 15455293533 (3 x 3 x 3 x 3 x 3 x 3 x 43 x 79 x 79 x 79) </br>
294581^2-13431013\*4\*49 = 84145487013 3 x 3 x 3 x 3 x 7 x 7 x 43 x 79 x 79 x 79 </br>

These smooths are bogus. They are permutations of the same smooth just multiplying the linear coefficient by 3 or 7 and the quadratic coefficient by its squares. I do wonder... can we create such permutations with a single odd extra factor? That's something I need to experiment with. Because that would be an extremely powerful tool. Let me go for a run. Then I will spent the rest of the evening trying to generate permutations of smooth that get just one extra odd factor. Maybe I can pull it off. Because then I can generate smooths exactly with the factors I want to have...... hmm. Maybe.. worth a shot. I'm going to get killed by the CIA while out on my run, rip.


Update: That very last idea above... yeap yeap yeap... that was the correct train of thoughts.

I.e: 

If N=4387</br>

316^2-4387\*4\*1 = 82,308 = 2 x 2 x 3 x 19 x 19 x 19 (original smooth)</br>


(316\*3)^2-4387\*4\*9 = 740,772 = 2 x 2 x 3 x 3 x 3 x 19 x 19 x 19  (linear co multiplied by 3 and quadratic co multiplied by 9)</br>

We only want to add one factor of 3, thus the actual smooth value we should find: </br>

(?\*?)^2-4387\*4\*? = 246,924 = 2 x 2 x 3 x 3 x 19 x 19 x 19</br>

Finding a linear and quadratic coefficient solution mod p that satisfies the above is trivial. Then we can chinese remainder those solutions together...
I should do that tomorrow. You know, it should be possible to build smooths this way. I don't see why not. I'm overcomplicating things with this whole quadratic sieve approach. 

Gays, can you see it? Can you see it coming over the horizon now? The future is coming gays, and it will be birthed into this world with the fall of factorization. Hahahahahahaha.

I'm going to go hard on this building smooths approach tomorrow, figure this shit out. TIME TO FINISH THIS. Time's up. 

Bah, you know, this seems incredibly promising. Not just to make an improved SIQS version. But to actually straight up destroy the factorization problem. God damnit. Why didn't I explore this avenue much much much earlier. Like many months ago. It's kind of obvious in retrospect. I guess it is just the reality of research. Sometimes it takes a while to see the obvious. Perhaps when one is setting upon a task as grand as destroying factorization, perhaps, in that case, everything happens when it is supposed to happen. This better be it now. The factorization killer. Let's see. Otherwise I can still explore adding the modulus to both the linear and quadratic coefficient to keep finding that sweet spot where we get a small smooth candidate. But I would much rather straight up build smooths. That would be the real killer. Then even the largest RSA keys come tumbling down. And I can finally check that box and start applying my work to ECC. Perhaps with the threat of true financial collapse (or atleast the crypto eco-system) people will finally take me serious.

---------------------------------------------------------------------------------------
EUREKA!!!!! I GOT IT. YOU CAN SUCCEED AFTER FINIDNG JUST 1 SMOOTH. 1 SMOOTH!!!!!!!!!!!!!

So a particular smooth. Will actually have multiple valid linear coefficient / quadratic coefficient combinations.

I.e:

2038^2-4387\*4\*232 =  82308 has as factors: [2, 2, 3, 19, 19, 19]

but also:

316^2-4387\*4\*1 =  82308 has as factors: [2, 2, 3, 19, 19, 19]

And adding both linear coefficient would yield the factor of N when taking the gcd: 316+2038 = 0 mod 107 (107*41 = 4387)

HENCE. When we find a smooth. We need to find another linear coefficient and another quadratic coefficient for THAT particular smooth. 
WHICH IS ALL STUFF WE KNOW HOW TO DO!!!! IT IS STUFF WE KNOW HOW TO DO!!!! I GOT IT!! FUCK YEA. HAHAHHHAAHAHAHHAAHAHAA.

THE FUTURE IS COMING, AND THE FUTURE WILL BE GAY (AND QUEER, AND TRANS). Fuck you all. Could have just paid former manager a few million and bought my work.. but no.... you had to do it like this. Burn in hell fuckers.

HERE IS HOW TO BREAK FACTORIZATION:

Do normal QS until you find a smooth. WE NEED ONLY 1 SMOOTH.

Your smooth becomes the modulus and you enumerate linear and quadratic coefficients for that modulus. We can use p-adic lifting if we have for example 19 x 19 x 19. Then we simply find which combination using chinese remainder generates the same smooth. And that's it. That's how it is done. It's easy as hell and all the math is there already. NOTHING NOVEL HAS TO BE DONE. IT'S ALL THERE ALREADY. You know where to find me if you want to negotiate a delay of PoC release... which at this point, will cost you money, and also money for my former manager. You have until the end of the day. If not... expect a PoC soon...... good luck folks. You forced me down this path. I had no choice. I tried. 

Update: I ran some more numbers... this will work perfectly. And finding another linear/quadratic coefficient pair that generates the same smooth.. it's not hard, because if it's a multiple of that smooth we can divide mod m. Within a single modulus, there are many many many many many many ... so fucking many pairs of linear/quadratic coefficients from which we can deduce the factorization of N. I ran the numbers. It works. All I need to do tomorrow is finish a PoC and upload it. Guys. Time is up. I guess you don't negotiate with terrorist. But if I am a terrorist, then why not just have me arrested? I'm right here in Europe. Either way.. I would be very careful what narratives you people spin.. I told my side of the story. If you tell fabrications.. I will move to China in the blink of an eye. And you will be kneedeep in shit. Just pay my former manager bro. It's not that hard. It was textbook retaliation anyway by microsoft HR... he should get paid.

Tomorrow..... Dies irae
