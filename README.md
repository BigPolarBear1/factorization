DISCLAIMER: AT NO POINT IN MY 2.5 YEARS OF RESEARCH DID I USE AI OR HAS AI IN ANY WAY HELPED WITH MY RESEARCH. I'M PUTTING THIS HERE BECAUSE I KNOW PEOPLE WILL ACCUSE ME OF IT. 

Note: Remember, as this now comes to completion in the coming days. Both the FBI and Microsoft knew I was working on factorization when I had only been working on it for a few months and nobody else in the entire world knew I was working on it. I tried desperately to avoid the path of public disclosure as I would have much rather sold my work. But at every turn people made it impossible. This would not have happened 10 years ago. An oppurtunity presented itself, it was literally handed to you on a silver platter.. and instead you people all treated me like shit and gambled on me failing with my project. I ask you, where has your common sense gone? I hate what this will do to everyone around me, this is never what I wanted.

Note: I am antifa, leader of the cipherpunk, fuck the FBI department. 

#### To run from folder "QS_variant" (Standard SIQS with our number theory as backend):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 200 -base 2000 -debug 1 -lin_size 100_000 -quad_size 1</br></br>

See below for an improved way of performing what this PoC does.. I'll delete this Proof of Concept once the PoC for Improved_QS_Variant matures a little more<br><br>
#### To run from folder "Improved_QS_Variant" (Implements more of my number theory to also take advantage of quadratic coefficients):</br></br>

To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 100 -base 1000 -debug 1 -lin_size 10_000 -quad_size 1</br></br>

Update: I have began implementing a strategy that would allow the algorithm to finish much sooner. See the find_similar() function. I need to improve it... just uploading my work in progress for the day before I tae a break now.

Strategy is a follows:

1. Use square moduli to find a smooth at quadratic coefficient  = 1
2. When we find a smooth, divide by the square modulus and using the factors with negative exponents, construct a new modulus. Because of using square moduli in step 1, the bit size of this one will be relatively low.
3. Now find smooths with these factors at different quadratic coefficients.
   

To do:

1. In step two, if the new modulus is too small, we can us p-adic lifting to increase the bit length (edit: but this rarely seems to case, usually it is slightly too large... in which case the solution is also simple, we can divide the modulus by a prime... and then rerun again dividing by a different prime.. so all primes end up covered... but I'll do this when my algorithm is complete)
2. ~~In find_similar, we should use a sieve interval~~
3. Once the sieve interval is implemented, we should mark the sieve interval with p^2 for primes that arn't in our modulus, so those end up being ignored in the linear algebra step... this will allow the linear algebra portion to potentially succeed much earlier.
4. Also use exponent 2 in the original sieve interval at quadratic coefficient = 1... for the primes outside the modulus we use to fill out the interval... this will help generating smooths with fewer negative exponent factors. Additionally I need to optimize that lifting code... because I'm fairly sure I don't need to brute force it like I'm doing now.. 

Update: I have quickly implemented a sieve interval.
The next important thing to do is in the construct_interval_similar() function, the primes outside the modulus we should only use even exponents there to mark the interval. This will increase the odds that we end up with smooths having just the factors we want. And if we're going to do a lot of p-adic lifting.. improving that code also becomes vital. In addition... the size of the sieve interval should be fitted properly according to the length of the new modulus. Eventually all that find_similar stuff... it needs to execute fast.. because there we are just interested in smooths that potentially complete the orginal smooth. I don't want to waste too much time in there otherwise.. so I need to think how to improve that.

Anyway.. enough work for today. My head is just not in the right space today to do math or coding. I'll go for a run soon. Next week things should start to happen quickly... there is still a lot of optimization work to be done, but we'll start seeing the results of this aproach quickly now...

Update: So the whole gist of what I'm trying to do is not to find smooths faster, its to succeed at the linear algebra step with much fewer smooths.  I tried this back in September, but now that I've made more progress with the number theory I see now how it can be done. Anyway.. just depressed today. Still haunted by everything that has happened. And I know there is a very high chance that this will actually work, I'm feeling very confident about it this time... and I also know the absolute shitshow it will unleash. I would have sold my work... wouldn't even have asked for much. Probably just a salaried job would have been enough if I couldn't straight up sell it. But now, it is all far too late. It really does feel like this was one big set-up... that someone wanted this outcome. Because my life got so suddenly and violenty thorn to shreds, and not just my life. also that of one of my best friends and one of the few people in this industry who believed in me and supported me.... that especially is what keeps tormenting me. And I think it will torment me until the day I die. On the macro level, it's just a drop in the ocean of human drama... people go through things much worse... I guess thats why I'm depressed, the world just feels like a relentless meatgrinder. I'm glad I had a chance to stand in the Arctic, experience that place.. had I only known this miserable human world, I wouldn't have survived this long. I want to see more of the Arctic, even if I have no more ambitions in this human world.

Update: Alright, goal for today is to improve the lifting function. Just found an algorithm for hensel's lifting for quadratic roots on stack overflow which appears to be better then my bruteforce method... so lets see if I can get it implemented. Then I can just lift my entire factorbase once while precomputing everything. Then tomorrow I'll finish the algorithm by marking the sieve interval with squares in the find_similar functions... and after that all that is left is optimizing everything.

Update: Alright, implemented hensel's lifting lemma in lift(), wasn't too hard.. wonder why I hadnt realized this could be simplified sooner.. guess it just never was a bottleneck up until now. And also just generally math maturity .. I understand things a little bit more now I guess. Ok, I'll go run and take a break. Tomorrow I'll mark the sieve interval in find_similar() with just even exponents ... so we have a high chance of hitting just the factors we need. After that.. optimize optimize optimize. 
You know.. this really needs to come to an end soon. Because I really need to have some time where I can afford to just study math the formal way. Learning stuff by messing around like I do now.. it's great, because it helps to develop a really deep and intuitive understanding.. but I need to study for math too so I can consume advanced literature and learn a these little tricks that could help me. I just need to finish this first, otherwise my life is just stuck. 
