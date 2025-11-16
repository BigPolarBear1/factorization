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
2. In find_similar, we should use a sieve interval
3. Once the sieve interval is implemented, we should mark the sieve interval with p^2 for primes that arn't in our modulus, so those end up being ignored in the linear algebra step... this will allow the linear algebra portion to potentially succeed much earlier.
4. Also use exponent 2 in the original sieve interval at quadratic coefficient = 1... for the primes outside the modulus we use to fill out the interval... this will help generating smooths with fewer negative exponent factors. Additionally I need to optimize that lifting code... because I'm fairly sure I don't need to brute force it like I'm doing now.. 

I'll begin attacking that to do list tomorrow. But first, I need to clean up all of the code... it's becoming convoluted with a lot of unneeded code. Frankenstein code due to making too many edits.

The code I uploaded today is really bad. I need to change and add a lot of things in the coming days. But I know I'm finally on the right track and that my work's completion is imminent now. 

Update: Also quickly added the trick the generate permutations of linear coefficients in the find_similar() function. Next this needs a sieve interval.. that will be my goal for tomorrow. Then the algorithm should really start to ramp up....

Update: Note-to-self even for the sieve intervals at quadratic coefficient = 1... we should mark p^2 in the sieve interval aswell.... the trick is really to generate reasonably small moduli made up from smooth factors with negative exponents. And seeing as we're going to spent alot of time in that lifting code, I should optimize that as well instead of doing that dumb bruteforcing. 

Update: I guess tomorrow once that sieve interval is done... then we're nearly there. Deep sadness, for all that was taken from me. All the injustice that was done to friends. Finish this... then just keep going.. this won't be the end of my work, this will be just the beginning. I still believe there is a much greater discovery to be found, and even now I am barely scratching the surface. I will violently throw myself against this, I will throw away my life in the pursuit of solving this. Nothing else matters anymore. People have ruined all that was good in my life. And it cannot be fixed anymore. 

Update: High levels of stress today. It's funny. I have been unemployed for 2 years. Yet I desperately feel like I need a vacation. And with vacation I mean, literally spent months alone in the arctic living in a tent. I miss the arctic. It's giving me serious stress... because I know the consequences of finishing my work. But this far in... it's too late now.... my range of choices has been limited to disclosure or drop everything I'm doing and stop working on it. I know the consequences of this. It's not so much myself I worry about, it's everyone around me. I ever meet the people who forced me down this path, knowing full well I was working on this, I'll fucking kill them. 

It's funny, some people even have the gal to portray me as incompetent because I don't appear entirely mentally stable. Fucking losers. I push myself harder then any of them. I dropped out of highschool. The only way to do the things I do, is to push myself harder and sacrifice more then everyone else. All you fucking losers don't understand shit what it means to push yourself this hard. You're all too busy with your social games and politics, while I'm here, day after day, grinding math, putting in the hours. I don't need shitty comments from lesser humans who don't know a thing about what I do. Fuck you all losers. Quit your jobs, you all suck.

Update: Finding it incredibly difficult to focus today. I'll just implement the sieve interval and go running. Then next week I can do all the final optimizations. If people knew how much, in the last 10 years I have sacrificed, to do the things I do... just to end up getting treated like complete shit by microsoft. 10 years of pushing myself mentally to my limits, and nothing to show for it. Everyone else my age is further ahead in life. I have no savings. I have no driver's license. I have no social life, due to being displaced and losing work visas. I can't get a job anymore, because infosec is a loser industry. They dont care about talent. They are not playing this game to win. It's all ego and social connections. The west has stopped being a meritocracy a long time ago.
