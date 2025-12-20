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

UPDATE Saturday 20 DEC: I'll upload a better PoC shortly as I work out the details.

After spending the afternoon thinking about this, I think the following must be done now:

We must combine quadratic polynomials, such that z\*x\*2 = y.
I need to work out the exact math for this, so that when we combine them, no extra factors that don't already occur in the quadratic outputs that we are combining get added. 
Ideally rather than using quadratic polynomials to perform sieving such that the discriminant is smooth, we should combine them such that we generate a square and can take the GCD.. so somehow that combining of quadratic polynomials must be fed into guassian elimination over GF(2). So that's what I'll begin exploring now. I'm fairly sure it's do-able though. This does feel like the right approach. The main issue will be approaching it in a way that makes sense.. rather then desperately trying to mimic steps in the number field sieve algorithm.. because that's just going to confuse me.

I'm 99% sure I can pull it off now. In the last week, when I was messing around with adding factors to the quadratic output and doing modular division on that quadratic coefficient... now that I understand all those moving parts, I should have the understanding to pull it off now, if it is at all possible...  Anyway.. time to dig in now. I hate that it is taking this long... but then again, all these in retrospect idiotic detours I've taken in my research.. they helped to shape the understanding I have today.. and perhaps a lot of it, I could have learned by reading some math textbooks, but it wouldn't have given me the intuitive understanding I have now.. and I know from experience and a long career doing vulnerability research... that intuitive understanding triumphs everythings.

Update: Here is an example of what needs to be done:

EDIT: BTW, THIS IS A BAD EXAMPLE, SINCE THIS DOESN'T HAVE THAT VALID MAPPING IN THAT QUADRATIC TO LINEAR COEFFICIENT MAPPING. JUST AN EXAMPLE TO ILLUSTRATE WHAT I'M TRYING TO DO I GUESS.

1\*40^2-150\*40+4387 = -13</br>
1\*36^2-158\*36+4387 = -5

If we perform sieving (just iterating x and y), lets say we find the above smooths. And in both cases the linear coefficient has a valid mapping to the quadratic coefficient in the hashmap we calculated for mod 13 and mod 5. Can we combine the two above polynomials to find one that generates -65? 

If we chinese remainder residues 150,-150 mod 13 and  158,-158 mod 5 and look for nearby linear coefficients and roots that generate -65 we find for example:

1\*42^2-148\*42+4387 = -65

However, is the existence of some nearby polynomial always garantueed that produces the product? And that question is probably solved with the work I did last week, where I was adding factors by multiplying/dividing coefficients and roots. So my goal now is to write some function, given two quadratics whose output factorizes over the factor base and has a valid mapping in the hashmap, combine then and return a new linear coefficient and root that produces the product. If I can solve this one thing, then I absolutely know how to finish everything else and my work should be completed within days... final sprint I guess.. lets hope this is it for real now.

You know what would be great? Be able to mentally take some vacation and just do things that arn't stressing me out. Like being in a place where I know things will turn out alright. But I'm not there and aside from those short few years when I was living in Vancouver Canada and had a job and friends, I just felt miserable and angry, isolated myself and grinded out 0days. People have no idea about the sacrifices I made to achieve all the things I did. And what do I have to show for it? I have nothing right now. Literally everyone else my age is further ahead in life. I'm literally at fucking rock bottom. The only thing that could make my life worse is becoming homeless. I've worked so hard my adult life, just to get treated like shit and see the people who support me get treated like shit (microsoft firing my manager because he once gave me a promotion). I will defeat factorization, no matter the cost. I succeed or my life is over. There is nothing else anymore. People havnt got the faintest clue how hard I work and how much I sacrifice, and this is why I despise humans, because its those same humans, who are too much of a coward to put everything on the line, who are always pointing fingers. I despise these lesser humans. People who point fingers, harass and/or laugh at others, they are lesser humans. I look down on these humans. 
