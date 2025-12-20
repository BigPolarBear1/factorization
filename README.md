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

EDIT2: Woops, I already see it now... I confused myself by looking at an example that doesn't have a valid mapping in the hashmap.. when both do, then its very easy and straight forward to combine then. Anyway.. I'll continue tomorrow.

1\*40^2-150\*40+4387 = -13</br>
1\*36^2-158\*36+4387 = -5

If we perform sieving (just iterating x and y), lets say we find the above smooths. And in both cases the linear coefficient has a valid mapping to the quadratic coefficient in the hashmap we calculated for mod 13 and mod 5. Can we combine the two above polynomials to find one that generates -65? 

If we chinese remainder residues 150,-150 mod 13 and  158,-158 mod 5 and look for nearby linear coefficients and roots that generate -65 we find for example:

1\*42^2-148\*42+4387 = -65

However, is the existence of some nearby polynomial always garantueed that produces the product? And that question is probably solved with the work I did last week, where I was adding factors by multiplying/dividing coefficients and roots. So my goal now is to write some function, given two quadratics whose output factorizes over the factor base and has a valid mapping in the hashmap, combine then and return a new linear coefficient and root that produces the product. If I can solve this one thing, then I absolutely know how to finish everything else and my work should be completed within days... final sprint I guess.. lets hope this is it for real now.

You know what would be great? Be able to mentally take some vacation and just do things that arn't stressing me out. Like being in a place where I know things will turn out alright. But I'm not there and aside from those short few years when I was living in Vancouver Canada and had a job and friends, I just felt miserable and angry, isolated myself and grinded out 0days. People have no idea about the sacrifices I made to achieve all the things I did. And what do I have to show for it? I have nothing right now. Literally everyone else my age is further ahead in life. I'm literally at fucking rock bottom. The only thing that could make my life worse is becoming homeless. I've worked so hard my adult life, just to get treated like shit and see the people who support me get treated like shit (microsoft firing my manager because he once gave me a promotion). I will defeat factorization, no matter the cost. I succeed or my life is over. There is nothing else anymore. People havnt got the faintest clue how hard I work and how much I sacrifice, and this is why I despise humans, because its those same humans, who are too much of a coward to put everything on the line, who are always pointing fingers. I despise these lesser humans. People who point fingers, harass and/or laugh at others, they are lesser humans. I look down on these humans. 

I'm going to get some sleep. Lets see if I can get a PoC ready tomorrow that combines quadratic polynomials. Once that works... thats basically everything, thats literally my last obstacle.. perhaps soon I will know some peace again mentally for a little while. 

I was thinking back about the past. The very first time I entered the US was probably for a job interview in Austin, Texas. I remember at the time there was this woman at the hotel who was like a government liason, who funnily enough was doing business with the same company I was interviewing at, and when we got to the lobby, the security guard at that company got all excited because that government liason woman had the same name as some character in the tv show "stargate" and started nerd-ing out. Americans were really fucking weird and goofy. Back then I had the impression it was some weird place where stuff actually "happened" unlike snobby bureaucratic europe.  And even back in 2019 when I went to interview in Redmond for Microsoft, I had so many conversations with random strangers. Someone even walked up to me like 5km outside of Redmond, saying she recognized me from the hotel and asked if I wanted to get drinks later... although she then hooked me up with a cyber security recruiter from Japan incase Microsoft didn't end up hiring me.. which made me a bit suspicious... but fuck, it was America, some goofy place where shit happened, who cares if shit is suspicious. But then, after living in Canada, when I finally moved to the US in 2023.. I didn't recognize that country.. it wasn't that goofy place anymore. Culture wars had destroyed that country. It wasn't a goofy place anymore. It was a very hostile place. And after what happened in 2023... even if I wasn't banned from the country and on a 100 terrorism related lists, I still wouldn't for the rest of my life ever set foot into that country again.  America has fallen. And it is a damn shame, because the rest of the world sucks as well. I never could get employment in Europe. I only ever got  paid money by american companies. I wish there was a goofy place where people could be free to be who they are... but there is no such place in the world anymore. And even in Europe, maybe its a safe country.. but for someone like me, there are no oppurtunities in Europe.. Europe is not a place where someone without degrees and education can achieve anything. There is no such place anymore in the world. 

I think if I were given the chance to get a visa and move to Russia or China and be allowed to do my work without constraints... I would do it.. without doubt. I dont care if they arn't very lgbtq friendly countries. To me, being respected for my work is more important then anything. In europe, people just look at me with suspicion. People here think I'm psychotic and I cant find a job. I would go literally anywhere if they could arrange a visa and employment. Fuck, I'll go work for Venezuela and help them in their fight against america. One of the few people to every really respect me for who I was as a person and for what I was capable of, work-wise, was my former manager whom microsoft retaliated against. Sometimes it fills me with deep sadness to not have that in my life anymore. Now I'm just some unemployed loser again, living with my parents, unable to get a job... which is funny, because people used to say I was one of the best in the world at what I did.. now the only people who still seem to respect me are the Chinese and Russians and I wish I had the courage to just take an airplane and go where I'm still respected. It surely is not here.
