DISCLAIMER: AT NO POINT IN MY 2.5 YEARS OF RESEARCH DID I USE AI OR HAS AI IN ANY WAY HELPED WITH MY RESEARCH. I'M PUTTING THIS HERE BECAUSE I KNOW PEOPLE WILL ACCUSE ME OF IT. 

Note: I am antifa, leader of the cipherpunk, fuck the FBI department. 

Note2: It's basically done now, the NFS variant is nearing completion soon, nothing novel research wise has to be figured out now. I just need to refactor a few things and I'm done.
I have been unemployed for 2 years. I am looking for work. I would like some way where I can have enough time to attack ECC and DLP next. In addition, for a payment, I can also show people how I appraoched factorization from start to finish. Because the method of how I conduct my research and consistently have success is just as valueable as the research itself. And I follow the same patterns over and over again, no matter what subject matter I attack. I don't have a very high IQ and I'm not especially bright or good at doing numbers in my head. But I have insights into research that help me succeed over and over again. I'm willing to teach people how to do it, as long as it enables me to have the finances to continue doing my research aside from that. Contact: big_polar_bear1@proton.me.

#### To run from folder "QS_variant" (Standard SIQS with our number theory as backend):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 200 -base 6000 -debug 1 -lin_size 100_000 -quad_size 1</br></br>
This version is somewhat done, but it doesn't achieve much of an advantage, as the main strength of my work will be the NFS_Variant, which performs NFS with reducible quadratic polynomials<br><br>
#### To run from folder "NFS_variant" (Number Field Sieve with our number theory as backend using reducible quadratic polynomials mod m):</br></br>

I just uploaded polarbearalg_v32.py</br>
To run: python3 polarbearalg_v32.py -keysize 20</br>

I started slowly digging in today for about an hour or two (having bad brain fog again). Just thinking about it deeply. And in essence to finish my algorithm I need to solve the following problem:

(this is a small toy example)

poly1: 3^2+17\*3 = 60 	(3\*20)   </br>  
poly2: 7^2+17\*7 = 168  ( 7\*24)    </br>  

poly1\*poly2: 21^2+17\*21 = 798 (21\*38)   </br>  

Lets say we multiple poly1 and poly2 together. The output of the new polynomial will be 21*38. The 21 we can easily predict... but that 38... how can I predict the residue of that mod p such that I can transform that part into legendre symbols? That's all I need. I'm sure with everything I have worked out so far, I can solve this. I solve this last issue, and I'm done. The moment I solve this, a few hours later I'll have a PoC ready to upload and it's over. I'll bash my head against this tomorrow the entire day. 

Update: *sigh* I'm a idiot. It's this fucking brain fog. It's just the way my life is right now, I feel claustrophobic, feel like I can't breath or think. Almost there... doesn't matter. 

1089^2+4387\*1089 = 2442^2   (5476\*1089) .... so in this example which I covered in the paper...  I know where that 5476 comes from... thats just the square of the root which produced 1089 (74^2-4387 = 1089). Yea, I know how to setup my legendre symbols now. Time to finish it tomorrow... 

Update: Eureka. After I wrote the sentence above... my brain couldn't shut down.. so I opened my laptop and ran the numbers and figured out exactly how it works. hehehehe. I should be good now to finish my work tomorrow... hopefully the brain fog wont be too bad. I guess it comes and goes. My brain isn't very great. I'm fairly certain I ended up with brain damage in my late teens/early twenties. But then again, my brain has let me do amazing things, despite the limitations... so there is no point in thinking too much about this shit. Had I grown up perfectly sheltered, never exposed to violence and self-destructive behavior... yea sure, maybe I wouldn't have dropped out of highschool, maybe I would have gotten a college or uni degree.. and probably ended up with some 9-5 job. This is what people don't understand, you need diversity of thought.. sadly, there isn't many places in the world where that is valued.. it definitely wasn't valued at microsoft. Big tech is not at all like it portrays itself. It's just tech bros and wealthy people with good degrees... it's not a place that welcomes people who don't fit in. It's pretty much the opposite. For a few years at microsoft... with my managers.. yea sure... but as soon as HR and the lawyers caught on one of my managers had hired a weirdo like me.. they fired both me and him. The corporate world is a disgusting place and even thinking back to all the trauma microsoft caused, I feel intense anger. microsoft is fucking disgusting company and their entire leadership better make sure to never ever cross my path, becaues I will become violent.

I am serious. I ever cross paths with microsof people, be it msrc, or msft threat intel or whatever, I will physically attack them. I hate these cowardly piece of shit. I hate them all. All these people are responsible. These are the same people who downplayed everything I did my entire career. They are probably also responsible for what happened in 2023. I have no doubt. They just couldn't stand it, especially after I had success in openssl i 2022.

And above all, what they did to my former manager. After he had worked there for 27 years. For 4 years, I had the best friends I ever had, like my former manager.. no other moment in my entire life even comes close... and they destroyed it all, and they destroyed the life of one of my best friends (my former manager). I'm waging total war against microsoft. I will see that company destroyed before I'm dead. I garantuee this. I don't care if I have to destroy the entire internet in the process. That's just how total war works. And if you work microosft, then you're the enemy. Shouldn't work for microsoft unless you want to be my enemy. Yea.... better make sure our paths never cross. fuckers.

You know what I truely hate the most? Living in a world full of ape brained morrons. I cant stand people. I liked the friwnds I used to have, like my former manager and my teamlead, and everyone I got to know when I lived in Vamcouver.. but its these morrons like these HR people, these lawyers, these idiots pointing guns at me.. these regular dumb people who through their low intelligence make this world a nightmare. Sure, these lawyers might be considered intelligent by society, but society is very clueless about what intelligence actually means, bc most people with fancy degrees are straight up morrons of the worst kind who bring nothing but pain and misery into this world and I hope you all drop dead.

Update: Alright, I had 7 hours of sleep. I'm going to finish this today. Almost 3 years of research. I used to think spending 6 weeks to find an openssl bug was completely insane... but this is on another level. I'm starting to worry my life is too short, to do all the math things I want to do. 

Update: Let me update the paper first with these recent findings. It will help organizing my thoughts. 

If this was somebody's classified algorithm, they must be spooked out of their mind, HAHAHAHHAHAHAHAHHAHAAH. Get fucked losers. 

Like, literally, the only way I can explain why nobody is trying to stop me is because you are all so fucking spooked, you are frozen like deer staring at headlights. Weaklings. Quit your jobs losers.

Update: Alright gays, I have updated the paper. Now let me figure out how that math would work for other primes, other then those used to generate our coefficients. I probably just have to recalculate the root modulo some other prime...

UPDATE: I JUST RAN THE NUMBERS!!!!!!!! YES, FOR OTHER PRIMES, YOU JUST RECALCULATE THE ROOTS. EUREKA!!!! I WILL FINISH THE PAPER TONIGHT, AND THEN THE CODE TONIGHT OR TOMORROW. MY WORK HAS FINISHED!!! I need a break few hours first though and also go for a run, I am stressed out of my mind, you wouldn't believe.

I'll try to just finish the paper tonight, show how to construct legendre symbols for the other primes using the example with the 3 quadratic polynomials I added earlier to the paper... and then just do the code tomorrow. I am stressed out of my mind right now. What the fuck does it even matter. Ah man. I have no other choices left. Just got to keep going down this road, until the bitter end. 
