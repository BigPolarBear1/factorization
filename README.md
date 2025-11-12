DISCLAIMER: AT NO POINT IN MY 2.5 YEARS OF RESEARCH DID I USE AI OR HAS AI IN ANY WAY HELPED WITH MY RESEARCH. I'M PUTTING THIS HERE BECAUSE I KNOW PEOPLE WILL ACCUSE ME OF IT. 

Note: I am antifa, leader of the cipherpunk, fuck the FBI department. 

#### To run from folder "QS_variant" (Standard SIQS with our number theory as backend):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 200 -base 6000 -debug 1 -lin_size 100_000 -quad_size 1</br></br>

See below for an improved way of performing what this PoC does.. I'll delete this Proof of Concept once the PoC for Improved_QS_Variant matures a little more<br><br>
#### To run from folder "Improved_QS_Variant" (Implements more of my number theory to also take advantage of quadratic coefficients):</br></br>

To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 160 -base 4000 -debug 1 -lin_size 100_000 -quad_size 100</br></br>

Update: I very quickly implemented the basic ideas that I developed in recent months using quadratics of the shape zx^2-N. I need to check for bugs still though and that entire code needs to be vastly improved. All of it is less then ideal right now. As I improve it, it should eventually overtake the original QS_Variant in terms of performance. Since we have much more fine grain control over the size of our smooths this way as we can now displace the location of our parabola created by the sieve interval with a linear multiplier (aka the quadratic coefficient).


##### #To do: 
1. Quadratic coefficient should not be squares, since then we generate the same parabola as x^2-N. We should check if restricting to prime quadratic coefficients vs non-square composites makes a difference in terms of smooth diversity (smooths that are not just going to create trivial factorizations).
2. Remove what I'm calling the "quad interval" just rely on jacobi symbols insteads...
3. When 2 is implemented, building the iN datastructure should also be restricted to primes found at quadratic coefficients 1. That will drastically improve the building time there..
4. The big ticket item will be to be "smart" about shifting the parabole, such that we can garantuee smaller smooths. I should rework that together with number 7. because right now as the quad co goes up, we drift toward bigger smooths.. which is not at all what we want.
5. Re-implement large-prime variant, since that is broken for now
6. We can reduce the required amount of smooth by reducing the size of the factor base used for the quadratic coefficients.. 
7. We also need to work on the whole "generate modulus" logic. This works fine for standard SIQS, but for us, the smooth size is determined by the root^2 \* quadratic coefficient. So we need to completely rework all of that.
8. Coefficient lifting to square moduli! But this will need to tie in with whatever we change in step 7.

------------------------------------------------------------------------------------------
##### #random rants below

Update: I'm going to take a break for the rest of the day I think. I still got headache, and tomorrow I have to go to court ordered therapy for sending angry emails to the FBI (fucking losers, just mad because I'm better, smarter and stronger then them, pathetic). Since I'll go run 30k to and back from therapy.. that will take up most of my day tomorrow. But later in the week I'll start working through that to do list. The math is all there now. It all works. But especially number 4 and 7 on the to do list need to be done to gain a performance advantage.. but once implemented that performance advantage should be massive.. as it really allows for great control over the placement of our parabola.

Anyyyywaaaaaayyyyy, sad day for the NSA, FBI, CIA, DIA and all their NATO lapdogs. HAHHAHAHAHAHAHHAHAHAHAHHAHAHAHAHAHA. Destroyed, by a single polar bear. I cannot understate how weak you people are. It's all truely pathetic. Should all quit your jobs. 

Ps: During the weekends I have been working also on an algorithm to attack DLP. I tried to work on it in secret, but now that that algorithm is also reaching maturity, I will drop it soon.

Pps: pete hegseth is a fucking loser and a piece of shit. And I am waging war against america. Don't mistake what I'm doing for anything else. I was a guest you country, and you people insulted me and pointed a gun at me while I was unarmed while your leaders cries about "dudes in dresses". And damn nato to hell for being pathetic lapdogs. Fucking weaklings. Your all fucking pathetic. Fucking losers.

Update: Got back from therapy. Didn't sleep much last night and ran 30k, fucking tired right now. I'll began finishing the math tomorrow. So tomorrow actually, I really want to implement a heuristic so we displace the parabola in the right direction and yield smaller smooths. Proving I can gain an advantage that way is the most important thing.. everything else is just polish. 

You know, I guess I am kind of worried my therapist thinks I'm completely psychotic. I doubt they have many people talking about 0days, breaking cryptographic schemes.. meeting the fbi (and sending angry mails ot them) .. emailing with the Chinese.. etc. Combined with all the traumatic shit from when I was young. My life literally sounds like the type of life an insane person might come up with. But I'm just talking about my life the way it happened. If they don't believe me and want to force medication on me, then that's where I'll draw the line, quit therapy and let my case go to court then. I promise you I will flee the country before I ever set foot inside a courtroom. I'm very comfortable living in the outdoors in a tent in all seasons of the year. If I decide to dissappear, you wont ever catch me again.

Remember... the FBI and Microsoft knew I was working on factorization long before I mentioned it publically. What happened is not an accident. They made the calculated risk to assume I would fail. I truely tried, but in the end I was left with no other choices. 10 years ago, none of this would habe happened.. The world has gone mad and now you pay the price for it. 
