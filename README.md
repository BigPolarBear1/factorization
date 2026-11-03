The only thing that will fix this is if my father survives (and makes a reasonable recovery) and all these incels who send me hate mail over the past 8 years, if they disappear never to be found again. And if everyone around me that I care about, including friends from the past (like my former managers), if nothing bad happens to any of them and all these people get to live a life completely insulated from the horrors of this shit world. I'm done with people just dying all around me. If my father dies, then congrats, you have succesfully succeeded at destroying me, but I suspect that was the objective all along. People think I'm some type of idiot, that I'm blind. I see everything you fuckers. 

Note: This project is pro-lgbtq, so fck off haters.

#### References (re-used many of the core number theoretical functions from these PoCs to fit my own algorithm): 
https://stackoverflow.com/questions/79330304/optimizing-sieving-code-in-the-self-initializing-quadratic-sieve-for-pypy
https://github.com/basilegithub/General-number-field-sieve-Python 

#### Paper
Math paper is a work in progress.
Note: This paper was written at different stages in my  research project, and as such I will need to go back to refactor some of the older chapter as I now have a better understanding of the entire structure.

#### To run from folder "fu_algo":</br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 50 -base 5000 -debug 0 -lin_size 10 -quad_size 1

Pronounced the fuck you algorithm

Update: Will complete somewhere this week. It already works and showcases the math from the paper, including final chapter. Whats missing is hensel and proper sieve intervals. 

Its fucking hard to focus. This week will be grim. Discussions about whether or not to let me father die. I just want my father back. I dont understand why he keeps developing complications. It's statistically unlikely the bad luck he had, its like they are trying to kill him on purpose in that hospital. He should have never gone to a belgian hospital. Belgium is a country of sadists and incompetent bureaucrats, human life has no value in belgian society, unless youre friends with all the right people. I despise this country so much. I've wasted my entire life here and now it has also taken my father from me.

Update: Fairly sure you can use the factorization of the original root to figure out how often to add the modulus to the final square root step. Let me think..

Update: This also has to hold true, where plugging in that shared root, also reveals the k value of the other sibling solutions. Let me see if I can use that:

(2\*41\*26)^2-(1708\*2)\*(2\*41\*26)+4387\*665 = 41\*4387</br>
(2\*41\*26)^2-(1708\*2)\*(2\*41\*26)+4387\*156\*4 = 0</br></br>

(2\*41\*26)^2-(424\*2)\*(2\*41\*26)+4387\*41 = 665\*4387</br>
(2\*41\*26)^2-(424\*2)\*(2\*41\*26)-4387\*156\*4 = 0</br></br>

Having a hard time focusing today. But in the above example, the correct coefficients will solve for 0 mod 41 (here k = 41 just happens to be a factor of N, but its coincidence, I should probably change this example in my paper eventually), so you can use these divisor of the original k value whose sibling solution we are trying to find. So you got this mechanism to lift yourself from a local solution to a global solution. I can see, I know it will work. Just too tired right now. I'll do it tomorrow. My father was a better person then all these fucking ghouls in this world. I just feel rage. When this is done, I'm not ever sharing anything ever again. Hell, I'm not ever talking to anyone ever again. Except to tell them to fuck off. Just this fucking rage. Fucking ghouls and troglodytes destroying everything that is beautiful in this world. Drop dead you stupid fuckers. And I swear, if people dont fuck off and keep harassing me like those fucking incels have been doing for the last fucking 8 years, I'm going to fucking kill someone. Also why the fuck are these fucking incels who have been harassing me for years still alive huh? Do something then you fucking cunts. My father was the pacifist, I am not. I believe some people are so ghoulish they deserve to die.

#### To run from folder "psieve":</br>
To build: python3 setup.py build_ext --inplace</br>
To run:  python3 run_qs.py -keysize 150 -base 10_000 -lin_size 10_000 -quad_size 1</br>

note: Up to 120-bit is near instant now.

Some ongoing experimentation.. this mixes SIQS and Coefficient Sieve. I either find a way to improve this, or take elements from this for my hybrid_sieve.
This is my main research project, trying to open up a completely new way of achieving factorization following the literature from my paper.

Will try and accelerate this research with Claude while I spent time studying abstract algebra. If there is anything novel in my work and it is close to finished, then Claude may actually succeed in bridging the final gap. Worth a shot. 

Update: Managed to push it to 140-bit now..

Update: Claude wont be able to finish this. I need to go back to representing this the same way, as a system of quadratics, that led to all my insights in my first year of research, but look at it again with my more developed knowledge now.... and specifically go hunting for a way to achieve direction computation of a sibling root that psieve finds by sieving. Direct computation would mean a complete break in the factorization problem and im going to figure it out.. just need to go back to doing research the oldschool way like I did in my first year. If anyone wants access to my work, I urgently need $ and am willing to sell any breakthrough that may happen: big_polar_bear1@proton.me

#### To run from folder "Coefficient_Sieve" (For use with the paper):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run:  python3 run_qs.py -keysize 40 -base 50 -debug 1 -lin_size 10_000 -quad_size 100</br>

Just demonstrates the math from the paper using quadratics. For educational purposes. And rather then taking a square root over a large prime we can also just calculate the discriminant. But this demonstrates the interesting relation between these quadratics and the factors of N.

#### To run debug.py" (Prints the linear and quadratic coefficients to solve for 0 in the integers, for use with my paper):</br></br>

To run: python3 debug.py -keysize 12

This basically creates a system of quadratics. Solving them mod p is easy. But there is only one root solution (the factor of N) which solves the system for 0 for any mod p (aka solves it in the integers). Figuring out how to exactly do this quickly is still an ongoing area of research for me. And if a polynomial time algorithm for factorization exists, it is likely done by solving this system of quadratics. 

(PERFORMING HEATHEN CYBER RITUALS IN THE BELGIAN MOORLANDS TO DESTROY MY ENEMIES)

