The only thing that will fix this is if my father survives (and makes a reasonable recovery) and all these incels who send me hate mail over the past 8 years, if they disappear never to be found again. And if everyone around me that I care about, including friends from the past (like my former managers), if nothing bad happens to any of them and all these people get to live a life completely insulated from the horrors of this shit world. I'm done with people just dying all around me. If my father dies, then congrats, you have succesfully succeeded at destroying me, but I suspect that was the objective all along. People think I'm some type of idiot, that I'm blind. I see everything you fuckers. 

Update: What a shitshow. My father really ended up in the worst ICU of belgium. And if I wasnt broke I could atleast make sure everything is done for him. They are just giving up on my father without doing proper diagnostics like EEG of the brain. You know if anyone wanted to buy my loyalty, they could just make sure my father had access to world class care... like right this moment before he dies in that shithole belgian ICU. I cannot lose my father. It wasnt even his time yet. He was murdered by belgian incompetence.

Note: This project is pro-lgbtq, so fck off haters.

#### References (re-used many of the core number theoretical functions from these PoCs to fit my own algorithm): 
https://stackoverflow.com/questions/79330304/optimizing-sieving-code-in-the-self-initializing-quadratic-sieve-for-pypy
https://github.com/basilegithub/General-number-field-sieve-Python 

#### Paper
Math paper is a work in progress.
Note: This paper was written at different stages in my  research project, and as such I will need to go back to refactor some of the older chapter as I now have a better understanding of the entire structure.

#### To run from folder "fu_algo":</br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 40 -base 5000 -debug 0 -lin_size 10 -quad_size 1

Pronounced the fuck you algorithm

Update: I will fix it shortly. The final chapter isn't quite correct. It starts out correct but then overcomplicates everything.
So any root solution for a given coefficient (the way I calculate them) encodes some multiple of the factors of N. Ofcourse the problem is that this is unbounded. 
However if modulo every prime there is some root such that we subtract the same multiple of N it becomes 0 modulo some large modulus, which can be a finite field for example.. like mod p^e, to make it easier to search for those.. then we have a solution inside that finite field that we can take the GCD with. Actually its quite elegant and simple.. and I'm an idiot because I had come so close already many months ago. Especially when I was exploring lifting and singular roots. But here we are now... better late then never. Going for a run. I'll start fixing everything now.

If I get murdered on my run. You know I was right and I still want my father to be saved. Shithead cunts. 

Update: ill try and upload some code tomorrow. Its hard to focus. I think i may legit be losing my mind bc of my father being in the ICU for months now. My father really was the one person I cant miss in my life right now. Like some cruel twist of fate that ill never understand. Anyway, I understand the math now.. so this week for sure. Going to finish it.

Update: Today is the day. Will try and upload by the evening. These last few years have beem a shitshow. Ive learned a truth about this world. Its that humans as a collective, they really hate outcasts, especially when thosr outcasts are able to do things none of these shitheads could even if they spent their entire life trying. I have lived in paralell with human society for most of my life now. I may aswell be a different species, a bear. Good thing im tall and big amd people know to stay away from me. I feel sorry for people like me who arnt tall and big. Their lifes must not be great. I will fix this world. I will bring into this world things no normal human can ever conjure. And that will be my big fuck you.

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

