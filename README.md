#### References (re-used many of the core number theoretical functions from these PoCs to fit my own algorithm): 
https://stackoverflow.com/questions/79330304/optimizing-sieving-code-in-the-self-initializing-quadratic-sieve-for-pypy
https://github.com/basilegithub/General-number-field-sieve-Python 

No AI was used for this research. I tried using it from time to time, but it is always a waste of time. AI is only having success at math bc of context aware bruteforce for answers. It cannot generate novel thought. I hate AI, waste of time.

Note: Psieve PoC is incomplete, "a" needs to be non-square in psieve_process_interval() and I'm close to getting it to work.. effectively letting me derive the factorization of N without having a fully square B-smooth.

Update: Ok fuck it. Since I already wrote the solution in the readme before. You calculate potential square roots to "a" when "a" is not square in Z/p where p divides the linear coefficient. You can either lift p or use CRT to combine multiple p... as long as the square is larger then the square root of N.. then we are guaranteed to find as solution. So basically we can garantuee succesful factorization by only looking at polynomial solutions for a modulus of roughly the size N^1/4. I wonder if I can further divide that with quartics or if CRT is going to reveal a more favorable pattern as I've only experimented with hensel so far. I'll check it out tomorrow. 

I'm quitting uploading my work here so it doesn't just end up getting scraped by AI and someone stealing the credits.

I'm also quite tired of the desperation for years and unemployment, especially with my father battling cancer and being in the ICU right now.

Looking for a job or contract: big_polar_bear1@proton.me
