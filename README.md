Disclaimer: At no point did AI contribute anything to this research project. Copilot can't even calculate basic congruences without making mistakes. People who think AI can do novel math research are delusional. I did not use AI for the paper, and I did not use AI for any of the code. I attacked this problem, from the very basics 3 years ago, slowly burrowing deeper and deeper. 

The day I break factorization will be day 0 of the gay future. A better world will be born. A world of human creativity. A world of art, science and spirituality. A world without suffering. I will see that world become reality, no matter the cost. 

#### To run from folder "NFS_WIP2" (Experimental WORK IN PROGRESS):</br></br>

##NFS related code is borrowed from: https://github.com/basilegithub/General-number-field-sieve-Python (note: Very impressively written, helped me big time, thanks)

To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 50 -base 60 -debug 1 -lin_size 10_000_000 -quad_size 20</br></br>

To do: The next thing to implement now is to use a bigger initial zx+y value instead of just 1. This basically divides the distances between the factors (while the polynomial value acts as an offset) and makes it so we can easier find those case where the factors are far apart without increasing the sieve interval size. 

Once that is done, we also need to re-introduce using different quadratic coefficients. Because I noticed that using divisors of k (as in zx^2+yx-Nk) and moving them over to the quadratic coefficient instead, like I demonstrate in CUDA_QS_variant.. this actually changes the Legendre symbols. Which makes sense with what we see in debug.py for example. So deciding which factors should be part of the quadratic coefficient and which ones should be part of k, I'm hoping to solve that with linear algebra. That is the big ticket item that should result in a breakthrough algorithm.

And then there is also a bunch of smaller things I need to fix. I.e in initialize_3() we dont need to calculate those roots because we already have those in our precalculated hashmap.

Update: Going to take a break until Sunday. Using a zx+y value that is larger then 1 is quite important. We really want that to be reasonably large. But implementing it is a bit tedious. If I want zx+y to be square and the polynomial value to be square.. then to get it to work with a sieve interval I'll have to lift my hashmap to an even power. Its going to be quite some tedious number theory.. definitely possible.. but tedious and with an ultra marathon coming up, I want to focus on that rather then implement some tedious math. I am curious though.. once it is implemented what performance that will yield since it reduce the size the sieve interval will need to be. And ofcourse... all of this is still without linear algebra step. But I want to implement everything else first before I attack that problem. I have a pretty good idea how to do it though.

Update: And to clarify.. the reason I'm waiting to implement any form of linear algebra until everything else is done is bc it may not even be needed. I like the idea of combining smooths with linear algebra.. but there is a chance that this ends up being faster, lets see. Its just a matter of using that precalculated hashmap and getting all my legendre symbols to be square.

Update: Depression man. I know what I got. I know it is a matter of time. I started with 0 math skills 3 years ago. Never finished highschool. Couldn't even do basic algebra. I understand factorization in an unique way. I know I do. And I know what the consequences will be when I succeed. Even just threatening RSA1024 would be a disaster, because it is probably still used in a lot of places. Especially industrial shit with memory constrained hardware. The irony is, the only reason I keep working on this, aside from it being a response to traumatic experiences, is that I have no other choices. Nobody wants to buy my 0day exploits. People are being chased away. What am going to do? Yes the Chinese and Russians would buy my work, but everything that is happening, it indicates I'm probably under surveillance (and I dont care how schizo that sounds, there is just too many indicators that dont make any sense and have all my alarm bells going off, I'm not an idiot and I'm not blind). So whats the point of even selling to the Russians or Chinese? They would just have microsoft patch my work immediatly or use SIGINT to gain insight into whoever is buying my work, I cant protect myself against those things while stuck in this Belgian shithole (yes it is a shithole, my entire life I could never find employment and not for a lack of trying). And honestly, the way I'm treated.. I probably would sell to the Chinese if I could. For LGBTQ people, everywhere in the world is shit right now, but atleast the Chinese respect me for my talents.. unlike whatever the fuck is going on here in the west. It is starting to feel like people in the west are trying to make me as desperate as possible and push me into committing suicide.. like those fucking stalker creeps kept emailing me for years. Wouldnt be suprised if those stalkers where part of the same group telling 0day buyers not to buy my work lol... shithead western infosec techbro losers with 0 talent. 

An intermediate step between QS and NFS (representing chapter VII in the paper) can be found here: https://github.com/BigPolarBear1/factorization/tree/7deba681fc78c349ea70e514a36ab723399f8e96/NFS_Variant_simple

#### To run from folder "CUDA_QS_variant":</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run: To run:  python3 run_qs.py -keysize 240 -base 100_000 -debug 1 -lin_size 100_000_000 -quad_size 100</br></br>
 
Prerequisites: </br>
-Python (tested on 3.13)</br>
-Numpy (tested on 1.26.2)</br>
-Sympy</br>
-cupy-cuda13x</br>
-cython</br>
-setuptools</br>
-h5py</br>
(please open an issues here if something doesn't work)</br></br>

Additionally cuda support must be enabled. I did this on wsl2 (easy to setup), since it gets a lot harder to access the GPU on a virtual machine.

This was an attempt at finding smooths with similar factorization using an SIQS variant. By using quadratic coefficients. But it didnt end up working as I had hoped so I abondoned this approach, but perhaps someone will get some use out of it.

#### To run debug.py" (Prints the linear and quadratic coefficients to solve for 0 in the integers, for use with my paper):</br></br>

To run: python3 debug.py -keysize 12

This basically creates a system of quadratics. Solving them mod p is easy. But there is only one root solution (the factor of N) which solves the system for 0 for any mod p (aka solves it in the integers). Figuring out how to exactly do this quickly is still an ongoing area of research for me. And if a polynomial time algorithm for factorization exists, it is likely done by solving this system of quadratics. Finding a polyomial time algorithm is my ultimate goal, as this would make progress toward solving p = np as well. 

