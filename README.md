Note: I am antifa, leader of the cipherpunk, fuck the FBI department. 

#### To run from folder "QS_variant" (Standard SIQS with our number theory as backend):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 200 -base 6000 -debug 1 -lin_size 100_000 -quad_size 1</br></br>
This version is somewhat done, but it doesn't achieve much of an advantage, as the main strength of my work will be the NFS_Variant, which performs NFS with reducible quadratic polynomials<br><br>
#### To run from folder "NFS_variant" (Number Field Sieve with our number theory as backend using reducible quadratic polynomials mod m):</br></br>
To run: python3 polarbearalg_debug.py -key 4387 

I just uploaded polarbearalg_10.py</br>
To run: python3 polarbearalg_10.py -key 4387  (or -keysize 16 to generate a random 16 bit modulus)</br>

The quadratic character base is now succesfully implemented. </br>

For example: </br>

python3 polarbearalg_10.py -key 57599</br>

We get three different types of square relations. The one without wrapping around of the modulus:</br> 

[Debug Info]Adding to left mod N: 4 adding to right mod N: 4 y0: 175410 z: 1 discriminant: 4 qchar: [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0] root: 1 eq1: 18</br>
sqrt_right:  2</br>
sqrt_left:  480</br>
[i]Debug info: disc%N: 4 (disc): 4 jacobi_symbols: [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] prod_left: 230400 prod_right: 4 z: 1 ja: -1 y0: 480</br>
[i]Found a square in the integers, the quadratic character base worked</br>
[SUCCESS]Factors are: 239 and 241</br>

The one with wrapping around of the modulus but where one discriminant is a multiple of another:</br>

[Debug Info]Adding to left mod N: 25125 adding to right mod N: 49984 y0: 95085 z: 1 discriminant: -205271 qchar: [-1, -1, -1, -1, -1, -1, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0] root: 127632 eq1: 0</br>
[Debug Info]Adding to left mod N: 42901 adding to right mod N: 27139 y0: 161610 z: 4 discriminant: -821084 qchar: [-1, -1, -1, -1, -1, -1, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0] root: 223350 eq1: 6</br>
sqrt_right:  410542</br>
sqrt_left:  50250</br>
[i]Debug info: disc%N: 37538 (disc): 168544733764 jacobi_symbols: [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] prod_left: 2525062500 prod_right: 9993601024 z: 25182 ja: 1 y0: 50250</br>
[i]Found a square in the integers, the quadratic character base worked</br>

And then there is another one (the really big one) which doesn't generate a square in the integers but isntead a huge enourmous discriminant and has wrapping around of the modulus.</br>

It does tend to find squares in the integers, and avoiding the case where we find squares where the discriminant is a multiple of the other, is as easy as restricting the quadratic coefficient to negative numbers.</br>

The most important thing is that the quadratic character base works. Its not finding quadratic residues now, it attempts to find squares in the integers. </br>
One appraoch would be number field sieves appraoch, where we build up a huge modulus with the linear algebra step and then take a square root over a finite field.</br>
Another appraoch would be to do something with square multiples of smooths, since the PoC does seem to be succesfully at finding those... I'll need to restructure a few things though.</br>

I think the square multiples of smooths would be the most interesting to test out first, as the PoC demonstrates it can already handle that case, I just need to change a few things for it to actually work with distinct linear coefficients mod N. I'll try to upload a PoC that does that soon. I also already have all the math worked out for that earlier. 

Update: Yea, square multiples of a smooth/discriminant. That would work brilliantly. The PoC is already finding those even with wrapping around of the modulus (just the trivial case for now where the linear coefficient is also a multiple and not distinct, bc I need to change some things, but I know exactly why this happens and what to do) .. so it demonstrates 100% that that would work. I need some time now to refactor my PoC for v11 to achieve this. But I'll release before the end of the week for sure. And then we're done. That's it. What remains then is simplifying and improving everything so we end up with an actual algorithm. I mean... it's already done now I guess.... since I don't need to figure out anything new anymore... all that square multiple shit is stuff I already figured out a few weeks ago. I guess now we just do our victory lap and finish our work. No rush, no sweat. At this point, I know I've won already. And I know people know, because if people don't know, they are the biggest idiots in the world considering what is at stake here.
