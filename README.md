DISCLAIMER: AT NO POINT IN MY 2.5 YEARS OF RESEARCH DID I USE AI OR HAS AI IN ANY WAY HELPED WITH MY RESEARCH. I'M PUTTING THIS HERE BECAUSE I KNOW PEOPLE WILL ACCUSE ME OF IT. 

Note: I am antifa, leader of the cipherpunk, fuck the FBI department. 

#### To run from folder "QS_variant" (Standard SIQS with our number theory as backend):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 200 -base 6000 -debug 1 -lin_size 100_000 -quad_size 1</br></br>

See below for an improved way of performing what this PoC does.. I'll delete this Proof of Concept once the PoC for Improved_QS_Variant matures a little more<br><br>
#### To run from folder "NFS_Variant" (Sort of a QS/NFS hybrid that improves on Quadratic Sieve using the number theory I developped):</br></br>

I just uploaded polarbearalg_v44.py</br>
To run: python3 polarbearalg_v44.py -keysize 50</br>

I removed the jacobi symbols again. I will experiment with them later. But I need to change some things. I do remember a few days ago, when setting y1 to be the square we found and generating y1^2-N*y1... then then other square mod n would be some multiple of that.. and that multiple.. there was something going on with legendre symbols there... I should dig deeper into that.

First I want to try some heuristic thats going to find an optimal amount to add the modulus to x and z. We need some x and z combination that:

1. Keeps z as small as possible
2. and zx^2-n should also be as small as possible. 

Let me bash my head against that first today... because solving that is the path of least resistance... then later I can explore some more complicated stuff again.
