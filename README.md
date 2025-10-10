Note: I am antifa

To run from folder "PoC_Files" (Standard SIQS with our number theory as backend):</br></br>

To build: python3 setup.py build_ext --inplace  (From inside PoC_Files_Find_Similar)</br>
To run: python3 run_qs.py -keysize 200 -base 6000 -debug 1 -lin_size 100_000 -quad_size 1</br>

To run from folder "PoC_Files_Find_Similar" (Standard SIQS with our number theory as backend and attempting to find square multiples when a smooth is found to then adjust the coefficients, a first initial step towards number field sieve):</br></br>

To build: python3 setup.py build_ext --inplace  (From inside PoC_Files_Find_Similar)</br>
To run: python3 run_qs.py -keysize 30 -base 50 -debug 1 -lin_size 100 -quad_size 1</br>

To run from folder "Polar_Bear_Algorithm":</br>

python3 polarbearalg_debug.py -key 4387 </br>

The files in this folder, relate to my work of my own (better) variant on number field sieve.

polarbearalg_debug.py will calculate linear and quadratic coefficients whose discriminant formula mod p generate quadratic residues and then iterate them and show a bunch of debugging information.
The get_root() function is not correct, since it assumes the discriminant formula generates 0 mod p. Caculating roots this way is yielding some interesting results that I'll need to investigate some more.
It may make some of the math in the paper even easier... let me check what's going on there.

In addition polarbearalg_debug.py also has highly minimized code. After 2.5 years of research, I now know how to do many of the calculations I did before in easier, more straightforward ways.
