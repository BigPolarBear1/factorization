Disclaimer: No AI was used for any of this, except for reviewing my paper these last 2 weeks, but AI has not written a single sentence.
None of the code is written by AI either, except for one or two functions like lift_root2(),
which is just hensel, something I had already implemented before but with a coefficient list as input.
I don't believe AI is quite there yet to do math research. It's very rigid and can't think outside the bounds of existing literature
and often just makes very dumb conceptual mistakes and it has a total lack of creating abstractions.
It is a tool, well suited for basic tasks, nothing more.
This research is also still ongoing, and especially some of the things stated in the last chapter might be missing the mark. 
I'll also be properly learning about number fields now, and see how all of that can be fit into my work.

#### References (re-used many of the core number theoretical functions from these PoCs to fit my own algorithm): 
https://stackoverflow.com/questions/79330304/optimizing-sieving-code-in-the-self-initializing-quadratic-sieve-for-pypy
https://github.com/basilegithub/General-number-field-sieve-Python 

#### About the paper
Math paper is a work in progress. The final chapters are a bit rushed and building an algorithm around p-adic lifting isnt as straight forward as I had assumed. I do think there is an angle there I can exploit, but I'll do some further experimentation first and get a working PoC before I make edits to the paper again.

#### To run from folder "polysieve" WORK IN PROGRES...extremely early version:</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run:    python3 run_qs.py -keysize 100 -base 500 -debug 0 -lin_size 10_000 -quad_size 1

Actually, with my two-sided setup using binomial expansions, perhaps both sides can have numberfields rather then pinning one side in the rationals (edit: intuition was right, but that's not really how numberfields work, see update below). I'll learn what I have to in Iceland and start coding when I'm back. Nature helps to think about this. I'll update this with a numberfield implementation shortly.. ignore this uploaded version for now.

Update: Started digging into abstract algebra and number fields. So from what I understand so far, the fact that kleinjung uses a polynomial f(x) of degree d and linear polynomial g(x), where the resultant = N, is less then optimal theoretically. As having f(x) and g(x) be the same degree would be better. I want to start bashing this possible improvement using my binomial expansion work. Things are quickly starting to click inside my head. But I'm still early in my abstract algebra journey so this will take some time..

Note: My next update to this repo will likely be some implementation of this, if I can get it to work. On a side note, I am still looking for work: big_polar_bear1@proton.me. I have no problems taking my research private or pivot back into VR if it implies salaried work, or some type of hybrid. I will climb this math hill, and I will eventually succeed at some type of breakthrough, minor or major. I don't think people fully understand how determined I am and currently, everything happening in life is just pushing me deeper into this math obsession. Although I am sure some elements are hoping I go into the woods to **** myself out of hopelessness.

Update: Actually, I may have an initial version up and running next week already in sage. Because I really just need to start with binomial expansion to create f(x) and g(x) of say, both degree 3, optimize the coefficients, use the residue math and stuff to make sure it has a high murphy-E score and run NFS over it. That I can do fairly quickly with the help of sage. Then what remains is research and getting deeper into the literature.. easy enough. Alright.. lets go.

Update: Wont upload anything for a while. Going to stop uploading experimentations until I have something that is a prove-able breakthrough with code. My offer to take my work private also still stands.

#### To run from folder "Coefficient_Sieve" (For use with the paper):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run:  python3 run_qs.py -keysize 40 -base 50 -debug 1 -lin_size 10_000 -quad_size 100</br>

Just demonstrates the math from the paper using quadratics. For educational purposes. And rather then taking a square root over a large prime we can also just calculate the discriminant. But this demonstrates the interesting relation between these quadratics and the factors of N.

#### To run from folder "CUDA_QS_variant" (Failed Experiment):</br></br>
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

This basically creates a system of quadratics. Solving them mod p is easy. But there is only one root solution (the factor of N) which solves the system for 0 for any mod p (aka solves it in the integers). Figuring out how to exactly do this quickly is still an ongoing area of research for me. And if a polynomial time algorithm for factorization exists, it is likely done by solving this system of quadratics. 


