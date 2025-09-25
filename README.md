# factorization

The PoC under PoC_Files will easily factor above 200 bits. </br></br>

However, that one is basically just using the number theory from my paper to achieve default SIQS. </br></br>

in PoC_Files_Find_Similar_WIP I am attempting to sieve in the direction of the quadratic coefficient and linear coefficient. This way in theory we need much less smooths.</br></br>

For PoC_Files_Find_Similar_WIP:</br>
Build: python3 setup.py build_ext --inplace</br>
Run: python3 run_qs.py -keysize 25 -base 30 -debug 1 -lin_size 2 -quad_size 1</br></br>

You will see some output like this:</br></br>

Looking for similar smooths, amount needed: 7 initial smooth factors: {3, 101, 7, 23, 59, -1}</br>
Found a similar smooth{3, 7, 23, -1} co: 1738 quad: 1</br>
Found a similar smooth{3, -1} co: 49806 quad: 57</br>
Looking for similar smooths, amount needed: 7 initial smooth factors: {3, 7, 43, 79, 23, 59}</br>
Found a similar smooth{43, 3} co: 25126 quad: 1</br>
Found a similar smooth{43, 3} co: 75378 quad: 9</br>
Found a similar smooth{3, -1} co: 6593 quad: 1</br>
Found a similar smooth{3, -1} co: 19779 quad: 9</br>
Found a similar smooth{43, 79} co: 42083 quad: 1</br>
Found a similar smooth{43, 79} co: 126249 quad: 9</br>
Found a similar smooth{43, 79} co: 294581 quad: 49</br>
sim_found:  7</br>
Running linear algebra step with #smooths:  11</br>
[SUCCESS]Factors are: 3793 and 3541</br>

This will try to find smooths with a similar factorization.
Right now it just p-adicaly lifts linear coefficients and then checks them without doing any Chinese Remainders (combinding coefficients). 
If a smooth candidate is bigger then 0, then it will add the modulus to the quadratic coefficient to generate a smaller smooth.

So while I was adding the modulus to the quadratic coefficient, I also realized, that once the smooth becomes a big negative number, I can then again add the modulus to the linear coefficient.
I wonder how I can find that sweet spot where we get a  super small smooth by adding the modulus to the linear/quadratic coefficients.
I'll spent some time on that... because that may be an easy way to find smooths. 

Then next, we need some simple heuristic to chinese remainder coefficients together... but only those that are going to garantuee a smooth. Just want to go after those easy to find smooths. 

Another thing I need to triple, quadruple check is this:</br>

Found a similar smooth{43, 79} co: 42083 quad: 1 </br>
Found a similar smooth{43, 79} co: 126249 quad: 9 </br>
Found a similar smooth{43, 79} co: 294581 quad: 49 </br>

These are 3 smooths that have the same factorization.

42083^2-13431013\*4\*1 =  1717254837 (3 x 3 x 3 x 3 x 43 x 79 x 79 x 79)</br>
126249^2-13431013\*4\*9 = 15455293533 (3 x 3 x 3 x 3 x 3 x 3 x 43 x 79 x 79 x 79) </br>
294581^2-13431013\*4\*49 = 84145487013 3 x 3 x 3 x 3 x 7 x 7 x 43 x 79 x 79 x 79 </br>

This is just what happens if we multiple the linear coefficient, 42083, by 3 and the quadratic coefficient by 9, or by 7 and 49.
This is an easy way to generate smooths instantly from a starting smooth. 
But all of those should be congruent mod N. These should be bad smooths.
These alone should never, never, never, yield the factorization of N. 
Or would they? OR WOULD THEY? I need to double check this. Just in the odd case I'm dismissing this as something that won't work, while it will. Because generating these types of smooths by multiplying coefficients, it's instantly. It lets me create smooths out of thin air. I need to make sure there really is never ever any use for these types of smooth... fucking brain fog. My intuition tells me these are just bogus smooths, just permutations of the original smooths. But in the above example, it is mainly just bogus smooths and still succeed when taking the GCD. So I need to spent some good amount of time on that, just in case my intuition is wrong. 



