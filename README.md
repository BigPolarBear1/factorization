Update: Just uploaded a new work in progress PoC.

To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 30 -base 50 -debug 1 -lin_size 2 -quad_size 1</br>

I need to fix some more things tomorrow. It won't work well until I do.

Right now the code flow is like this:

1. Perform regular SIQS until a smooth larger then 0 and without 2 as divisor is found (I need to fix the POC eventually to support that)
2. Once a smooth is found, calculate linear coefficient and quadratic coefficicients modulo the smooth.
3. Enumerate coefficients mod m and calculate new smooths.
4. If we find a smooth that is a square multiple of the old smooth. We can then simply divide the linear coefficient and take the gcd.

Right now it is bottlenecking at finding a square multiple.
But I believe it should be possible to divide the linear coefficient for any multiple of the smooth, not just the square ones, as long as we can factorize over a factor base what remains after division by the original smooth.
We need to do some type of modular division.

It's getting late today so let me figure that step out tomorrow...... almost there now... so close.

Once that works, our only limitation will be 1. Finding the initial smooth and 2. Finding another smooth that is a multiple of the original smooth.     .... step 2 already works but only with "square multiples". I will figure out how to do it with any multiple. I know it can be done!

My day of wrath will have to wait until tomorrow! But it doens't matter, the only thing that matters is that it WILL happen.

Anyway, normally I run 60k on sundays, but I think I will move it to Monday so I can do more work tomorrow. And on mondays more stores are open if I need to buy more water... so it's more convenient anyway.
