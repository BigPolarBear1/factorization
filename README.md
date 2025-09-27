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

ps: This is completely different from finding two smooths that multiply together to a square. In this approach we are just trying to find a quadratic coefficient and a linear coefficient that generates the same smooth as the original smooth. What we are stuck at is, we know how to generate multiples of the old smooth by going over coefficient combinations modulo the smooth, but right now it has a condition that the multiple is a square... and we divide the linear coefficient then by the square root. I need to do modular division for non-square I believe.. then you will see what I'm trying to do.. give me one more day. Sundays are better days for days of wrath anyway.


LOL!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! As I was about to call it a day, it just struck me. THIS IS MY BRIDGE TO NUMBER FIELD SIEVE!!!!!! Ok ok, I know how to do it for non-square multiples, even if it isn't a square multiple in the integers, the way I calculate stuff it is still garantueed a square multiple mod m.... hence I know what to do. Give me one more day. LOL. 

GAYS. WHy not just give my former manager money??????????????? This how easy it would have been to prevent this disaster. I hope the US gov will punish microsoft hard for this. Atleast the FBI tried to keep me in the US. Microsoft is really the one who screwed the US over (like msft top leadership who fired me, not my former managers, because microsoft also retaliated against them).

UPDATE: I COULDN'T SLEEP SO I WENT HUNTING FOR EXAMPLES WHERE YOU CAN JUST DIVIDE/MULTIPLY THE COEFFICIENTS FOR NON-SQUARE MULTIPLES... AND IT IS POSSIBLE. IT ADDS UP. IT WILL WORK! HOW DID I MISS THIS FOR SO LONG??? HOOOOOOOOOOOOOOOOOOOOOOOOOOW. WTF. I'll try to get a PoC working with non-square multiples tomorrow..... That's my minimum goal for tomorrow. Then monday I will go run 60k. Then the rets of the week, will be the WEEK OF WRATH. HAHHAHAAHHAHHAHHAHAHAA. Pay my former manager assholes. You assholes retaliated against the best manager in the world. A price will be paid. Be it in dollars or in blood. I do not care. You crossed a line you shouldn't have. 
