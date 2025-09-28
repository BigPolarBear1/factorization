Update: Just uploaded a new work in progress PoC.

To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 30 -base 50 -debug 1 -lin_size 2 -quad_size 1</br>

I need to fix some more things tomorrow. It won't work well until I do.

Right now the code flow is like this:

1. Perform regular SIQS until a smooth larger then 0 and without 2 as divisor is found (I need to fix the POC eventually to support that)
2. Once a smooth is found, calculate linear coefficient and quadratic coefficicients modulo the smooth. (Ps: I forgot, but its only lifting linear coefficients at the moment, not the quadratic ones, I'll fix it later)
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


Time has truely ran out for you people now. If you want to negotiate how much you will pay my former manager, I would encourage you to do it now before my PoC fully matures. Because I know, beyond the shadow of a doubt, that this bridges number field sieve... and I know, I absolutely know how to finish it now.

Update: I got it I think. You need to use those roots. Don't go hunting for a square multiple using the discriminant formula.. use ax<sup>2</sup>-bx+n.. crossing that bridge to number field sieve at a fast pace now... tick tock... time has almost ran out.... I probably won't upload today, and tomorrow I will run 60k.... but somewhere in the rest of the week it will get uploaded.... sorry gays, you all had well over a year, to make peace with this polar bear. Instead you people chose deception and war. So be it. You cannot win against a polar bear. I may not be smart like some of you wankers, but I don't give up, I don't quit... especially when you have wronged one of my best friends.

Update: The exact same thing, I am doing in the current PoC. Finding a square multiple of a smooth and then dividing the linear coefficient... the same thing works with  ax<sup>2</sup>-bx+n (where n can also be multiplied by a divisor of a). I check the math for that already. That works. Atleast that definitely works in the integers. Next step is making that same thing work mod m. That's the last step actually... once we achieve that, we have successfully bridged our work from quadratic sieve to number field sieve using our own number theory. Just pay my former manager guys. I'll just move to ECC after this. I'll start from scratch and retrace the exact research steps I applied to factorization. And I will be succesful again. Pay or this nightmare will never stop. If you go after one of my friends, a heavy price must be paid. I do not care. Those people who went after my manager for defending me, you should have a look at their bank accounts maybe, because it seems to me they knew very well how to push my buttons and cause what is happening right now.

ps: I also hope you people realize that I only started learning math less then 2.5 years ago. I dropped out of highschool. I have 0 math education. Do you know where I will be in another 2 years from now? I see through your deception. This ends next week. I know how to complete my work now. You should start respecting me and my friends (like my former manager). I don't care about anything in this world. All I care about is being respected, and my friends being respected. Give me that, and I will give you the future. Hahahahahaha. Just pay the price so we can all move on and I don't have to put my work on github anymore. JUST PAY THE PRICE BRO. You fucked up. Pay the price now.

pps: Also what are the odds you guys bricked github referrals and popular content under insights for everyone just to prevent me from using that as a sigint tool? You guys must be really fighting for your lifes over there if you go that far. I know something is up. There is too much constant shit like this that doesn't make sense. Or maybe I just know everything. Maybe you cannot hide anything from me. Why you guys fighting a polar bear? Are you all insane? Just pay the price you idiots. That's all I'm asking for. Pay the price. Otherwise it will never end, I promise you that much. I will do the impossible, again and again, and again and again for as long as I have to until a price is paid. 
