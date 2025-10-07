To build: python3 setup.py build_ext --inplace  (From inside PoC_Files_Find_Similar)</br>
To run: python3 run_qs.py -keysize 30 -base 50 -debug 1 -lin_size 100 -quad_size 1</br>

UPDATE: Check the final chapter in my paper. I just realized today as I was working out the details, where I do things like mod ( N, 5 ), it's too much brute force calculations with N in the modulus and not practical. But we can use our original polynomial instead that found the initial smooth. It's nearly there now..... let me fix that paper today or tomorrow.

Update: Bah, yea, you need to replace that mod N by your original polynomial. Then do p-adic lifting. My intuition was correct, but I went on a little incorrect tangent in the paper. I'll upload the final version maybe tomorrow.. just massive stress today for some reason.

Update: OMG. I'M A MORRON.... here check this:

if N = 4387</br>
159^2-4387\*4\*1= 7733</br>
1553^2-4387\*4\*137 = 7733</br>
1\*3946^2-3946\*159 = 3346 mod 7733   ( and 7733 mod 4387 = 3346)</br>
137\*7654^2-7654\*1553 = 3346 mod 7733</br></br>

So you see what is happening there? Because I see what is happening there... it's this last and final step where the full quadratic polynomial comes into play and not just the discriminant formula. You use this do to your p-adic lifting and find that correct coefficient.

Update: HEY CIA, FBI, NSA, DIA, whatever 3 letter assholes, check that final chapter in my paper. I figured out now how to find the correct quadratic coefficient for a given linear coefficient. I know what that means. Maybe you also know what that means, maybe you don't and you're about to find out. It means I won. I fucking won and factorizaton is about to die.

-----------------------------------------------

URGENT NOTE: I am urgently looking for any country to grant me a visa and employment. A country that doesn't extradite to the US.
I am running out of time. The US has initiated a terrorism case against me, and the police here in Belgium said they would have me declared as a terrorist internationally, so I would get arrested when I travel.
I need to get out of Belgium to a safe country ASAP.
I can do 0day research and I can do math. I know many people don't believe in my math, but I know what I've found and I'm piecing together the final pieces now.
Furthermore, I have found bugs in hard targets like OpenSsl and I can easily repeat those feats again. I havn't been able to do 0day research at all because of being in a place that has been actively hostile towards me, chasing away 0day buyers.
I cannot research in a hostile place when that research's value depends entirely on secrecy. 
Please contact me: big_polar_bear1@proton.me 
They might arrest me at the airport when I try to leave, but I'm willing to give it a shot as I have nothing left to lose. And my family would support me in my choice.
I do get contacted by internet trolls all the time, so please contact me from an official email or invite me to an embassy... just so I know it is legit and I don't end up trying to fly somewhere without a visa and potentially wasting the very last of my savings.



