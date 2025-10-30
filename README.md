DISCLAIMER: AT NO POINT IN MY 2.5 YEARS OF RESEARCH DID I USE AI OR HAS AI IN ANY WAY HELPED WITH MY RESEARCH. I'M PUTTING THIS HERE BECAUSE I KNOW PEOPLE WILL ACCUSE ME OF IT. 

Note: I am antifa, leader of the cipherpunk, fuck the FBI department. 

Note2: It's basically done now, the NFS variant is nearing completion soon, nothing novel research wise has to be figured out now. I just need to refactor a few things and I'm done.
I have been unemployed for 2 years. I am looking for work. I would like some way where I can have enough time to attack ECC and DLP next. In addition, for a payment, I can also show people how I appraoched factorization from start to finish. Because the method of how I conduct my research and consistently have success is just as valueable as the research itself. And I follow the same patterns over and over again, no matter what subject matter I attack. I don't have a very high IQ and I'm not especially bright or good at doing numbers in my head. But I have insights into research that help me succeed over and over again. I'm willing to teach people how to do it, as long as it enables me to have the finances to continue doing my research aside from that. Contact: big_polar_bear1@proton.me.

#### To run from folder "QS_variant" (Standard SIQS with our number theory as backend):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 200 -base 6000 -debug 1 -lin_size 100_000 -quad_size 1</br></br>
This version is somewhat done, but it doesn't achieve much of an advantage, as the main strength of my work will be the NFS_Variant, which performs NFS with reducible quadratic polynomials<br><br>
#### To run from folder "NFS_variant" (Number Field Sieve with our number theory as backend using reducible quadratic polynomials mod m):</br></br>
To run: python3 polarbearalg_debug.py -key 4387 (for use with the final chapter of the paper... but ignore this and use PoC below until I refactor the paper one final time)

I just uploaded polarbearalg_v23.py</br>
To run: python3 polarbearalg_v23.py -keysize 14</br>

Made some further improvements compared to v22. Now instead of trying coefficients for a large modulus. We test coefficients mod p^2 (small square moduli with a single prime factor)... and when certain conditions are met, we save them and multiply them together at the end. This is simply an intermediate step I'm doing to make the transition to linear algebra... I'm going to spent a little more time on v23 and see how a quadratic character base would work, so we only multiple coefficients mod p^2 together that are going to garantuee a square relation mod n. I want to have that figured out first, before I add trail division and linear algebra.



UPDATE: I GOT IT! YOU CAN MOVE ANY FACTOR FROM THE LINEAR COEFFICIENT TO THE QUADRATIC COEFFICIENT, HENCE REDUCING THE SIZE OF THE LINEAR COEFFICIENT. I shouldn't use x^2+yx-zN, I should only use zx^2+yx-N. So that no matter what quadratic coefficient we have... that polynomial will always add or subtract N one single time. THIS IS IMPORTANT. Because if we know N is going to be added or subtract only once... then we can modulo redue N and we don't need to worry about multiples of N. I SEE IT IN MY HEAD!!! I SEE THE SOLUTION! EUREKA!
