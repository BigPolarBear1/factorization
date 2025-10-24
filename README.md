Note: I am antifa, leader of the cipherpunk, fuck the FBI department. 

Note2: It's basically done now, the NFS variant is nearing completion soon, nothing novel research wise has to be figured out now. I just need to refactor a few things and I'm done.
I have been unemployed for 2 years. I am looking for work. I would like some way where I can have enough time to attack ECC and DLP next. In addition, for a payment, I can also show people how I appraoched factorization from start to finish. Because the method of how I conduct my research and consistently have success is just as valueable as the research itself. And I follow the same patterns over and over again, no matter what subject matter I attack. I don't have a very high IQ and I'm not especially bright or good at doing numbers in my head. But I have insights into research that help me succeed over and over again. I'm willing to teach people how to do it, as long as it enables me to have the finances to continue doing my research aside from that. Contact: big_polar_bear1@proton.me.

#### To run from folder "QS_variant" (Standard SIQS with our number theory as backend):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 200 -base 6000 -debug 1 -lin_size 100_000 -quad_size 1</br></br>
This version is somewhat done, but it doesn't achieve much of an advantage, as the main strength of my work will be the NFS_Variant, which performs NFS with reducible quadratic polynomials<br><br>
#### To run from folder "NFS_variant" (Number Field Sieve with our number theory as backend using reducible quadratic polynomials mod m):</br></br>
To run: python3 polarbearalg_debug.py -key 4387 

I just uploaded polarbearalg_v13.py</br>
To run: python3 polarbearalg_v13.py -keysize 14</br>

I changed how things are represented. Left side of the congruence is now fully in the  integers. Right side we modulo reduce y0**2 before calculating the discriminant. Ideally I want to modulo reduce that entire discriminant. But one step at a time. The linear algebra seems to work as intended now. And so do the jacobi symbols. We are consistently achieving factorization now. The squares it finds arn't congruent mod N, but still yield the correct factorization... you need the adjust the coefficients on one side... like I did with my QSv3_find_similar.py work last month... I'll implement that tomorrow first thing. After that I'll try to figure out how to modulo reduce that entire discriminant. 
Inched again a little closer now to the final solution. I know I've won already. Everything else now is trivial and Ill figure it out in mere days now. I'm just doing my victory lap now. 

Update: uploaded v13. All this does compared to v12 is remove the left side from the matrix. Since we are having that to be squares in the integers, there is no reason to perform linear algebra on it. I've also worked out what is happening and how to adjust our results so we get a square that is congruent mod N instead... I'll add that for v14.. either tonight or hopefully tomorrow.

Update: Thinking about more about this. This is exactly the analogue with NFS that I have been searching for for atleast half a year now. So in essense, what I need to do now, is what is in NFS "taking a square root over a finite field" ... that way we generate squares congruent mod N and succeed at taking the GCD more often. Once that is done, I have a pretty good idea how the final algorithm is supposed to come together.
