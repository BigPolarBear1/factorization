DISCLAIMER: AT NO POINT IN MY 2.5 YEARS OF RESEARCH DID I USE AI OR HAS AI IN ANY WAY HELPED WITH MY RESEARCH. I'M PUTTING THIS HERE BECAUSE I KNOW PEOPLE WILL ACCUSE ME OF IT. 

Note: I am antifa, leader of the cipherpunk, fuck the FBI department. 

#### To run from folder "QS_variant" (Standard SIQS with our number theory as backend):</br></br>
To build: python3 setup.py build_ext --inplace</br>
To run: python3 run_qs.py -keysize 200 -base 6000 -debug 1 -lin_size 100_000 -quad_size 1</br></br>

See below for an improved way of performing what this PoC does.. I'll delete this Proof of Concept once the PoC for Improved_QS_Variant matures a little more<br><br>
#### To run from folder "Improved_QS_Variant" (Sort of a QS/NFS hybrid that improves on Quadratic Sieve using the number theory I developped):</br></br>

I just uploaded polarbearalg_v37.py</br>
To run: python3 polarbearalg_v37.py -keysize 50</br>

UPDATE: I NEED TO FIX ONE MORE IMPORTANT STEP IN THE POC TOMORROW. WE DONT NEED zx^2 TO BE A SQUARE IN THE INTEGERS. WE CAN TAKE A SQUARE ROOT OVER A FINITE FIELD THERE INSTEAD. PPFFFFT... sure took me long enough to connect the dots there. I know how to do it.. just worked out the math. I need some sleep and I'll upload the fixed POC tomorrow.
