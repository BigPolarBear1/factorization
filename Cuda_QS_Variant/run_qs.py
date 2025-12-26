import QSv3_simd
import argparse
import pstats, cProfile
import pyximport
pyximport.install()
keysize=100
key=0
workers=1
debug=0
base=500
lin_size=100_000
quad_size=100
def print_banner():
    print("Polar Bear was here       ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀                       ")
    print("⠀         ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀ ⣀⣀⣀⣤⣤⠶⠾⠟⠛⠛⠛⠛⠷⢶⣤⣄⣀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀   ")
    print("⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣀⣤⣴⠶⠾⠛⠛⠛⠛⠉⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⠙⠛⢻⣿⣟ ⠀⠀⠀⠀      ")
    print("⠀⠀⠀⠀⠀⠀⠀⢀⣤⣤⣶⠶⠶⠛⠋⠉⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠙⠳⣦⣄⠀⠀⠀⠀⠀   ")
    print("⠀⠀⠀⠀⠀⣠⡾⠟⠉⢀⣀⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠹⣿⡆⠀⠀⠀   ")
    print("⠀⠀⠀⣠⣾⠟⠀⠀⠀⠈⢉⣿⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⢿⡀⠀⠀   ")
    print("⢀⣠⡾⠋⠀⢾⣧⡀⠀⠀⠈⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣄⠈⣷⠀⠀   ")
    print("⢿⡟⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣀⠀⢹⡆⣿⡆⠀   ")
    print("⠈⢿⣿⣛⣀⣀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠹⣆⣸⠇⣿⡇⠀   ")
    print("⠀⠀⠉⠉⠙⠛⠛⠓⠶⠶⠿⠿⠿⣯⣄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢠⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢠⣿⠟⠀⣿⡇⠀   ")
    print("⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠻⣦⡀⠀⠀⠀⠀⠀⠀⠀⠠⣦⢠⡄⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣿⡞⠁⠀⠀⣿⡇⠀   ")
    print("⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠘⣿⣶⠄⠀⠀⠀⠀⠀⠀⢸⣿⡇⢸⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⣴⠇⣼⠋⠀⠀⠀⠀⣿⡇⠀   ")
    print("⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⡿⣿⣦⠀⠀⠀⠀⠀⠀⠀⣿⣧⣤⣿⡄⠀⠀⠀⠀⠀⠀⠀⠀⣿⣾⠃⠀⠀⠀⠀⠀⣿⠛⠀   ")
    print("⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣾⠀⠘⢿⣦⣀⠀⠀⠀⠀⠀⠸⣇⠀⠉⢻⡄⠀⠀⠀⠀⠀⠀⡘⣿⢿⣄⣠⠀⠀⠀⠀⠸⣧⡀   ")
    print("⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢠⣿⠀⠀⠀⠙⣿⣿⡄⠀⠀⠀⠀⠹⣆⠀⠀⣿⡀⠀⠀⠀⠀⠀⣿⣿⠀⠙⢿⣇⠀⠀⠀⠀⠘⣷   ")
    print("⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣸⡏⠀⠀⢀⣿⡿⠻⢿⣷⣦⠀⠀⠀⠹⠷⣤⣾⡇⠀⠀⠀⠀⣤⣸⡏⠀⠀⠈⢻⣿⠀⠀⠀⠘⢿   ")
    print("⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⣴⠿⠁⠀⠀⢸⡿⠁⠀⠀⠙⢿⣧⠀⠀⠀⠀⠠⣿⠇⠀⠀⠀⠀⣸⣿⠁⠀⠀⢀⣾⠇⠀⠀⠀⠀⣼   ")
    print("⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣠⡾⡁⠀⠀⠀⠀⣸⡇⠀⠀⠀⠀⠈⠿⣷⣤⣴⡶⠛⡋⠀⠀⠀⠀⢀⣿⡟⠀⠀⣴⠟⠁⠀⣀⣀⣀⣠⡿   ")
    print("⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢠⣿⣿⣿⣤⣾⣧⣤⡿⠁⠀⠀⠀⠀⠀⠀⠀⠈⣿⣀⣾⣁⣴⣏⣠⣴⠟⠉⠀⠀⠀⠻⠶⠛⠛⠛⠛⠋⠉⠀   ")
    return

def parse_args():
    global keysize,key,workers,debug,base,lin_size,quad_size
    parser = argparse.ArgumentParser(description='Factor stuff')
    parser.add_argument('-key',type=int,help='Provide a key instead of generating one') 
    parser.add_argument('-keysize',type=int,help='Generate a key of input size')    
    parser.add_argument('-workers',type=int,help='# of cpu cores to use')
    parser.add_argument('-debug',type=int,help='1 to enable more verbose output')
    parser.add_argument('-base',type=int,help='Size of the factor base')
    parser.add_argument('-lin_size',type=int,help='Size of the factor base')
    parser.add_argument('-quad_size',type=int,help='Size of the factor base')
    args = parser.parse_args()
    if args.keysize != None:    
        keysize = args.keysize
    if args.key != None:    
        key=args.key
    if args.workers != None:  
        workers=args.workers
    if args.debug != None:
        debug=args.debug  
    if args.base != None:
        base=args.base  
    if args.lin_size != None:
        lin_size=args.lin_size  
    if args.quad_size != None:
        quad_size=args.quad_size   
    return

if __name__ == "__main__":
    parse_args()
    print_banner()
    #cProfile.runctx("QSv3_simd.main(keysize,workers,debug,base,key,lin_size,quad_size)", globals(), locals(), "Profile.prof")

    #s = pstats.Stats("Profile.prof")
    #s.strip_dirs().sort_stats("time").print_stats()
    QSv3_simd.main(keysize,workers,debug,base,key,lin_size,quad_size)