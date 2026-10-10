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
mode="psieve"
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
    global keysize,key,workers,debug,base,lin_size,quad_size,mode
    parser = argparse.ArgumentParser(description='Factor stuff')
    parser.add_argument('-key',type=int,help='Provide a key instead of generating one') 
    parser.add_argument('-keysize',type=int,help='Generate a key of input size')    
    parser.add_argument('-workers',type=int,help='# of cpu cores to use')
    parser.add_argument('-debug',type=int,help='1 to enable more verbose output')
    parser.add_argument('-base',type=int,help='Size of the factor base')
    parser.add_argument('-lin_size',type=int,help='Size of the factor base')
    parser.add_argument('-quad_size',type=int,help='Size of the factor base')
    parser.add_argument('-mode',type=str,choices=['psieve','nfs'],help='psieve (default): SIQS + psieve.  nfs: SIQS and a degree-d number field sieve (gnfs.pyx) in separate processes feeding one matrix, the NFS aimed at the singleton factors')
    parser.add_argument('-nfs_degree',type=int,help='nfs mode: degree of the number field sieve polynomial (default 3)')
    parser.add_argument('-nfs_lines',type=int,help='nfs mode: lines added to the ordinary sieve region by each NFS call (default 100)')
    parser.add_argument('-nfs_t',type=int,help='nfs mode: half-width of the ordinary region (default: skew of the polynomial * lines so far)')
    parser.add_argument('-nfs_want',type=int,help='nfs mode: an NFS call stops its ordinary region once it has made this many b-smooths (default 600)')
    parser.add_argument('-nfs_lp',type=int,help='nfs mode: large prime bound in bits; a relation may keep one prime up to 2^N outside the factor base (default 24, 0 = off)')
    parser.add_argument('-nfs_base',type=int,help='nfs mode: algebraic factor base of the NFS = the first N primes of the SIQS factor base (default: all of it)')
    parser.add_argument('-la_every',type=int,help='nfs mode: try the matrix each time this many new b-smooths have arrived (default 1000)')
    parser.add_argument('-nfs_force_t',type=int,help='nfs mode: half-width of the strip per line on the forced lattice (default: the largest factor base prime)')
    parser.add_argument('-nfs_force_lines',type=int,help='nfs mode: lines to sieve on the forced lattice per call (default 1000)')
    parser.add_argument('-nfs_sing_qr',type=int,help='nfs mode: above -nfs_sing_small, primes up to this bound join the rational factor base only if n is a square modulo them, the only primes a SIQS relation can hold (default -1 = up to the end of the SIQS factor base, 0 = off)')
    parser.add_argument('-nfs_sing_small',type=int,help='nfs mode: every prime up to this bound counts as core: it joins the rational factor base and never makes a relation a singleton when targets are chosen (default 200)')
    parser.add_argument('-nfs_sing_force',type=int,help='nfs mode: number of singleton targets to force per run (default 1)')
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
    if args.mode != None:
        mode=args.mode
    if args.nfs_degree != None:
        QSv3_simd.NFS_DEGREE=args.nfs_degree
    if args.nfs_lines != None:
        QSv3_simd.NFS_LINES=args.nfs_lines
    if args.nfs_t != None:
        QSv3_simd.NFS_T=args.nfs_t
    if args.nfs_want != None:
        QSv3_simd.NFS_MIX_WANT=args.nfs_want
    if args.nfs_lp != None:
        QSv3_simd.NFS_LP_BITS=args.nfs_lp
    if args.nfs_base != None:
        QSv3_simd.NFS_BASE=args.nfs_base
    if args.la_every != None:
        QSv3_simd.LA_EVERY=args.la_every
    if args.nfs_force_t != None:
        QSv3_simd.NFS_FORCE_T=args.nfs_force_t
    if args.nfs_force_lines != None:
        QSv3_simd.NFS_FORCE_LINES=args.nfs_force_lines
    if args.nfs_sing_small != None:
        QSv3_simd.NFS_SING_SMALL=args.nfs_sing_small
    if args.nfs_sing_qr != None:
        QSv3_simd.NFS_SING_QR=args.nfs_sing_qr
    if args.nfs_sing_force != None:
        QSv3_simd.NFS_SING_FORCE=args.nfs_sing_force
    return

if __name__ == "__main__":
    parse_args()
    print_banner()
    #cProfile.runctx("QSv3_simd.main(keysize,workers,debug,base,key,lin_size,quad_size)", globals(), locals(), "Profile.prof")

    #s = pstats.Stats("Profile.prof")
    #s.strip_dirs().sort_stats("time").print_stats()
    QSv3_simd.main(keysize,workers,debug,base,key,lin_size,quad_size,mode)
