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
    parser.add_argument('-nfs_want',type=int,help='no longer used (accepted so that old command lines still run)')
    parser.add_argument('-nfs_lp2',type=int,help='nfs mode: an NFS relation may keep a large prime on both sides at once if both are below 2^nfs_lp2 (default 0 = never both)')
    parser.add_argument('-nfs_procs',type=int,help='nfs mode: number of processes that sieve for the NFS (default 1)')
    parser.add_argument('-nfs_siqs',type=int,help='nfs mode: 1 = the SIQS process runs (default), 0 = NFS only')
    parser.add_argument('-nfs_save',type=str,help='nfs mode: directory where the NFS rows are kept after every round; a later run with the same directory goes on from there')
    parser.add_argument('-la_maxw',type=int,help='nfs mode: ...or until the lightest column left is in more rows than this (default 40)')
    parser.add_argument('-la_fill',type=int,help='nfs mode: ...or until the sparse rows hold this many entries in all; bounds the memory of the matrix job (default 20000000, about 4.5 GB at 300 bits)')
    parser.add_argument('-la_dense',type=int,help='nfs mode: the matrix job eliminates light columns until this many are left, then solves densely (default 3000)')
    parser.add_argument('-nfs_lp',type=int,help='nfs mode: large prime bound in bits; a relation may keep one prime up to 2^N outside the factor base (default 24, 0 = off)')
    parser.add_argument('-nfs_base',type=int,help='nfs mode: algebraic factor base of the NFS = the first N primes of the SIQS factor base (default: all of it)')
    parser.add_argument('-la_every',type=int,help='nfs mode: try the matrix each time this many new b-smooths have arrived (default 1000)')
    parser.add_argument('-nfs_force_t',type=int,help='nfs mode, forced primes: half-width in a of the part of each line that is looked at (default 0 = ten times the current width of the region)')
    parser.add_argument('-nfs_force_lines',type=int,help='nfs mode, forced primes: number of lines looked at for each (default 150)')
    parser.add_argument('-nfs_sing_qr',type=int,help='nfs mode: above -nfs_sing_small, primes up to this bound join the rational factor base only if n is a square modulo them, the only primes a SIQS relation can hold (default -1 = up to the end of the SIQS factor base, 0 = off)')
    parser.add_argument('-nfs_sing_small',type=int,help='nfs mode: every prime up to this bound counts as core: it joins the rational factor base and never makes a relation a singleton when targets are chosen (default 200)')
    parser.add_argument('-nfs_share_lp',type=int,help='nfs mode: 1 = the SIQS keeps partial relations (one large prime up to 2^nfs_lp) and the two sieves share those large primes as matrix columns (default 1), 0 = off')
    parser.add_argument('-nfs_sing_force',type=int,help='no longer used (accepted so that old command lines still run)')
    parser.add_argument('-nfs_pool',type=int,help='nfs mode: large primes of SIQS partial relations forced on the NFS per round, smallest first (default 0 = off)')
    parser.add_argument('-nfs_slack',type=int,help='nfs mode: the NFS sieve checks a position when its logs fall short by at most this many bits, plus the room for large primes (default 10)')
    parser.add_argument('-nfs_lp2_same',type=int,help='nfs mode, with -nfs_lp2: 1 = two large primes below that bound may also sit on the same side of a relation (default), 0 = only one per side')
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
    if args.nfs_lp2 != None:
        QSv3_simd.NFS_LP2_BITS=args.nfs_lp2
    if args.nfs_procs != None:
        QSv3_simd.NFS_PROCS=args.nfs_procs
    if args.nfs_siqs != None:
        QSv3_simd.NFS_SIQS=args.nfs_siqs
    if args.nfs_save != None:
        QSv3_simd.NFS_SAVE=args.nfs_save
    if args.la_maxw != None:
        QSv3_simd.LA_MAXW=args.la_maxw
    if args.la_fill != None:
        QSv3_simd.LA_FILL=args.la_fill
    if args.la_dense != None:
        QSv3_simd.LA_DENSE=args.la_dense
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
    if args.nfs_share_lp != None:
        QSv3_simd.NFS_SHARE_LP=args.nfs_share_lp
    if args.nfs_pool != None:
        QSv3_simd.NFS_POOL=args.nfs_pool
    if args.nfs_slack != None:
        QSv3_simd.NFS_SLACK=args.nfs_slack
    if args.nfs_lp2_same != None:
        QSv3_simd.NFS_LP2_SAME=args.nfs_lp2_same
    return

if __name__ == "__main__":
    parse_args()
    print_banner()
    #cProfile.runctx("QSv3_simd.main(keysize,workers,debug,base,key,lin_size,quad_size)", globals(), locals(), "Profile.prof")

    #s = pstats.Stats("Profile.prof")
    #s.strip_dirs().sort_stats("time").print_stats()
    QSv3_simd.main(keysize,workers,debug,base,key,lin_size,quad_size,mode)
