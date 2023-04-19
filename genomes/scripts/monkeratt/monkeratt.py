import os
import subprocess
import argparse
import sys


# argument parser
parser = argparse.ArgumentParser()

# positional
parser.add_argument('assembly',
                    metavar='<fasta>',
                    help='assembly to annotate')

# required flags
required = parser.add_argument_group('required flags')

required.add_argument('-e',
                      '--embl',
                      dest='embl_dir',
                      metavar='<dir>',
                      help='dir with annotation(s) to export; embl format',
                      required=True)

required.add_argument('-o',
                      '--outdir',
                      dest='outdir',
                      metavar='<dir>',
                      help='output directory',
                      required=True)

required.add_argument('-p',
                      '--prefix',
                      dest='prefix',
                      metavar='<prefix>',
                      help='result prefix; str',
                      required=True)

required.add_argument('-t',
                      '--type',
                      dest='type',
                      metavar='<Assembly|Strain|Species|etc>',
                      help='transfer-type',
                      required=True)


def subprocess_ratt(args):
    '''
    Runs ratt as a subprocess.
    '''
    print()
    print('\033[0;35m+++++++++++++++++++++++++++++++++++++++++++++++++++++++++\033[0;0m')
    print('\033[0;35m MonkeRatt 🐒 🐀 \033[0;0m')
    print()
        
    # make outdir if it doesnt exist
    if not(os.path.exists(args.outdir)):
        os.mkdir(args.outdir)
        
    # build ratt command
    ratt_cmd = f'ratt {args.embl_dir} {args.assembly} -t {args.type} -p {args.prefix}'

    # run ratt in outdir
    result = subprocess.run([ratt_cmd],
                            shell=True, 
                            cwd=args.outdir)

    # print outcome, exit if failed, else return 0
    if result.returncode == 0:
        print(f'\033[0;32m>     ratt for {args.prefix} finished.\033[0;0m')
        print('\033[0;35m+++++++++++++++++++++++++++++++++++++++++++++++++++++++++\033[0;0m')
        return 0
    
    else:
        print(f'\033[1;31m>     ratt for {args.prefix} failed.\033[0;0m')
        print('\033[0;35m+++++++++++++++++++++++++++++++++++++++++++++++++++++++++\033[0;0m')
        sys.exit(1)



def main():
    args = parser.parse_args()
    
    # converted passed paths to absolute paths
    args.assembly = os.path.abspath(args.assembly)
    args.embl_dir = os.path.abspath(args.embl_dir)
    args.outdir = os.path.abspath(args.outdir)
    
    
    subprocess_ratt(args)

# ..........................................................................

if __name__ == '__main__':
    main()