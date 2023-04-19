import os
import subprocess
import argparse
import sys
import pandas as pd


# ..........................................................................

parser = argparse.ArgumentParser()

# positional
parser.add_argument('kraken_file',
                    metavar='<kraken file>',
                    help='kraken output file')

# required flags
required = parser.add_argument_group('required flags')

required.add_argument('-t',
                      dest='taxids',
                      metavar='<taxids file>',
                      help='taxids to extract',
                      required=True)

required.add_argument('-1',
                      dest='r1',
                      metavar='<fastq>',
                      help='reads 1',
                      required=True)

required.add_argument('-o1',
                      dest='out_r1',
                      metavar='<fastq>',
                      help='output reads 1',
                      required=True)

# optional flags
optional = parser.add_argument_group('optional flags')

optional.add_argument('-2',
                      dest='r2',
                      metavar='<fastq>',
                      help='reads 2, if providing paired-end reads',
                      required=False)

optional.add_argument('-o2',
                      dest='out_r2',
                      metavar='<fastq>',
                      help='output reads 2',
                      required='-2' in sys.argv)

optional.add_argument('-s',
                      dest='summary',
                      metavar='<summary_tsv>',
                      help='produces a summary table')

# ..........................................................................

def get_matching_results(kraken_file, taxid_list):
    '''
    extract dataframe from kraken results file that includes
    only classified reads within a taxid_list
    '''
    
    # read kraken results
    df = pd.read_csv(kraken_file,
                     sep='\t',
                     header=None,
                     names=['classified', 'id', 'taxid', 'len', 'lca'])
    
    # get only lines that match taxids
    df_results = df[df['taxid'].isin(taxid_list)]
    
    # get reads list
    reads_list = df_results['id'].to_list()
    
    
    return df_results, reads_list, df['classified'].value_counts()


def extract_reads(reads_list, fastq):
    '''
    given a list of read names, extract individual fastq records
    '''
    
    # pattern to match for grep
    if len(reads_list) != 0:
        
        read_match = '\|'.join(f'@{x} ' for x in reads_list)

        # grep command
        cmd = ['grep',
               '--no-group-separator',
               '-A 3',
               read_match,
               f'{fastq}']


        process = subprocess.run(cmd,
                                 universal_newlines=True,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE)
    
        return process.stdout
    
    else:
        
        return ''
        
    


def make_summary(df_results, class_stats, out_tsv):
    '''
    produces simple summary tsv including counts for 
    unclassified, classified and matching taxids
    '''
    
    unclass_n = class_stats['U']
    class_n = class_stats['C']
    match_n = len(df_results)

    df_summary = pd.DataFrame([['unclassified', unclass_n],
                           ['classified', class_n],
                           ['match_taxid', match_n]])
    
    df_summary.to_csv(out_tsv,
                      sep='\t',
                      header=None,
                      index=None)
    
    
    return df_summary


def main():
    
    args = parser.parse_args()
    
    print()
    print('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
    print(' get_kraken_reads.py 🦑')
    print()
    # make taxid list
    print(f'>  Reading taxids...')
    try:
        taxid_list = pd.read_csv(args.taxids,
                                 sep='\t',
                                 header=None)[0].to_list()
        
        print(f'   Done.')
        
    except:
        print(f'\033[1;31m>  Could not fetch taxid from {args.taxids}.\033[0;0m')
        print('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
        sys.exit(1)
    
    
    # get df of only matching results
    print(f'>  Reading kraken results...')
    
    try:
        df_results, reads_list, class_stats = get_matching_results(args.kraken_file, taxid_list)
        
        print(f'   Done.')
        
        if args.summary:
            
            print(f'>  Saving summary table...')
            
            make_summary(df_results, class_stats, args.summary)
        
            print(f'   Done. Wrote in: {args.summary}')
        
    except:
        print(f'\033[1;31m>  Could not read {args.kraken_file}.\033[0;0m')
        print('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
        sys.exit(1)
    
    # check if its using single-end or paired end
    # write read files
    if not args.r2:
        
        print(f'   Extracting reads...')
        
        # extract reads
        process = extract_reads(reads_list, args.r1)
        
        # write into file
        f = open(args.out_r1, 'w')
        f.write(process)
        f.close()
        
        print(f'   Done. Wrote in" {args.out_r1}')
    
    else:
        print(f'>  Extracting r1 reads...')
        
         # extract reads 1
        process1 = extract_reads(reads_list, args.r1)
        
        # write into file
        f = open(args.out_r1, 'w')
        f.write(process1)
        f.close()
        
        print(f'   Done. Wrote r1 in: {args.out_r1}')
        
        print(f'>  Extracting r2 reads...')
        
        # extract reads 2
        process2 = extract_reads(reads_list, args.r2)
        
        # write into file
        f = open(args.out_r2, 'w')
        f.write(process2)
        f.close()
        
        print(f'   Done. Wrote r2 in: {args.out_r2}')
    
    print('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
        
# ..........................................................................

if __name__ == '__main__':
    main()
