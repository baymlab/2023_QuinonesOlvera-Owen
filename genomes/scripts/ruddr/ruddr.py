# ╭───────────────────────────────────────────────────────────────────────────╮
#   ⎈ ruddr
#
#   author: Natalia Quinones-Olvera
#   email: nquinones@g.harvard.edu
# ╰───────────────────────────────────────────────────────────────────────────╯

import argparse
import subprocess
import os
import pandas as pd
from io import StringIO
from Bio import SeqIO
import random
import string
import shutil


def main_argparser():
    '''
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('assembly',
                        metavar='<assembly>',
                        help='Assembly to orient, in fasta format.')
    parser.add_argument('start_genes',
                        metavar='<start_genes>',
                        help='Protein to use to orient assembly, in fasta format.')
    parser.add_argument("-o", "--out",
                        metavar='<out>',
                        help='output fasta',
                        type=str)
    parser.add_argument('-s', '--suffix',
                        metavar='<suffix>',
                        help='suffix to add to processed assembly files',
                        type=str,
                        default='')

    return parser.parse_args()



def find_start_gene(assembly_seq, start_genes, tmp_dir, keep_tmp=False):
    '''
    '''
        
    # build db for each sequence with makeblastdb
    path_base = os.path.basename(assembly_seq)
    path_nofasta = path_base.split('.fasta')[0]
    db_out = os.path.join(tmp_dir, path_nofasta)

    cmd_db = f'makeblastdb -dbtype nucl -in {assembly_seq} -out {db_out}'

    result_db = subprocess.run(cmd_db, 
                               shell=True, 
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE, 
                               universal_newlines=True)
    
    # handle error with makeblastdb
    if result_db.stderr:
        print(f'makeblastdb encountered an error:\n {result_db.stderr}\n')
        
        return None

    # search for hit with tblastn
    else:
        cmd_tblastn = f'tblastn -db {db_out} -query {start_genes} -outfmt 7 -evalue 0.01'

        result_tblastn = subprocess.run(cmd_tblastn, 
                                        shell=True, 
                                        stdout=subprocess.PIPE,
                                        stderr=subprocess.PIPE,
                                        universal_newlines=True)
        
        # delete temporary files
        if not keep_tmp:
            shutil.rmtree(tmp_dir)

        # handle error with tblastn
        if result_tblastn.stderr:
            print(f'tblastn encountered an error:\n {result_tblastn.stderr}\n')
            
            return None
        
        # if everything works
        else:
            # make dataframe from results
            try:
                df = pd.read_csv(StringIO(result_tblastn.stdout),
                                 sep='\t',
                                 comment='#',
                                 header=None)
                
                # blast outfmt 7
                df.columns = ['query_acc', 'subj_acc', 'percent_id', 
                              'ali_len',  'missmatches', 'gap_opens', 
                              'query_start', 'query_end', 'subj_start', 
                              'subj_end', 'evalue', 'bit_score']
                
                
                # sort table to keep on top long hits with good evalue
                # TODO: Add option to modify these params
                df = df.sort_values(by=['ali_len', 'evalue'], ascending=[False, True])
                
                # keep only one choice per contig blasted
                df_choice = df.drop_duplicates(subset=['subj_acc'], keep='first')
                
                # convert all names to string for downstream comparison
                df_choice['subj_acc'] = df_choice['subj_acc'].apply(str)
                
                return df_choice
            
            # handle table with no hits
            except pd.errors.EmptyDataError:
                print(f'Did not find any hits to: {start_genes}')
                
                return None



def rotate_contig(assembly_seq, start_genes,
                  out_path, out_suffix='.oriented', circular=False,
                  keep_tmp=False, force_write=True):
    '''
    
    '''
    
    # log
    log = []
    separator = '-' * 69
    log.append(f'-{separator}-\n')
    log.append('⎈ ruddr - log\n')
    log.append(f'-{separator}-\n')
    log.append(f'\n file: {assembly_seq}\n\n')
    log.append(f'+{separator}+\n')
    log.append('| contig'.ljust(40))
    log.append('match'.ljust(20))
    log.append('rotate'.ljust(10))
    log.append('|')
    log.append('\n')
    log.append(f'+{separator}+\n')
    
    out_dir = os.path.dirname(out_path)
    
    # make name for temporary folder for blast files
    random_string = ''.join(random.choices(string.digits, k=5))
    tmp_path = os.path.join(out_dir, f'tmp_{random_string}')
    
    # find matches
    start_gene_df = find_start_gene(assembly_seq, start_genes, tmp_dir=tmp_path, keep_tmp=keep_tmp)
    
    # read fasta file
    records = list(SeqIO.parse(assembly_seq, 'fasta'))
    
    if start_gene_df is not None:
        
        # added this as a check that it is actually finding the match
        # from the table
        match_count = 0
    
        # process contigs accordingly

        # go through contigs and check if they need to be rotated
        record_write = []
        
        for contig in records:
            
            # if the contig is in the table, means it has a match
            if str(contig.id) in start_gene_df['subj_acc'].to_list():
                
                match_count =+ 1

                # pattern to query dataframe
                pattern = (start_gene_df.subj_acc == contig.id)

                # store start of the hit for the log
                start = int(start_gene_df['subj_start'][pattern])

                # find the (end - start) for the hit for current contig
                diff = int(start_gene_df['subj_end'][pattern] - 
                           start_gene_df['subj_start'][pattern])

                # if the difference is positive, the end is 'after' the start
                # it is in the correct orientation, write as is
                if diff > 0:
                    record_write.append(contig)
                    # log
                    log.append(f'| {contig.id}'.ljust(40))
                    log.append(f'match at {start}'.ljust(20))
                    log.append(f'no'.ljust(10))
                    log.append('|')
                    log.append(f'\n')


                # if the difference is negative, the end is 'before' the start
                # it is in the incorrect orientation, needs reverse complement
                else:
                    contig.seq = contig.seq.reverse_complement()
                    record_write.append(contig)
                    # log
                    log.append(f'| {contig.id}'.ljust(40))
                    log.append(f'match at {start}'.ljust(20))
                    log.append(f'yes'.ljust(10))
                    log.append('|')
                    log.append(f'\n')

            # if it not in the table, it has no match, write as is
            else:
                # append as is
                record_write.append(contig)

                # log
                log.append(f'| {contig.id}'.ljust(40))
                log.append(f'-'.ljust(20))
                log.append(f'no'.ljust(10))
                log.append('|')
                log.append(f'\n')
        
        if match_count == 0:
            raise Exception('I found matches, but failed at processing. Check names.''')


        # save
        SeqIO.write(record_write, out_path, 'fasta')
        
        # log
        log.append(f'+{separator}+')
        log.append(f'\n\nsaved at: {out_path}\n')
        print(''.join(log))
    
    
    # if it doesn't find any matches
    else:
        # write it anyways
        if force_write:
            path_base = os.path.basename(assembly_seq)
            out_path = os.path.join(out_dir, path_base)
            SeqIO.write(records, out_path, 'fasta')
            log.append(f'\nsaved at: {out_path}\n')


# ─────────────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    
    args = main_argparser()
    
    rotate_contig(args.assembly, args.start_genes, out_path=args.out, out_suffix=args.suffix)
