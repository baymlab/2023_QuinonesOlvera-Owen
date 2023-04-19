import argparse
from io import StringIO 
from Bio import SeqIO
from Bio.Data.CodonTable import TranslationError


def main_argparser():
    parser = argparse.ArgumentParser()

    # positional
    parser.add_argument('input_embl',
                        metavar='<embl>',
                        help='input embl file to convert')

    parser.add_argument('output_gb',
                        metavar='<gb>',
                        help='output genbank file')

    # required flags
    required = parser.add_argument_group('required arguments')
    required.add_argument('-p',
                          dest='prefix',
                          metavar='<prefix>',
                          help='locus name in genbank file',
                          required=True)


    optional = parser.add_argument_group('optional flags')
    optional.add_argument('-t',
                          dest='translate',
                          default=False,
                          action='store_true',
                          help='add translation to CDS features',
                          required=False)

    optional.add_argument('-m',
                          dest='molecule',
                          default='DNA',
                          metavar='<DNA|RNA>',
                          help='molecule type to add to genbank, default: DNA',
                          required=False)

    optional.add_argument('-d',
                          dest='division',
                          default='UNA',
                          metavar='<BCT|PHG|etc>',
                          help='GenBank division, default: UNA',
                          required=False)
    
    return parser.parse_args()


def add_translation(embl_file_in, force_translation=True):
    '''
    Grabs a genbank file and adds translation qualifier to feature.
    force_translation: will try to make it a valid CDS, if it fails, it will 
    translate anyways
    '''
    # list for translated records, to support multi-genbank
    new_records = []
    
    # iterate over records in genbank
    for record in SeqIO.parse(embl_file_in, 'embl'):
        
        # get nucleotide seq object
        nt_seq = record.seq
        
        # iterate over features and get CDSs
        for feature in record.features:
            if feature.type == 'CDS':
                # grab location object
                location = feature.location
                # get nt sequence in that location
                feature_seq = location.extract(nt_seq)
                # translate, with cds rules
                # https://biopython.org/docs/latest/api/Bio.Seq.html?highlight=translate#Bio.Seq.translate
                try:
                    aa = feature_seq.translate(cds=True)
                
                # if it fails, check if force_translation
                except TranslationError:
                    
                    if force_translation:
                        aa = feature_seq.translate()
                    else:
                        print(f'Translation error in {feature}')
                        aa = ''
                
                # add translation qualifier to feature
                feature.qualifiers['translation'] = aa
                
        # add record with translations to list
        new_records.append(record)
    
    return new_records
    

def embl_to_gb(args):
    '''
    Read an embl file and make it into a genbank file,
    can add translations to the CDS features
    '''
    
    # stupid header to make embl file valid
    header = 'ID                   XXXXXXXX standard; DNA; UNA; 1 BP\n'
    
    # open file to replace header
    with open(args.input_embl) as f:
        
        lines = f.readlines()
        
        # replace header
        lines[0] = header
        
        # read string as io
        io_str = StringIO(''.join(lines))
        
        # add translation if required
        if args.translate:
            record_to_save = add_translation(io_str)
        else:
            record_to_save = SeqIO.parse(io_str, 'embl')
        
        for record in record_to_save:
            record.name = args.prefix
            record.annotations['molecule_type'] = args.molecule
            record.annotations['data_file_division'] = args.division
            record.annotations['accessions'] = ''
    
    write_n = SeqIO.write(record_to_save, args.output_gb, 'genbank')
    return write_n, args.output_gb
    

def main():
    args = main_argparser()
    write_n, out_path = embl_to_gb(args)
    
    if write_n == 0:
        print(f'\033[1;31m  >     something went wrong. no files written :( \033[0;0m')
        sys.exit(1)
        
    else:
        print(f'\033[0;32m  >     success!\033[0;0m wrote {write_n} record to: {out_path}')
    

# ..........................................................................

if __name__ == '__main__':
    main()