# ╭───────────────────────────────────────────────────────────────────────────╮
#   FIGURE - Nucleotide diversity
#
#   author: Natalia Quinones-Olvera
#   author: Maximillian G. Marin
#   email: nquinones@g.harvard.edu
# ╰───────────────────────────────────────────────────────────────────────────╯


import pandas as pd
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import numpy as np


# SETUP
# -----------------------------------------------------------------------------

MAIN_DIR = 'data'
GENOMES_DIR = '../../../genomes/data/'

def main_dir(child):
    return os.path.join(MAIN_DIR, child)

def genomes_dir(child):
    return os.path.join(GENOMES_DIR, child)

# -----------------------------------------------------------------------------

DEBUG = True
DAAG = True

# sample metadata
df = pd.read_table(genomes_dir('metadata/phage_metadata.tsv'), sep='\t')

# reference metadata
df_ref = pd.read_table(genomes_dir('metadata/reference_genomes.tsv'), sep='\t')

# -----------------------------------------------------------------------------

if DEBUG:
    # run a subset of samples
    samples = ['PRD1', 'PRDcerulean', 'PRDobsidian']
    
    if not DAAG:
        # this is to print a debug message
        # when building dag
        print('\033[1;31m +--------------------------------------------- \033[0;0m')
        print(f'\033[1;31m | DEBUG MODE: only {samples} \033[0;0m')
        print('\033[1;31m +--------------------------------------------- \033[0;0m')
    
else:
    # get sample names
    samples = df[df['group'] == 'P']['name']


    
# Define window sizes for nucleotide diversity calculation
NucDiv_WindowSizes_bp = ['25', '50', '100']


# RULE ALL
# -----------------------------------------------------------------------------

rule all:
    input:
        consensus = expand(main_dir('consensus/{sample}.consensus.fasta'), sample=samples),
        renamed = expand(main_dir('vcf/{sample}.VCFRenamedPATH.txt'), sample=samples),
        PIs = expand(main_dir("NucDiv/NucDiv.{WindowSize}bp.windowed.pi"), WindowSize = NucDiv_WindowSizes_bp),
        NucDiv_TSVs = expand(main_dir("NucDiv/NucDiv.{WindowSize}bp.windowed.tsv"), WindowSize = NucDiv_WindowSizes_bp),
        NucDiv_SlideBy1_TSVs = expand(main_dir("NucDiv/NucDiv.{WindowSize}bp.slideby1.windowed.pi"), WindowSize = NucDiv_WindowSizes_bp)
        

# RULE: Align genomes to reference with minimap2
# -----------------------------------------------------------------------------

rule align_to_reference:
    input:
        asssembly = genomes_dir('assemblies_oriented/{sample}.fasta')
    output:
        sam = main_dir('aligned_assemblies/{sample}.sam'),
        bam = main_dir('aligned_assemblies/{sample}.bam'),
        bam_bai = main_dir('aligned_assemblies/{sample}.bam.bai')
    params:
        reference_genome = genomes_dir('reference/genomes/PRD1.fasta')
    conda:
        'envs/PRD_NucDiv_V1.yml'
    shell:
        'minimap2 '\
            '-ax asm20 '\
            '-B2 '\
            '-O6,26 '\
            '--end-bonus 100 '\
            '--cs '\
            '{params.reference_genome} '\
            '{input.asssembly}  '\
            '> {output.sam} && '\
        'samtools view '\
            '-bS {output.sam} | '\
        'samtools sort '\
            '- '\
            '> {output.bam} && '\
        'samtools index '\
            '{output.bam}'


# RULE: Make vcf file from sam file with paftools
# -----------------------------------------------------------------------------

rule make_paftools_vcf:
    input:
        sam = main_dir('aligned_assemblies/{sample}.sam')
    output:
        vcf = main_dir('vcf/{sample}.vcf')
    params:
        reference_genome = genomes_dir('reference/genomes/PRD1.fasta'),
        min_aln_len_cov = 1000,
        min_aln_len_var = 1000
    conda:
        'envs/PRD_NucDiv_V1.yml'
    shell:
        'paftools.js sam2paf '\
            '{input.sam} | '\
        'sort '\
            '-k6,6 '\
            '-k8,8n | '\
        'paftools.js call '\
            '-s {wildcards.sample} '\
            '-L {params.min_aln_len_var} '\
            '-l {params.min_aln_len_cov} '\
            '-f {params.reference_genome} '\
            '- '\
            '> {output.vcf}'

        
# RULE: Get a SNPs only vcf
# -----------------------------------------------------------------------------

rule snp_only_vcf:
    input:
        vcf = main_dir('vcf/{sample}.vcf')
    output:
        snp_vcf = main_dir('vcf/{sample}.snp.vcf')
    conda:
        'envs/PRD_NucDiv_V1.yml'
    shell:
        'bcftools view '\
            '--types snps '\
            '{input.vcf} '\
            '> {output.snp_vcf}'



# RULE: Propagate SNPs from to reference sequence to get "consensus"
# -----------------------------------------------------------------------------

rule propagate_snps:
    input:
        snp_vfc = main_dir('vcf/{sample}.snp.vcf')
    output:
        snp_vfc_gz = main_dir('vcf/{sample}.snp.vcf.gz'),
        snp_vfc_gz_idx = main_dir('vcf/{sample}.snp.vcf.gz.csi'),
        consensus_tmp = temp(main_dir('consensus/{sample}.consensus.tmp.fasta'))
    params:
        reference_genome = genomes_dir('reference/genomes/PRD1.fasta')
    conda:
        'envs/PRD_NucDiv_V1.yml'
    shell:
        'bcftools view '\
            '{input.snp_vfc} '\
            '-Oz '\
            '-o {output.snp_vfc_gz} &&'\
        'bcftools index '\
            '{output.snp_vfc_gz} &&'
        'bcftools consensus '\
            '-f {params.reference_genome} '\
            '-o {output.consensus_tmp} '\
            '{output.snp_vfc_gz}'


# RULE: Rename "consensus" fasta records to sample name
# -----------------------------------------------------------------------------

rule rename_fasta:
    input:
        consensus_tmp = main_dir('consensus/{sample}.consensus.tmp.fasta')
    output:
        consensus = main_dir('consensus/{sample}.consensus.fasta')
    run:
        record = list(SeqIO.parse(str(input.consensus_tmp), 'fasta'))[0]
        
        record.id = wildcards.sample
        record.name = ''
        record.description = ''
        
        with open(str(output.consensus), 'w') as f:
            SeqIO.write([record], f, 'fasta')

            

# RULE: Reformat the paftools VCF for merging into single VCF
# -----------------------------------------------------------------------------
       
rule rename_VCF:
    input:
        snp_vfc_gz = main_dir('vcf/{sample}.snp.vcf.gz'),
        snp_vfc_gz_idx = main_dir('vcf/{sample}.snp.vcf.gz.csi'),
    output:
        SampleID_TXT = main_dir('vcf/{sample}.name.txt'),
        vcf_renamed = main_dir('vcf/{sample}.snp.renamed.vcf.gz'),
        vcf_renamed_idx = main_dir('vcf/{sample}.snp.renamed.vcf.gz.csi'),
        vcf_renamed_PATH_TXT = main_dir('vcf/{sample}.VCFRenamedPATH.txt'),
    conda:
        'envs/PRD_NucDiv_V1.yml'
    shell:
        "echo {wildcards.sample} > {output.SampleID_TXT} \n"        
        "bcftools reheader -s {output.SampleID_TXT} {input.snp_vfc_gz} -o {output.vcf_renamed} \n"
        "bcftools index {output.vcf_renamed} \n"
        "echo '{output.vcf_renamed}' > {output.vcf_renamed_PATH_TXT} "
     
         
# RULE:
# -----------------------------------------------------------------------------       

rule merge_VCF_PATHs:
    input:
        expand(main_dir('vcf/{sample}.VCFRenamedPATH.txt'), sample = samples),
    output:
        main_dir("MergeVCFs/AllPATHs.VCFs.Renamed.txt"), 
    shell:
        "cat {input} > {output}"


# RULE:
# -----------------------------------------------------------------------------  

rule merge_All_VCFs:
    input:
        main_dir("MergeVCFs/AllPATHs.VCFs.Renamed.txt")
    output:
        main_dir("MergeVCFs/MM2.Paftools.MergedSNPs.vcf")
    conda:
        'envs/PRD_NucDiv_V1.yml'
    shell:
        'bcftools merge -i "-" -l {input} -o {output} -O v -0'


# RULE:
# -----------------------------------------------------------------------------  
rule calculate_NucDiversity:
    input:
        main_dir("MergeVCFs/MM2.Paftools.MergedSNPs.vcf")
    output:
        main_dir("NucDiv/NucDiv.{WindowSize}bp.windowed.pi")
    params:
        NucDiv_Output_Prefix = main_dir("NucDiv/NucDiv.{WindowSize}bp")
    shell:
        "vcftools --vcf {input} --window-pi {wildcards.WindowSize} --out {params.NucDiv_Output_Prefix}"

        
# RULE:
# -----------------------------------------------------------------------------  
rule calculate_NucDiversity_SlideBy1:
    input:
        main_dir("MergeVCFs/MM2.Paftools.MergedSNPs.vcf")
    output:
        main_dir("NucDiv/NucDiv.{WindowSize}bp.slideby1.windowed.pi")
    params:
        NucDiv_Output_Prefix = main_dir("NucDiv/NucDiv.{WindowSize}bp.slideby1")
    shell:
        "vcftools --vcf {input} --window-pi {wildcards.WindowSize} --window-pi-step 1 --out {params.NucDiv_Output_Prefix}"


# RULE:
# -----------------------------------------------------------------------------  
rule IncludeAllWindows_NucDiv:
    input:
        pi = main_dir("NucDiv/NucDiv.{WindowSize}bp.windowed.pi")
    output:
        div = main_dir("NucDiv/NucDiv.{WindowSize}bp.windowed.tsv")
    run:
        pi_PATH = str(input.pi)
        
        NucDiv_TSV_PATH = str(output)

        pi_DF = pd.read_csv(pi_PATH, sep = "\t")
        pi_DF["BIN_START"] = pi_DF["BIN_START"] - 1


        i_listOfTuples_NucDiv = []

        NucDiv_WindowSize = int(wildcards.WindowSize)

        Approx_GenomeSize = 15000

        for i_start in (np.arange(0, Approx_GenomeSize + NucDiv_WindowSize, NucDiv_WindowSize) ): 
            
            i_end = i_start + NucDiv_WindowSize
            i_PI_Series = pi_DF.query(f"(BIN_START == {i_start}) & (BIN_END == {i_end})")["PI"]
            
            if len(i_PI_Series) == 1: i_PI = i_PI_Series.values[0] 
                
            elif len(i_PI_Series) == 0: i_PI = 0
            
            i_row = ("NC_001421.2", i_start, i_end, i_PI)    
            i_listOfTuples_NucDiv.append(i_row)

        # Create dataframe that contains NucDiv for ALL windows
        NucDiv_DF = pd.DataFrame(i_listOfTuples_NucDiv)
        NucDiv_DF.columns = ["chrom", "start", "end", "NucDiv"]
        NucDiv_DF["Middle"] = (NucDiv_DF["start"] + NucDiv_DF["end"]) / 2

        NucDiv_DF.to_csv(NucDiv_TSV_PATH, sep = "\t", index = False)
