"""
for use with:
- alu_consensus_gc.ipynb (partitioining alu consensus sequence regions)
- MSE_extremes_vs_Alu_consensus.ipynb (comparing top/bottom scoring Alus 
    against their deviation from their consensus sequence 

"""

import os
import gzip

from Bio import SeqIO, pairwise2
from Bio.Seq import Seq
import numpy as np 
import pandas as pd
import pysam

# Paths - edit for your environment
SUPREMO_DIR = os.environ.get("SUPREMO_DIR", "/pollard/home/szhang20/akita_variant_scoring")
AKITA_DIR = os.environ.get("AKITA_DIR", "/pollard/home/szhang20/akita")

#from Bio.Blast import NCBIWWW, NCBIXML

#NCBIWWW is a module to submit BLAST queries to NCBI servers on internet, retrieves as XML
#NCBIXML parses XML-formatted BLAST results 
#wait jk we don't need to blast we just want to align 

#####################################
# reading utils

fasta_path=f'{SUPREMO_DIR}/data/hg38.fa'
hg38_fasta = pysam.Fastafile(fasta_path)

def get_seq_hg38(CHR, start, end):
    '''
    extract sequence from hg38 reference given genomic coordinates
    '''
    
    seq=hg38_fasta.fetch(CHR, start, end).upper()
    return(seq)

def read_fasta(fa, is_gzip=True):
    '''
    read in fasta file and return list with header name of sequence, and sequence list  
    '''   
    seq_list=[]

    if is_gzip: 
        with gzip.open(fa, 'rt') as f:
            for record in SeqIO.parse(f, 'fasta'):
                
                #split the ID and name 
                header = record.description
                header_split = header.split(' ', 1)
              
                seq_list.append([header_split[0], header_split[1], str(record.seq)])

        f.close()
        
    else: 
        #fill this in later 
        pass
        with gzip.open(fa, "r") as fa_file:
            pass
        fa_file.close()

    return(seq_list)   

def reverse_complement(sequence):
    """
    get reverse complement of sequence
    """
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N':'N'}
    reverse_comp_sequence = ''.join(complement_dict[base] for base in reversed(sequence))
    return reverse_comp_sequence

#####################################
# For finding regions of consensus 

def find_all_subseq(seq, sub_seq):
    """
    get indices of sub sequences of 
    """
    subseq_idx_list = []
    start_idx = 0

    while start_idx < len(seq):
        idx = seq.find(substring, start_idx)

        if idx == -1:
            break

        subseq_idx_list.append(idx)
        start_idx = idx + 1

    return subseq_idx_list

def find_substring_with_exception(main_string, substring):
    len_substring = len(substring)

    for i in range(len(main_string) - len_substring + 1):
        match = True

        for j in range(len_substring):
            if j == 1:  # Skip the check for the second character
                continue

            if main_string[i + j] != substring[j]:
                match = False
                break

        if match:
            return i  # Return the starting index of the first match
        
    return -1  # Return -1 if no match is found

#for counting the polyAtail 
def count_polyA(seq):
    count = 0
    index = len(seq) - 1  # start from the end of the string

    while index >= 0 and seq[index] == 'A':
        count += 1
        index -= 1

    return count

######################
# GC related things 

#for counting GC
def calculate_gc(seq, start, end):
    seq = seq.upper()
    seq_substr=seq[start:end]
    gc_count = seq_substr.count('G') + seq_substr.count('C')
    return (gc_count / len(seq_substr)) 

###########################
# For BLAST/Alignment related things 

def align_alu_consensus(alu_new, alu_consensus, global_local='global'):
    """
    
    Defaults for both global and local:
    Gap open penalty: -0.5
    Gap extension penalty: -0.1
    Mismatch penalty: 0 (matches contribute +1, and mismatches contribute 0)

    returns:
        alignment score: alignment score between the two sequences, scaled by the length of the alu_new sequence 
        identity_percentage: percent of alignment that has perfect match 
    """

    if global_local=='global':
        alignments = pairwise2.align.globalxx(alu_new, alu_consensus, one_alignment_only=True)
        #print(alignments)
    
        #scale alignment score by length of new alu sequence 
        alignment_score=alignments[0].score/len(alu_new)
    
        #get just percentage of alignment  
        #identity_percentage=pairwise2.format_alignment(*alignments[0]).count('|') / len(alu_new)

    if global_local=='local':
        alignments = pairwise2.align.localxx(alu_new, alu_consensus, one_alignment_only=True)
        #print(alignments)
    
        #scale alignment score by length of new alu sequence 
        alignment_score=alignments[0].score/len(alu_new)
    
        #get just percentage of alignment  
        #identity_percentage=pairwise2.format_alignment(*alignments[0]).count('|') / len(alu_new)

    #alignment score and identity percentage end up being the same 
    return alignment_score

