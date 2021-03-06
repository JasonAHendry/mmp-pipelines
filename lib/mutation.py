# MMPP: Mobile Malaria Project Pipelines
# ----------------------------------------
# Functions used with:
# - `run_mutation-search.py`
# - `run_mutation-scan.py`
# - `run_characterize-error.py`
# ----------------------------------------
# JHendry, 2019/03/27


import re

complement_map = {"A": "T", "T": "A", 
                  "G": "C", "C": "G", 
                  "-": "-", "+": "+" }

genetic_code = { 
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W', 
    }


def process_pileup(pileup, ref):
    """
    Process a read pileup from
    `samtools mpileup` into
    A,T,C,G,+,-
    
    Note, the solution is ugly, but
    re.sub([-+][0-9]+[ATCGatcg]+, ... fails when a variant
    follows directly after an indel.
    
    params
        pileup: str
            String giving all nucleotides mapped to a specific
            position, with characters defined by 
            `samtools mpileup`.
        ref: str
            String giving the reference base at this position.
    returns
        pileup:
            String with length <= pileup. Characters given in
            `samtools mpileup` are converted to the appropriate
            A, T, C, and G values. Insertions are represented
            by a single '+', deletions by a single '-'.
    """
    # indels are initiated with [+-][0-9]+[ATCGatcg]
    for indel_size in set(re.findall("\d+", pileup)):
        indel_size = int(indel_size)
        pileup = re.sub(r"[+]%d[ATCGatcg]{%d}" % (indel_size, indel_size), "+", pileup)  # insertions
        pileup = re.sub(r"[-]%d[ATCGatcg]{%d}" % (indel_size, indel_size), "-", pileup)  # deletions
   
    pileup = re.sub("\*", "-", pileup)  # deletions cover multiple positions
    pileup = re.sub("\$|\^.", "", pileup)
    pileup = re.sub("\.|,", ref, pileup)
    pileup = pileup.upper()
    return pileup


def codon_to_amino(codon, genetic_code):
    """
    Convert a given `codon` to its corresponding
    amino acid
    
    params
        codon: str
            Length three string of 'ATCG' or other.
        genetic_code: dict
            Dictionary mapping all possible codons
            to their corresponding amino acids.
    returns
        amino: str
            Amino acid corresponding to input
            `codon`.
    
    """
    if codon in genetic_code.keys():
        amino = genetic_code[codon]
    else:
        amino = None
    return amino