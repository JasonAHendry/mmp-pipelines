# MMPP: Mobile Malaria Project Pipelines
# ----------------------------------------
# Functions used with `run_mutation_search.py`
# ----------------------------------------
# JHendry, 2019/03/27


import re

def process_pileup(pileup, ref):
    """
    Process a read pileup from 
    `samtools mpileup` into
    A,T,C,G,+,-
    
    Note, the solution is ugly, but
    re.sub(\+([0-9])[ATCGatcg]+, ... fails when a variant
    follows directly after an insertion deltion.
    
    """
    
    n_pileup = len(pileup)
    processed_pileup = ''
    i = 0
    while i < n_pileup:
        p = pileup[i]
        if p in ["+", "-"]:
            size = int(pileup[i + 1])
            j = size + 2
            processed_pileup += p
        elif p == "*":
            j = 1
            processed_pileup += "-"
        else:
            j = 1
            processed_pileup += p
        i += j
    
    processed_pileup = re.sub("\$|\^\]", "", processed_pileup)
    processed_pileup = re.sub("\.|,", ref, processed_pileup)
    processed_pileup = processed_pileup.upper()
    
    return processed_pileup


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

def codon_to_amino(codon, genetic_code):
    if codon in genetic_code.keys():
        amino = genetic_code[codon]
    else:
        amino = None
    return amino