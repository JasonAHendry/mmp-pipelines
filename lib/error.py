# MMPP: Mobile Malaria Project Pipelines
# ----------------------------------------
# Functions used to characterize the error
# rate across different amplicons
# ----------------------------------------
# JHendry, 2019/05/17

import re
import pandas as pd

def get_indels(pileup):
    """
    Return a list of all indels found
    in `pileup`
    
    params
        pileup: str
            The pileup is a string representing all
            of the bases covering a given position of
            the genome; generated using samtools mpileup.
            Indels are encoded as r"[+-]\d+[ATCG]+". e.g.
            +4AATT for an insertion of four bases.
    returns
        indels: list int
            A list of the sizes of all indels in `pileup`.
            Insertions are positive, deletions are negative.
    
    """
    inserts = [int(i) for i in re.findall(r"[+]\d+", pileup)]
    deletes = [int(d) for d in re.findall(r"[-]\d+", pileup)]
    return inserts + deletes


def convert_to_frequencies(df, total, exclude):
    """
    Convert a dataframe `df` containing counts
    which sum to `total` into their frequencies
    i.e. divide columns by `total`
    
    params
        df: DataFrame
            Each row is nucleotide position.
            Columns of the contain counts 
            for different amino acids or
            nucleotides; they sum to total.
        total: str
            Column name which contains
            the total count for each row.
        exclude: list, str
            List of column names which should
            NOT be normalized, e.g.
            ['ref', 'position'] &c.
    returns
        df_norm: DataFrame
            Data fame with the same shape as
            `df`, but where all columns except
            for `exclude` have been divided by
            `total`.
    
    """
    
    cols = df.columns
    cols_to_normalize = [c for c in cols if c not in exclude]
    df_norm = df[cols_to_normalize].div(df[total], 0)
    df_norm = pd.concat([df[exclude], df_norm], 1)
    return df_norm


def annotate_homopolymers(seq):
    """
    Annotate every homopolymer within a sequence `seq`
    
    TODO:
    - The last homopolymer is not correctly annotated
    
    """
    seq_len = len(seq)
    
    results_dt = {
        "position": np.arange(seq_len),
        "ref": [s for s in seq],
        "homo": np.zeros(seq_len, 'bool'),
        "homo_id": np.zeros(seq_len),
        "homo_len": np.zeros(seq_len),
    }
    
    l = 1
    previous = None
    for ix, current in enumerate(seq):
        if previous == current:
            l += 1
        else:
            results_dt["homo"][(ix-l):ix] = l > 1
            results_dt["homo_id"][(ix-l):ix] = "%d.%d" % (ix, l)
            results_dt["homo_len"][(ix-l):ix] = l
            previous = current
            l = 1
    return pd.DataFrame(results_dt)


def find_homopolymers(seq, min_length=1):
    """
    Find all of the homopolymers within
    a given `seq` & report their 
    length, start, stop, base
    
    """
    
    # To store
    hp_dt = {"length": [], 
             "base": [], 
             "start": [],  # first base before homopolymer, zero-indexed
             "end": []}  # first base after homopolymer, zero-indexed

    # Run
    l = 1
    previous = None
    for ix, current in enumerate(seq):
        if previous == current:  # grow
            l += 1
        elif l >= min_length:  # record & restart
            hp_dt["length"].append(l)
            hp_dt["base"].append(previous)
            hp_dt["start"].append(ix-l)
            hp_dt["end"].append(ix)
            previous = current
            l = 1
        else:  # restart
            previous = current
            l = 1
            
    return pd.DataFrame(hp_dt)