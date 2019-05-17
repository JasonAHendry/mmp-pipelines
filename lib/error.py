# MMPP: Mobile Malaria Project Pipelines
# ----------------------------------------
# Functions used to characterize the error
# rate across different amplicons
# ----------------------------------------
# JHendry, 2019/05/17


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