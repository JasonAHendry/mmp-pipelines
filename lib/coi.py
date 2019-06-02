# MMPP: Mobile Malaria Project Pipelines
# ----------------------------------------
# Functions used to infer COI from MSP2
# ----------------------------------------
# JHendry, 2019/05/02

from scipy.cluster.hierarchy import distance, linkage, dendrogram
from collections import Counter
import numpy as np


# Convert the DNA string into integers
# in order to store it as a numpy array
dna_to_ints = {
    "A": 1,
    "T": 2,
    "C": 3,
    "G": 4,
    "-": -1
}




def process_by_cigar(start, seq, cigar):
    """
    Given a read sequence `seq`
    and an associated CIGAR string `cigar`
    Modify the sequence such that each
    base corresponds to a mapped position
    from `start`
    
    Relevant CIGAR flags:
        M: match
        D: deletion
        I: insertion
        S: soft-clip. meaning bases *are* present
        in `seq`, but have been masked. As far as I 
        have been able to discern, this means they
        were not mappable, for e.g. they could be
        adapter sequence OR they could extend beyond
        the targeted mapping region, if 
        `samtools view -L region.bed` was used.
        H: hard-clip. meaning bases *are not* present
        in `seq`. In this case, meaning is the sequence
        wasn't mappable in current place, but mappable
        elsewhere. For e.g. if the read mapped as a 
        chimera.
        
    Note that in the case where a specific region is targeted,
    clipping will occur often at the start and end of a read.
    
    params
        start: int
            The start position of the read as defined
            by the SAM file.
        seq : str, len(read)
            Nucleotide sequence of read.
        cigar : str
            CIGAR string of read.
    
    returns
        A list of tuples containing alignment positions 
        from start and the corresponding base at that position.
    
    """
    
    ix = 0
    l = ''
    mapped_seq = ''
    for c in cigar:
        if c == "M":
            l = int(l)
            mapped_seq += seq[ix:(ix + l)]
            ix += l
            l = ''
        elif c == "D":
            l = int(l)
            mapped_seq += "-"*l
            l = ''
        elif c in ["I", "S"]:  # skip insertions or soft clippings
            l = int(l)
            ix += l
            l = ''
        elif c == "H":  # with hard clips, reset but keep index
            l = '' 
        else:
            l += c
            
    positions = np.arange(start, start + len(mapped_seq))
            
    return positions, mapped_seq




def filter_to_informative_sites(gene_array, min_frequency,
                                ignore_indels=True,
                                verbose=False):
    """
    Filter an array containing a pileup of
    reads `gene_array`
    to only those sites with the second highest
    within-sample allele frequency (WSAF)
    greater than `min_frequency`
    
    params
        gene_array : ndarray, shape(n_reads, n_sites), int
            An array that contains all the reads which completely
            span the given gene of interest, typical MSP2. Nucleotides
            will have been incoded as integers using the
            `dna_to_ints` dictionary.
        min_frequency : float
            Minimum frequency of the second most common allele at the site,
            in order for the site to be considered informative
    returns
        gene_informative : ndarray(n_reads, n_informative_sites), int
            Same as gene_array, but sites have been reduced to contain
            only those satisfying `min_frequency`.
    
    """
    
    assert gene_array.dtype == 'int'
    n_reads, n_sites = gene_array.shape
    min_counts = int(min_frequency*n_reads)
    
    if verbose:
        print("Number of reads:", n_reads)
        print("Number of sites:", n_sites)
        print("Frequency min.:", min_frequency)
        print("Count min.:", min_counts)
        print("Ignore indels?:", ignore_indels)
        
    informative_site_ix = []
    for ix in np.arange(n_sites):
        allele_counts = Counter(gene_array[:, ix])
        
        if ignore_indels:
            del allele_counts[-1]
        
        n_alleles = len(sorted(allele_counts))
        if n_alleles > 1:
            second_allele, second_count = allele_counts.most_common(2)[1]
            if second_count >= min_counts:
                informative_site_ix.append(ix)
    
    n_informative = len(informative_site_ix)
    if verbose:
        print("Number of informative sites:", n_informative)
        
    if n_informative > 0:
        gene_informative = gene_array[:, informative_site_ix]
    else:
        assert False, "No informative sites."
    
    return gene_informative




def calc_distance_matrix(gene_informative, ignore_indels=True, metric='jaccard'):
    """
    Calculate a pairwise distance matrix from a 
    pileup of reads across informative sites
    
    """
    n_reads, n_sites = gene_informative.shape
    if n_reads > 1000:
        print("Greater than 1000 reads!")
        print("... consider downsampling")
    
    if ignore_indels:
        f = np.copy(gene_informative).astype("float")
        f[f == -1] = np.nan
    else:
        f = gene_informative
    
    X = distance.pdist(gene_informative, metric=metric)
    return X 