# MMPP: Mobile Malaria Project Pipelines
# ----------------------------------------
# Functions used to facilitate analysis
# of read mapping.
# ----------------------------------------
# JHendry, 2019/05/02

import subprocess

def get_mapping_statistics(bam, verbose=True):
    """
    Get mapping statistics for a given `bam` file
    using samtools
    
    Note that the statistics are all calculated 
    using the 
    
    
    Some key concepts include...
    
    Linear alignment:
    An alignment of a read to a single reference sequence
    that may include insertions, deletions, skips and clipping,
    but may not include direction changes (i.e. part forward
    and part reverse strand). A linear alignment can be 
    represented in a single SAM record.
    
    Multiple mapping: 
    The alignment of a read may be ambiguous, e.g. as a 
    result of repeats. In this case the read will map to 
    multiple locations. One of these mappings will be 
    considered primary. All others have the
    *secondary* flag set.
    
    Chimeric alignment:
    An alignment of a read that cannot be represented 
    as a linear alignment. A chimeric alignment is
    represented as a set of linear alignments that do
    not have large overlaps. Typically, one alignment
    is defined as 'representative' and all the others
    have a *supplementary* flag set.
    
    Statistics are computed using the bit-wise flag (second
    field in a SAM record). Relevant bits are include...
    
    0x4  Segment is unmapped
    0x100  Secondary alignment
    0x800  Supplementary alignment
    
    The flag "-f" retrieves all SAM records where the bit(s)
    is/are set. "-F" retrieves all SAM records where they are
    not set.
    
    
    params
        bam : str
            Path to a bam file for which
            mapping statistics are desired.
        verbose : bool
            Print the mapping statistics to screen.
    returns
        stat_dt : dict
            A dictionary containing information on
            a set of mapping statistics.
    
    """
    samtools = "samtools view -c"
    
    flag_dt = {
        "Alignments": "",
        "Reads": "-F 0x900",
        "Mapped": "-F 0x904",
        "Unmapped": "-f 0x004",
        "Multi-mapped": "-f 0x100",
        "Chimera-mapped": "-f 0x800"
    }
    
    stat_dt = {}
    for stat, flag in flag_dt.items():
        result = subprocess.check_output(" ".join([samtools, flag, bam]), shell=True)
        stat_dt[stat] = int(result)
        
    if verbose:
        print("Mapping Statistics")
        print("  BAM: %s" % bam)
        n_reads = stat_dt["Reads"]
        for stat, value in stat_dt.items():
            print("  %s: %d (%.02f%%)" % (stat, value, 100*value/n_reads))
        
    return stat_dt