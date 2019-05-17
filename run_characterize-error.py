# MMPP: Mobile Malaria Project Pipelines
# ----------------------------------------
# Characterize the error
# relative to a reference genome sequence
# in a specified gene 
# amongst the reads of a sorted `.bam` file
# ----------------------------------------
# JHendry, 2019/05/09

import os
import sys
import subprocess
import configparser
import getopt
import numpy as np
import pandas as pd
from collections import Counter
from lib.mutation import *
from lib.error import *


# Parse user inputs
try:
    opts, args = getopt.getopt(sys.argv[1:], ":b:i:d:", ["bam=", "ini=", "downsample="])
	# -b or --bam : .sorted.bam file in which to search for mutations
	# -i or --ini : .ini file which contains gene name, location, and a listed set of mutations
	# -d or --downsample : int, if set, downsample to specified number of reads to accelerate analysis?
except getopt.GetoptError:
    print("Option Error.")
    
downsample = False
for opt, value in opts:
    if opt in ("-b", "--bam"):
        # define input .bam file and ouput directory
        input_bam = value
        output_path = input_bam.replace("results", "analysis")
        output_dir = os.path.dirname(output_path)
            
    elif opt in ("-i", "--ini"):
        gene_ini = value
        gene_ini_path = os.path.dirname(gene_ini)
        config = configparser.ConfigParser()
        config.read(gene_ini)
        
        # hold gene location information
        gene_dt = {}
        gene_dt["name"] = config.get("Parameters", "name")
        gene_dt["id"] = config.get("Parameters", "id")
        gene_dt["genome"] = config.get("Parameters", "genome")
        gene_dt["gff"] = config.get("Parameters", "gff")
        
        # hold mutation information
        mutations = config.get("Parameters", "mutations").split(", ")
        n_mutations = len(mutations)
        
    elif opt in ("-d", "--downsample"):
        # mostly for testing purposes
        downsample = True  
        n_reads = int(value)
        output_dir += "/downsampling"  # re-direct output
        output_path = os.path.join(output_dir, os.path.basename(input_bam))
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)
        
         
    else:
        print("Parameter %s not recognized." % opt)
        sys.exit(2)

        
# Create a GFF describing the exons of the target gene
exon_gff = "%s.exons.gff" % os.path.join(gene_ini_path, gene_dt["name"])
cmd = "grep -E 'exon.*%s' %s > %s" % (gene_dt["id"],
                                      gene_dt["gff"],
                                      exon_gff)
os.system(cmd)
gff_columns = ["seq", "source", "feature", "start", "end", "score", "strand", "phase", "attributes"]
exon_df = pd.read_csv(exon_gff, sep="\t", header=-1)
exon_df.columns = gff_columns
exon_df.start = exon_df.start - 1  # Off by one incongruity with mpileup, unfortunately.
exon_df.to_csv(exon_gff, sep="\t", index=False, header=False)


# Create an associated BED file
exon_bed = "%s.exons.bed" % os.path.join(gene_ini_path, gene_dt["name"])
cmd = "cut -f 1,4,5 %s > %s" % (exon_gff, exon_bed)
os.system(cmd)

gene_dt["chromosome"] = exon_df.seq[0]
gene_dt["start"] = exon_df.start.iloc[0]
gene_dt["end"] = exon_df.end.iloc[-1]
gene_dt["n_exons"] = len(exon_df)
gene_dt["strand"] = exon_df.strand.iloc[0]


# Let us begin...
print("====================================================================================================")
print("MMP Mutation Search Pipeline")
print("----------------------------------------------------------------------------------------------------")
print("Gene:", gene_dt["name"])
print("Chromosome:", gene_dt["chromosome"])
print("Start:", gene_dt["start"])
print("End:", gene_dt["end"])
print("Exons:", gene_dt["n_exons"])
print("Strand:", gene_dt["strand"])
print("")
print("Searching for %d mutations." % n_mutations)
print("")
print("Input BAM:", input_bam)
print("Output path:", output_path)
print("Reference genome:", gene_dt["genome"])
print("")
print("Downsampling?", downsample)
if downsample:
    print("..to %d reads." % n_reads)
print("====================================================================================================")


# Downsample if flagged
if downsample:
    print("Downsampling from SAM file...")

    # Counting total reads can be slow...
    # total_reads = int(subprocess.check_output("samtools view -c %s" % input_bam, shell=True))
    # print(" Total Reads:", total_reads)
    # print(" Downsampling to: %d (%.02f%%)" % (n_reads, 100*n_reads/total_reads))

    # You have to downsample from the `.sam` file
    input_sam = input_bam.replace("sorted.bam", "sam")  # start from `.sam`
    dwn_sam = output_path.replace("sorted.bam", "sam")
    dwn_bam = dwn_sam[:-3] + "bam"
    dwn_sorted_bam = dwn_bam.replace("bam", "sorted.bam")

    # Downsample by shuffling lines, first extract header
    # Get header
    os.system("grep '^@' %s > %s" % (input_sam, dwn_sam))

    # Downsample
    os.system("sed '/^@/d' %s > no_header.tmp.bam" % input_sam)
    os.system("gshuf -n %d no_header.tmp.bam >> %s" % (n_reads, dwn_sam))
    os.system("rm no_header.tmp.bam")

    # Prepare file for pileup
    print("Converting to BAM...")
    os.system("samtools view -S -b %s > %s" % (dwn_sam, dwn_bam))
    print("Sorting BAM...")
    os.system("samtools sort %s -o %s" % (dwn_sam, dwn_sorted_bam))
    print("Indexing BAM...")
    os.system("samtools index %s" % (dwn_sorted_bam))
    print("Done.")

    pileup_bam = dwn_sorted_bam
else:
    print("Proceeding with all reads.")
    pileup_bam = input_bam
    

# Generate read pileup for a specific gene using samtools
pileup_path = output_path.replace("sorted.bam", "%s.pileup" % gene_dt["name"])  # I want to add the gene name here

cmd = "samtools mpileup -f %s -l %s -Q 0 -aa -B %s > %s" % (gene_dt["genome"], 
                                                            exon_bed,  # Includes exon boundaries
                                                            pileup_bam, 
                                                            pileup_path)
print("Generating pileup...")
print(cmd)
os.system(cmd)
print("Done.")


# Make necessary prepartions if the gene is on the reverse strand
if gene_dt['strand'] == '-':
    print("Gene is on reverse strand, inverting pileup.")
    print("Note: still need to perform reverse complementation.")
    pileup_path_reverse = pileup_path.replace("pileup", "reverse.pileup")
    os.system('tail -r %s > %s' % (pileup_path, pileup_path_reverse))  # tail -r is not available on all systems
    pileup_path = pileup_path_reverse
else:
    print("Gene is on forward strand.")


# Common-sense test on the pileup
pileup_ref_seq = "".join([l.split("\t")[2] for l in open(pileup_path, "r")])
start_codon = pileup_ref_seq[:3]
stop_codon = pileup_ref_seq[-3:]
len_nt = len(pileup_ref_seq)

if gene_dt["strand"] == "+":
    assert len_nt % 3 == 0
    assert start_codon == "ATG"
    assert stop_codon in ["TAA", "TAG", "TGA"]
else:
    assert len_nt % 3 == 0
    assert start_codon == "TAC"
    assert stop_codon in ["ATT", "ATC", "ACT"]
  
    
# ANALYZE ERROR ON PER-NUCLEOTIDE BASIS
mutation_dt = {
    "position": [],
    "ref": [],
    "total": [],
    "A": [],
    "T": [],
    "C": [],
    "G": [],
    "-": [],
    "+": [],
    "SNV": [],
    "error": []
}

indel_dt = {
    "position": [],
    "ref": [],
    "length": []
}

# PILEUP ACROSS ALL NUCLEOTIDES
with open(pileup_path, "r") as fn:
    for i, line in enumerate(fn):
        
        # PREPARE PILEUP
        chrom, pos, ref, coverage, pileup, _ = line.split("\t")
        processed_pileup = process_pileup(pileup, ref)
        if gene_dt["strand"] == '-':
            ref = complement_map[ref]
            processed_pileup = "".join([complement_map[base] for base in processed_pileup])
            
        # GET MUTATION FREQUENCIES IN PILEUP
        nt_frequencies = Counter(processed_pileup)
        total = sum([count for nt, count in nt_frequencies.items()])
        error = sum([count for nt, count in nt_frequencies.items() 
                     if nt != ref])
        snv = sum([count for nt, count in nt_frequencies.items() 
                   if nt != ref and not nt in ["+", "-"]])

        # Store
        mutation_dt["position"].append(i)
        mutation_dt["ref"].append(ref)
        mutation_dt["total"].append(total)
        mutation_dt["A"].append(nt_frequencies["A"])
        mutation_dt["T"].append(nt_frequencies["T"])
        mutation_dt["C"].append(nt_frequencies["C"])
        mutation_dt["G"].append(nt_frequencies["G"])
        mutation_dt["-"].append(nt_frequencies["-"])
        mutation_dt["+"].append(nt_frequencies["+"])
        mutation_dt["SNV"].append(snv)
        mutation_dt["error"].append(error)
        
        # GET INDEL LENGTH DISTRIBUTION
        indels = get_indels(pileup)
        n_indels = len(indels)

        # Store
        indel_dt["position"].extend([i]*n_indels)
        indel_dt["ref"].extend([ref]*n_indels)
        indel_dt["length"].extend(indels)
        

# ANALYZE ERROR ON PER-AMINO ACID BASIS
amino_acids = set(genetic_code.values())
amino_dt = {
    "position": [],
    "ref": [],
    "total": []
}
amino_dt.update({a: [] for a in amino_acids})

# PILEUP ACROSS ALL CODONS
with open(pileup_path, "r") as fn:
    
    # want to loop until you have a full codon
    nts = 0
    codon_position = 0
    codon_ref_nts = []
    codon_pileup = []
    for i, line in enumerate(fn):
        nts += 1
        chrom, pos, ref, coverage, pileup, _ = line.split("\t")
        processed_pileup = process_pileup(pileup, ref)
        if gene_dt["strand"] == 'reverse':
            ref = complement_map[ref]
            processed_pileup = "".join([complement_map[base] for base in processed_pileup])
        codon_ref_nts.append(ref)
        codon_pileup.append(processed_pileup)

        if nts == 3:
            codon_position += 1
            # Search the codon for non-synonymous mutations above a specified frequency
            assert len(codon_ref_nts) == 3
            assert len(codon_pileup) == 3
            
            # Pre-process from lists into strings
            codon_ref = "".join(codon_ref_nts)
            amino_ref = codon_to_amino(codon_ref, genetic_code)
            codon_pileup = ["".join(c) for c in zip(*codon_pileup)]
            
            # Get frequencies of codons (i.e. nucleotide level)
            codon_frequencies = Counter(codon_pileup)
        
            amino_ref = codon_to_amino(codon_ref, genetic_code)
            amino_frequencies = Counter([codon_to_amino(c, genetic_code) for c in codon_pileup])
            # next line removes indels which have prevented making amino acid calls
            amino_frequencies = Counter(dict([(k, v) for k, v in amino_frequencies.items() 
                                              if k != None]))

            # Store
            amino_dt["position"].append(codon_position)
            amino_dt["ref"].append(amino_ref)
            for aa in amino_acids:
                amino_dt[aa].append(amino_frequencies[aa])
            amino_dt["total"].append(sum(amino_frequencies.values()))
            
            nts = 0
            codon_ref_nts = []
            codon_pileup = []

            
# WRITE OUTPUT            
mutation_df = pd.DataFrame(mutation_dt)
indel_df = pd.DataFrame(indel_dt)
amino_df = pd.DataFrame(amino_dt)

output_labels = {
    "nt_error": mutation_df,
    "aa_error": amino_df,
    "indel_counts": indel_df
}

if downsample:
    # e.g. DHFR.001000_reads.nt_error.02.csv
    write_path_template = pileup_path.replace("pileup", "%.06d_reads.%s.csv")
    for label, df in output_labels.items():
        path = write_path_template % (n_reads, label)
        fn = os.path.basename(path)
        n = sum([fn[:-4] in f for f in os.listdir(output_dir)])
        output_path = path.replace("csv", "%.02d.csv" % n)
        df.to_csv(output_path, index=False)
else:
    # e.g. DHFR.nt_error.csv
    write_path_template =  pileup_path.replace("pileup", "%s.csv")
    for label, df in output_labels.items():
        output_path = write_path_template % label
        df.to_csv(output_path, index=False)

        
print("----------------------------------------------------------------------------------------------------")
print("Error analysis complete.")
print("====================================================================================================")
