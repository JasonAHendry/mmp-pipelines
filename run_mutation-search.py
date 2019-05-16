# MMPP: Mobile Malaria Project Pipelines
# ----------------------------------------
# Search for the presence of 
# a listed set of non-synonymous mutations
# in a specified gene 
# amongst the reads of a sorted `.bam` file
# ----------------------------------------
# JHendry, 2019/05/09

import os
import sys
import configparser
import getopt
import numpy as np
import pandas as pd
from collections import Counter
from lib.mutation import *


# Parse user inputs
try:
    opts, args = getopt.getopt(sys.argv[1:], ":b:i:d", ["bam=", "ini=", "downsample"])
	# -b or --bam : .sorted.bam file in which to search for mutations
	# -i or --ini : .ini file which contains gene name, location, and a listed set of mutations
	# -d or --downsample : boolean, do you want to downsample reads to accelerate analysis?
except getopt.GetoptError:
    print("Option Error.")
    
downsample = False
for opt, value in opts:
    if opt in ("-b", "--bam"):
        # define input .bam file and ouput directory
        input_bam = value
        output_path = input_bam.replace("results", "analysis")
            
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
print("====================================================================================================")


# Downsample if flagged
if downsample:
    print("Downsampling from SAM file...")
    n_reads = 1000

    # You have to downsample from the `.sam` file
    input_sam = input_bam.replace("sorted.bam", "sam")
    dwn_sam = os.path.join(output_path, os.path.basename(input_sam))
    dwn_bam = dwn_sam.replace("sam", "bam")
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
    

# Generate read pileup using samtools
pileup_path = output_path.replace("sorted.bam", "pileup")

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


# Initialize a dictionary for output statistics 
mutation_dt = {
    "mutation": [],

    "total_codon_count": [],
    "ref_codon": [],
    "ref_codon_count": [],
    "major_codon": [],
    "major_codon_count": [],
    #"n_mutation_codons": [] #  it's possible MULTIPLE codons are underlying the signal
    #"major_mutation_codon": [],
    #"major_mutation_codon_count": [],
    "n_codon_types": [],
    "indel_count": [],
    "snv_count": [],

    "total_amino_count": [],
    "ref_amino": [],
    "ref_amino_count": [],
    "major_amino": [],
    "major_amino_count": [],
    "mutation_amino": [],
    "mutation_amino_count": [],
    "n_amino_types": []
}


# Serarch across every mutation
for mutation in mutations:
    
    # Parse Mutation Information
    print("====================================================================================================")
    print("Searching for...")
    print("  Mutation:", mutation)
    codon = int(mutation[1:-1])
    codon_nts = np.arange(3*(codon - 1), 3*codon)
    print("  Codon:", codon)
    print("  Corresponding bases:", codon_nts)
    amino_mutation = mutation[-1]
    
    with open(pileup_path, "r") as fn:
        
        # Extract pileups corresponding to codon position
        codon_ref_nts = []
        codon_pileup = []
        for i, line in enumerate(fn):
            if i in codon_nts:
                
                # Extract pileup information
                chrom, pos, ref, coverage, pileup, _ = line.split("\t")
                
                # Convert the pileup string to A, T, C, G, +, -
                processed_pileup = process_pileup(pileup, ref)
                
                # Reverse complement if necessary
                if gene_dt["strand"] == '-':
                    ref = complement_map[ref]
                    processed_pileup = "".join([complement_map[base] for base in processed_pileup])
                
                # Append 
                codon_ref_nts.append(ref)
                codon_pileup.append(processed_pileup)
                
        # Check complete pileup has been extracted
        assert len(codon_ref_nts) == 3
        assert len(codon_pileup) == 3
        
        # Pre-process from lists into strings
        codon_ref = "".join(codon_ref_nts)
        codon_pileup = ["".join(c) for c in zip(*codon_pileup)]
        
        # Get frequencies of codons (i.e. nucleotide level)
        codon_frequencies = Counter(codon_pileup)
        total_codon_count = sum(codon_frequencies.values())
        ref_codon_count = codon_frequencies[codon_ref]
        major_codon, major_codon_count = codon_frequencies.most_common(1)[0]
        indel_codon_count = sum([count for c, count in codon_frequencies.items() 
                                 if "-" in c or "+" in c])
        snv_codon_count = sum([count for c, count in codon_frequencies.items() 
                               if not "-" in c and not "+" in c and c != codon_ref])
        
        print("Discovered...")
        print("  Reference codon (from 3D7):", codon_ref)
        print("  Majority codon (from pileup):", major_codon) 
        n = total_codon_count
        print("  Total codon count: %d (%.00f%%)" % (total_codon_count, 100*total_codon_count/n))
        print("  Reference codon count: %d (%.00f%%)" % (ref_codon_count, 100*ref_codon_count/n))
        print("  Majority codon count: %d (%.00f%%)" % (major_codon_count, 100*major_codon_count/n))
        print("  Indel count: %d (%.00f%%)" % (indel_codon_count, 100*indel_codon_count/n))
        print("  SNV count: %d (%.00f%%)" % (snv_codon_count, 100*snv_codon_count/n))
        print("  Unique codons:", len(codon_frequencies))

        print("")

        # Get frequencies of amino acids
        amino_ref = codon_to_amino(codon_ref, genetic_code)
        amino_frequencies = Counter([codon_to_amino(c, genetic_code) for c in codon_pileup])
        # next line removes indels which have prevented making amino acid calls
        amino_frequencies = Counter(dict([(k, v) for k, v in amino_frequencies.items() if k != None]))
        total_amino_count = sum(amino_frequencies.values())
        ref_amino_count = amino_frequencies[amino_ref]
        major_amino, major_amino_count = amino_frequencies.most_common(1)[0]
        
        print("  Reference amino (from 3D7):", amino_ref)
        print("  Majority amino (from pileup):", major_amino) 
        na = total_amino_count
        print("  Total amino count: %d (%.00f%%)" % (total_amino_count, 100*total_amino_count/na))
        print("   ... this count excludes all codons containing indels.")
        print("  Reference amino count: %d (%.00f%%)" % (ref_amino_count, 100*ref_amino_count/na))
        print("  Majority amino count: %d (%.00f%%)" % (major_amino_count, 100*major_amino_count/na))
        print("  Unique aminos:", len(amino_frequencies))
            
        # Finally, check for non-synonymous change of interest
        if amino_mutation in amino_frequencies.keys():
            mutation_detected = True
            mutation_count = amino_frequencies[amino_mutation]
        else:
            mutation_detected = False
            mutation_count = 0
            
        print("Mutation detected?:", mutation_detected)
        print("Mutation count:", mutation_count)
        print("Percent of all codons: %.02f%%" % (100*mutation_count/n))
        print("Percent of all aminos (non-indel): %.02f%%" % (100*mutation_count/na))
        
        # Store output
        mutation_dt["mutation"].append(mutation)
        mutation_dt["total_codon_count"].append(total_codon_count)
        mutation_dt["ref_codon"].append(codon_ref)
        mutation_dt["ref_codon_count"].append(ref_codon_count)
        mutation_dt["major_codon"].append(major_codon)
        mutation_dt["major_codon_count"].append(major_codon_count)
        mutation_dt["n_codon_types"].append(len(codon_frequencies))
        mutation_dt["indel_count"].append(indel_codon_count)
        mutation_dt["snv_count"].append(snv_codon_count)
        
        mutation_dt["total_amino_count"].append(total_amino_count)
        mutation_dt["ref_amino"].append(amino_ref)
        mutation_dt["ref_amino_count"].append(ref_amino_count)
        mutation_dt["major_amino"].append(major_amino)
        mutation_dt["major_amino_count"].append(major_amino_count)
        mutation_dt["mutation_amino"].append(amino_mutation)
        mutation_dt["mutation_amino_count"].append(mutation_count)
        mutation_dt["n_amino_types"].append(len(amino_frequencies))
        
        print("====================================================================================================")
    
mutation_df = pd.DataFrame(mutation_dt)
mutation_df.to_csv(pileup_path.replace("pileup", "%s.search.csv" % gene_dt["name"]), index=False)
      
print("----------------------------------------------------------------------------------------------------")
print("Mutation search complete.")
print("====================================================================================================")
