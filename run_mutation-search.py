# MMPP: Mobile Malaria Project Pipelines
# ----------------------------------------
# Search for the presence of 
# a listed set of non-synonymous mutations
# in a specified gene 
# amongst the reads of a sorted `.bam` file
# ----------------------------------------
# JHendry, 2019/03/27


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
        input_path = os.path.dirname(input_bam)
        output_path = input_path.replace("data", "analysis")
        if not os.path.isdir(output_path):
            os.mkdir(output_path)
            
    elif opt in ("-i", "--ini"):
        gene_ini = value
        config = configparser.ConfigParser()
        config.read(gene_ini)
        
        # hold gene location information
        gene_dt = {}
        gene_dt["name"] = config.get("Location", "name")
        gene_dt["genome"] = config.get("Location", "genome")
        gene_dt["chromosome"] = config.get("Location", "chromosome")
        gene_dt["start"] = config.getint("Location", "start")
        gene_dt["end"] = config.getint("Location", "end")
        gene_dt["strand"] = config.get("Location", "strand")
        
        # hold mutation information
        mutations = config.get("Mutations", "listed").split(", ")
        n_mutations = len(mutations)
        
    elif opt in ("-d", "--downsample"):
        # mostly for testing purposes
        downsample = True
        
    else:
        print("Parameter %s not recognized." % opt)
        sys.exit(2)

print("================================================================================")
print("MMP Mutation Search Pipeline")
print("--------------------------------------------------------------------------------")
print("Gene:", gene_dt["name"])
print("Chromosome:", gene_dt["chromosome"])
print("Start:", gene_dt["start"])
print("End:", gene_dt["end"])
print("Strand:", gene_dt["strand"])
print("")
print("Searching for %d mutations." % n_mutations)
print("")
print("Input BAM:", input_bam)
print("Output path:", output_path)
print("Reference genome:", gene_dt["genome"])
print("")
print("Downsampling?", downsample)
print("================================================================================")

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
pileup_path = pileup_bam.replace("sorted.bam", "pileup")
position = gene_dt["chromosome"] + ":" + str(gene_dt["start"]) + "-" + str(gene_dt["end"])

cmd = "samtools mpileup -f %s -r %s -Q 0 -aa -B %s > %s" % (gene_dt["genome"], 
                                                            position, 
                                                            pileup_bam, 
                                                            pileup_path)
print("Generating pileup...")
print(cmd)
os.system(cmd)
print("Done.")


# Make necessary prepartions if the gene is on the reverse strand
if gene_dt['strand'] == 'reverse':
    print("Gene is on reverse strand, inverting pileup.")
    print("Note: still need to perform reverse complementation.")
    pileup_path_reverse = pileup_path.replace("pileup", "reverse.pileup")
    os.system('tail -r %s > %s' % (pileup_path, pileup_path_reverse))
    pileup_path = pileup_path_reverse
else:
    print("Gene is on forward strand.")

      
# Search for mutations
mutations = config.get("Mutations", "listed").split(", ")
n_mutations = len(mutations)

mutation_dt = {
    "mutation": [],
    "detected": [],
    "total_count": [],
    "major_amino": [],
    "major_count": [],
    "ref_major": [],
    "ref_count": [],
    "mutation_amino": [],
    "mutation_count": [],
    "n_aminos": [],
}

for mutation in mutations:
    
    # Parse Mutation Information
    print("====================================================================================================")
    print("Searching for...")
    print("  Mutation:", mutation)
    codon = int(mutation[1:-1])
    codon_nts = np.arange(3*(codon - 1), 3*codon)
    print("  Codon:", codon)
    print("  Corresponding bases:", codon_nts)
    amino_alt = mutation[-1]
    
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
                if gene_dt["strand"] == 'reverse':
                    complement_map = {"A": "T", "T": "A", 
                                      "G": "C", "C": "G", 
                                      "-": "-", "+": "+" }
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
        major_codon, major_codon_count = codon_frequencies.most_common(1)[0]
        ref_codon_count = codon_frequencies[codon_ref]
        total_codon_count = sum(codon_frequencies.values())
        
        print("Discovered...")
        print("  Reference codon (from 3D7):", codon_ref)
        print("  Majority codon (from pileup):", major_codon)
        print("  Number of unique codons discovered (including indels):", len(codon_frequencies))
        print("  Reference codon count:", ref_codon_count)
        print("  Majority codon count:", major_codon_count)
        print("  Total codon count:", total_codon_count)
        print("")

        # Get frequencies of amino acids
        amino_ref = codon_to_amino(codon_ref, genetic_code)
        amino_frequencies = Counter([codon_to_amino(c, genetic_code) for c in codon_pileup])
        # next line removes indels which have prevented making amino acid calls
        amino_frequencies = Counter(dict([(k, v) for k, v in amino_frequencies.items() if k != None]))
        major_amino, major_amino_count = amino_frequencies.most_common(1)[0]
        ref_major = amino_ref == major_amino	# is the majority amino acid reference?
        ref_amino_count = amino_frequencies[amino_ref]
        total_amino_count = sum(amino_frequencies.values())
        
        print("  Reference amino (from 3D7):", amino_ref)
        print("  Majority amino (from pileup):", major_amino)
        print("  Number of unique aminos discovered (including indels):", len(amino_frequencies))
        print("  Reference amino count:", ref_amino_count)
        print("  Majority amino count:", major_amino_count)
        print("  Total amino count:", total_amino_count)
            
        # Finally, check for non-synonymous change of interest
        if amino_alt in amino_frequencies.keys():
            mutation_detected = True
            mutation_count = amino_frequencies[amino_alt]
        else:
            mutation_detected = False
            mutation_count = 0
            
        print("Mutation detected?:", mutation_detected)
        print("Mutation count:", mutation_count)
        print("Percent of total: %.02f%%" % (100*float(mutation_count)/total_amino_count))
        
        
        # Store output
        mutation_dt["mutation"].append(mutation)
        mutation_dt["detected"].append(mutation_detected)
        mutation_dt["total_count"].append(total_amino_count)
        mutation_dt["major_amino"].append(major_amino)
        mutation_dt["major_count"].append(major_amino_count)
        mutation_dt["ref_major"].append(ref_major)
        mutation_dt["ref_count"].append(ref_amino_count)
        mutation_dt["mutation_amino"].append(amino_alt)
        mutation_dt["mutation_count"].append(mutation_count)
        mutation_dt["n_aminos"].append(len(amino_frequencies))
        
        print("====================================================================================================")
        
# A bit of cleaning, & some derived statistics
mutation_df = pd.DataFrame(mutation_dt)
mutation_df = mutation_df[["mutation", "detected", 
                           "total_count",
                           "major_amino", "major_count",
                            "ref_major", "ref_count",
                           "mutation_amino", "mutation_count",
                           "n_aminos"]]
mutation_df.to_csv(pileup_path.replace("pileup", "%s.csv" % gene_dt["name"]), index=False)
      
print("--------------------------------------------------------------------------------")
print("Mutation search complete.")
print("================================================================================")
