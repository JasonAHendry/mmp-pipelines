# MMPP: Mobile Malaria Project Pipelines
# ----------------------------------------
# Scan for the presence of 
# any non-synonymous mutations
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
    opts, args = getopt.getopt(sys.argv[1:], ":b:i:f:d", ["bam=", "ini=", "min_freq=", "downsample"])
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
        output_path = input_bam.replace("results", "analysis")
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
        
    elif opt in ("-f", "--min_freq"):
        min_freq = float(value)
        
    elif opt in ("-d", "--downsample"):
        # mostly for testing purposes
        downsample = True
        
    else:
        print("Parameter %s not recognized." % opt)
        sys.exit(2)

print("================================================================================")
print("MMP Mutation Scan Pipeline")
print("--------------------------------------------------------------------------------")
print("Gene:", gene_dt["name"])
print("Chromosome:", gene_dt["chromosome"])
print("Start:", gene_dt["start"])
print("End:", gene_dt["end"])
print("Strand:", gene_dt["strand"])
print("")
print("Searching for any non-synonymous mutations")
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
pileup_path = output_path.replace("sorted.bam", "pileup")
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

       
# Prepare to store mutation data
mutation_dt = {
    "position": [],
    "total_counts": [],
    "total_counts_noindel": [],
    "ref_codon": [],
    "ref_amino": [],
    "ref_freq": [],
    "major_codon": [],
    "major_amino": [],
    "major_freq": [],
    "n_mutation_types": [],
    "mutation_codon": [],
    "mutation_amino": [],
    "mutation_freq": [],
}


# Search for mutations
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
            print("--------------------------------------------------------------------------------")
            print("Checking Codon:", codon_position)
            # Search the codon for non-synonymous mutations above a specified frequency
            assert len(codon_ref_nts) == 3
            assert len(codon_pileup) == 3
            
            # Pre-process from lists into strings
            codon_ref = "".join(codon_ref_nts)
            amino_ref = codon_to_amino(codon_ref, genetic_code)
            codon_pileup = ["".join(c) for c in zip(*codon_pileup)]
            
            # Get frequencies of codons (i.e. nucleotide level)
            codon_frequencies = Counter(codon_pileup)
            major_codon, major_codon_count = codon_frequencies.most_common(1)[0]
            ref_codon_count = codon_frequencies[codon_ref]
            total_codon_count = sum(codon_frequencies.values())

            print("Discovered...")
            print("  Total codon count (including indels):", total_codon_count)
            print("  Reference codon (from 3D7):", codon_ref)
            print("  Reference codon count:", ref_codon_count)
            print("  Majority codon (from pileup):", major_codon)
            print("  Majority codon count:", major_codon_count)
            print("  Number of unique codons discovered (including indels):", len(codon_frequencies))
            print("")

            # Restrict to codon observations without indels
            # - these can be used to look for non-synonymous changes
            codon_noindel_frequencies = Counter(dict([(c, v) for c, v in codon_frequencies.items() 
                                                      if not "+" in c and not "-" in c]))
            total_noindel_count = sum(codon_noindel_frequencies.values())

            # Now get non-synonymous frequencies
            codon_nonsyn_frequencies = Counter(dict([(c, v) for c, v in codon_noindel_frequencies.items()
                                                     if codon_to_amino(c, genetic_code) != amino_ref]))
            
            
            
            
            print("  Total codon count (excluding indels):", total_noindel_count)
            print("    Unique types:", len(codon_noindel_frequencies))
            print("    Non-synonymous:", len(codon_nonsyn_frequencies))
            print("")
            if len(codon_nonsyn_frequencies) > 0:
                highest_nonsyn_codon, highest_nonsyn_count = codon_nonsyn_frequencies.most_common(1)[0]
                print("    Highest frequency non-synonymous: %s = %s" 
                      % (highest_nonsyn_codon, codon_to_amino(highest_nonsyn_codon, genetic_code)))
                print("                        at frequency: %.05f" 
                      % (highest_nonsyn_count/total_noindel_count))
                
                # Finally, check if any of the non-synonymous mutations are above
                # set threshold
                above = Counter(dict([(c, v) for c, v in codon_nonsyn_frequencies.items() 
                      if v/total_noindel_count > min_freq]))
                if len(above) > 0:
                    n_above = len(above)
                    print("Detected %d non-synonymous mutations above frequency threshold of %.02f."
                          % (n_above, min_freq))
                    print("Storing results.")
                    
                    for mut_codon, mut_counts in above.items():
                        mutation_dt["position"].append(codon_position)
                        mutation_dt["total_counts"].append(total_codon_count)
                        mutation_dt["total_counts_noindel"].append(total_noindel_count)
                        
                        mutation_dt["ref_codon"].append(codon_ref)
                        mutation_dt["ref_amino"].append(amino_ref)
                        mutation_dt["ref_freq"].append(ref_codon_count/total_codon_count)
                        
                        mutation_dt["major_codon"].append(major_codon)
                        mutation_dt["major_amino"].append(codon_to_amino(major_codon, genetic_code))
                        mutation_dt["major_freq"].append(major_codon_count/total_noindel_count)
                        
                        mutation_dt["n_mutation_types"].append(len(codon_nonsyn_frequencies))
                        mutation_dt["mutation_codon"].append(mut_codon)
                        mutation_dt["mutation_amino"].append(codon_to_amino(mut_codon, genetic_code))
                        mutation_dt["mutation_freq"].append(mut_counts/total_noindel_count)

            else:
                print("    No non-synonymous mutations discovered.")
            print("")
            
            # Once mutation search is complete, reset codon and continue
            nts = 0
            codon_ref_nts = []
            codon_pileup = []

mutation_df = pd.DataFrame(mutation_dt)
columns = ["position", "total_counts", "total_counts_noindel", 
           "ref_codon", "ref_amino", "ref_freq",
           "major_codon", "major_amino", "major_freq",
           "n_mutation_types",
           "mutation_codon", "mutation_amino", "mutation_freq"]
mutation_df = mutation_df[columns]
mutation_df.to_csv(pileup_path.replace("pileup", "%s.scan.csv" % gene_dt["name"]), index=False)
print("--------------------------------------------------------------------------------")
print("Mutation scan complete.")
print("================================================================================")


    
