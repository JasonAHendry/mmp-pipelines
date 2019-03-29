# MMPP: Mobile Malaria Project Pipelines
# --------------------
# Convert a set of SAM files
# to sorted, indexed
# BAM files
# --------------------
# JHendry, 2019/03/28


import getopt
import sys
import os
import numpy as np
import time


# Parse user inputs
try:
    opts, args = getopt.getopt(sys.argv[1:], ":t:", ["target"])
    # -t or --target : directory in which .sam files reside
except getopt.GetoptError:
    print("Option Error.")

for opt, value in opts:
    if opt in ("-t", "--target"):
        sam_dir = value
    else:
        print("Parameter %s not recognized." % opt)
        sys.exit(2)


# Run Program
sam_files = [f for f in os.listdir(sam_dir) if ".sam" in f]
n_sam_files = len(sam_files)
print("Target Directory:", sam_dir)
print("Discovered %d SAM files." % n_sam_files)
print("Converting to BAM, sorting, and indexing...")
print("================================================================================")

for i, sam_file in enumerate(sam_files):
    print("File", i)
    print("Name:", sam_file)
    bam_file = sam_file.replace("sam", "bam")
    
    # Convert to BAM
    print("Converting to BAM...")
    cmd = "samtools view -S -b %s > %s" % (os.path.join(sam_dir, sam_file),
                                           os.path.join(sam_dir, bam_file))
    print("  ", cmd)
    os.system(cmd)
    
    # Sort BAM
    print("Sorting BAM...")
    bam_sorted_file = bam_file.replace("bam", "sorted.bam")
    cmd = "samtools sort %s -o %s" % (os.path.join(sam_dir, bam_file),
                                      os.path.join(sam_dir, bam_sorted_file))
    print("  ", cmd)
    os.system(cmd)
    
    # Index BAM
    print("Indexing BAM...")
    cmd = "samtools index %s" % os.path.join(sam_dir, bam_sorted_file)
    print("  ", cmd)
    os.system(cmd)
    
    
    print("Done.")
    print("")      
        
print("================================================================================")
