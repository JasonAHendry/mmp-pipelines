# MMPP: Mobile Malaria Project Pipelines
# ----------------------------------------
# Search for the presence of 
# a listed set of non-synonymous mutations
# in a specified gene 
# amongst the reads of ALL sorted `.bam` files
# within a directory
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
	# -b or --bam : directory in which to look for .sorted.bam files
	# -i or --ini : .ini file which contains gene name, location, and a listed set of mutations
	# -d or --downsample : boolean, do you want to downsample reads to accelerate analysis?
except getopt.GetoptError:
    print("Option Error.")

downsample = False
for opt, value in opts:
    if opt in ("-b", "--bam"):
        # now this should be a directory
        bam_dir = value
        assert os.path.isdir(bam_dir)
            
    elif opt in ("-i", "--ini"):
        gene_ini = value
        
    elif opt in ("-d", "--downsample"):
        # mostly for testing purposes
        downsample = True
        
    else:
        print("Parameter %s not recognized." % opt)
        sys.exit(2)

             
bam_files = [b for b in os.listdir(bam_dir) if re.search("\.sorted\.bam$", b)]
bam_paths = [os.path.join(bam_dir, b) for b in bam_files]
n_bam_paths = len(bam_paths)

print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print("Run MMP Mutation Search Pipeline for ALL BAM files in a directory")
print("--------------------------------------------------------------------------------")
print("Target directory:", bam_dir)
print("BAM files discovered:", n_bam_paths)

std_cmd = "python run_mutation-search.py -b %s -i %s"
if downsample:
    std_cmd += " -d"
    
for i, bam_path in enumerate(bam_paths):
    print("Running BAM file:", i)
    cmd = std_cmd % (bam_path, gene_ini)
    print(cmd)
    os.system(cmd)
    
print("--------------------------------------------------------------------------------")
print("Done.")
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
