"""
MMPP: Mobile Malaria Project Pipelines
--------------------
Run P.falciparum
Adapter trimming,
Demultiplexing &
Mapping
in Real-time
--------------------
JHendry, 2019/03/28

This script checks a given source directory
for new .fastq files and copies them to a 
target directory, before performing adapter 
trimming and demultiplexing (both by porechop)
and then read mapping (minimap2).

"""

import getopt
import sys
import os
import numpy as np
import time


# Parse user inputs
try:
    opts, args = getopt.getopt(sys.argv[1:], ":s:t:w:", ["source", "target", "wait"])
    # -s or --source : directory in which MinKNOW deposits .fastq files
    # -t or --target : directory in which to copy .fastq files
    # -w or --wait : if no new files are discovered, how long before next check (in seconds)?
except getopt.GetoptError:
    print("Option Error.")

for opt, value in opts:
    if opt in ("-s", "--source"):
        source_dir = value
    elif opt in ("-t", "--target"):
        target_dir = value
    elif opt in ("-w", "--wait"):
        mean_wait = float(value)
    else:
        print("Parameter %s not recognized." % opt)
        sys.exit(2)
        
        
# P.f. Reference Path
pf_ref_path = "data/resources/plasmodb/39/PlasmoDB-39_Pfalciparum3D7_Genome.fasta"

# Prepare folder structures
fastq_dir = os.path.join(target_dir, "fastq")
trimmed_dir = os.path.join(target_dir, "fastq_trimmed")

if not os.path.isdir(fastq_dir):
    os.mkdir(fastq_dir)
if not os.path.isdir(trimmed_dir):
    os.mkdir(trimmed_dir)


# Run Program
seen_files = set() # To start, we have seen no files.
while True:
    print("Checking for new files...")
    current_files = set([f for f in os.listdir(source_dir) if ".fastq" in f])
    new_files = current_files.difference(seen_files)
    if len(new_files) > 0:
        seen_files.update(new_files)  # Add new files to those seen, as they are about to be processed.
        print("Discovered %d new files." % len(new_files))
        print("Beginning processing.")
        print("--------------------------------------------------------------------------------")
        for i, file in enumerate(new_files):
            ID, ext = os.path.splitext(file)
            print("File", i)
            print("  Name:", file)
            print("  ID:", ID)
            print("  Copying to /fastq directory...")
            cmd = "cp %s %s" % (os.path.join(source_dir, file), fastq_dir)
            print("  ", cmd)
            os.system(cmd)
            print("  Done.")
            print("")
            print("  Beginning adapter trimming...")
            trimmed_output_dir = os.path.join(trimmed_dir, ID)
            if not os.path.isdir(trimmed_output_dir):
                os.mkdir(trimmed_output_dir)
            print("  Output dir:", trimmed_output_dir)
            cmd = "porechop -i %s -b %s" % (os.path.join(fastq_dir, file), trimmed_output_dir)
            print("  ", cmd)
            os.system(cmd)
            print("  Done.")
            print("")
            print("  Mapping to P.f.")
            barcode_files = [f for f in os.listdir(trimmed_output_dir) if ".fastq" in f]
            n_barcodes = len(barcode_files)
            print("  Number of barcodes discovered:", n_barcodes)
            if n_barcodes > 0:
                for j, barcode in enumerate(barcode_files):
                    print("    Mapping Barcode ", j)
                    print("    Name:", barcode)
                    barcode_fastq = os.path.join(trimmed_output_dir, barcode)
                    barcode_sam = os.path.join(target_dir, barcode.replace("fastq", "sam"))
                    cmd = 	"minimap2 -ax map-ont %s %s >> %s" % (pf_ref_path, 
                                                                  barcode_fastq,
                                                                  barcode_sam)
                    os.system(cmd)
                    print("    Done.")
                    print("")
            else:
                print("No barcodes discoverd.")
    else:
        print("No new files found, waiting...")
        time.sleep(wait_time)
    
    






