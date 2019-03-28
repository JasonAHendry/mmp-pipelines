"""
MMPP: Mobile Malaria Project Pipelines
--------------------
Simulate
the generation of
.fastq files
from a MinION
--------------------
JHendry, 2019/03/28
"""

import getopt
import sys
import os
import numpy as np
import time


# Parse user inputs
try:
    opts, args = getopt.getopt(sys.argv[1:], ":s:t:w:r:", ["source", "target", "wait", "reset"])
except getopt.GetoptError:
    print("Option Error.")

for opt, value in opts:
    if opt in ("-s", "--source"):
        source_dir = value
    elif opt in ("-t", "--target"):
        target_dir = value
    elif opt in ("-w", "--wait"):
        mean_wait = float(value)
    elif opt in ("-r", "--reset"):
        reset = bool(value)
    else:
        print("Parameter %s not recognized." % opt)
        sys.exit(2)

print("User Inputs:")
print("  Source:", source_dir)
print("  Target:", target_dir)
print("  Wait Time:", mean_wait)
print("  Reset?:", reset)

# Run file moving       
n_source = len(os.listdir(source_dir))
n_target = len(os.listdir(target_dir))
while n_source > 0:
    print("--------------------------------------------------------------------------------")
    print("Number of files in source: %d" % n_source)
    print("Number of files in target: %d" % n_target)
    print("Drawing wait time...")
    wait_time = np.random.exponential(scale=mean_wait)
    print("...wait time (s): %.02f" % wait_time)
    print("Begining wait...")
    time.sleep(wait_time)
    print("Done waiting.")
    print("Moving file:")
    source_file = os.path.join(source_dir, os.listdir(source_dir)[0])
    cmd = "mv %s %s" % (source_file, target_dir)
    print("  %s" % cmd)
    os.system(cmd)
    print("Done.")
    n_target = len(os.listdir(target_dir))
    n_source = len(os.listdir(source_dir))
    print("--------------------------------------------------------------------------------")
    
    
# Optionally reset
if reset:
    cmd="mv %s %s" % (os.path.join(target_dir, "*.fastq"), source_dir)
    os.system(cmd)
n_target = len(os.listdir(target_dir))
n_source = len(os.listdir(source_dir))
print("Final State:")
print("Number of files in source: %d" % n_source)
print("Number of files in target: %d" % n_target)
