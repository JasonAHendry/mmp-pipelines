#!/bin/bash

# Download data resources for ont analysis
# JHendry, 2019/03/06


# REFERENCE GENOME URLS
# Pf
pf_ref_url="https://plasmodb.org/common/downloads/release-39/Pfalciparum3D7/fasta/data/"
pf_ref_fn="PlasmoDB-39_Pfalciparum3D7_Genome.fasta"

# Hs
hs_ref_url="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/"
hs_ref_fn="GRCh38_full_analysis_set_plus_decoy_hla.fa"
hs_ref_fn_fai="GRCh38_full_analysis_set_plus_decoy_hla.fa.fai"

# Ag
ag_ref_url=""
ag_ref_fn=""


# MAKE RESOURCE DIRECTORY
resource_dir="data/resources"
if [ ! -d $resource_dir ]; then
  mkdir data/resources
fi


# DOWNLOAD REFERENCE GENOMES
echo "Downloading reference genomes into $resource_dir"
echo "--------------------------------------------------------------------------------"

# Pf.
pf_dir="$resource_dir/plasmodb"
if [ ! -d $pf_dir ]; then
  mkdir $pf_dir
fi
echo "  Downloading P.f. reference genome..."
curl "$pf_ref_url/$pf_ref_fn" > "$pf_dir/$pf_ref_fn"
echo "  Done"

# Hs.
hs_dir="$resource_dir/1000g"
if [ ! -d $hs_dir ]; then
  mkdir $hs_dir
fi
echo "  Downloading H.s. reference genome..."
curl "$hs_ref_url/$hs_ref_fn" > "$hs_dir/$hs_ref_fn"
curl "$hs_ref_url/$hs_ref_fn_fai" > "$hs_dir/$hs_ref_fn_fai"
echo "  Done."

