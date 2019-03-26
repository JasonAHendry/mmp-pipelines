# MMP P. falciparum Amplicon Analysis Pipeline
# Run Locally
# JHendry, 2019/02/18

SECONDS=0

# Parse CLI
while [[ $# -ge 1 ]]; do
tag=$1
case $tag in
-f)
	input_path=$2
	input_dir=${input_path%/*}
	input_file=${input_path##*/}
	input_id=${input_file%%.*}
	input_filetype=${input_path##*.}
	if [ $input_filetype != "fq" ] && [ $input_filetype != "fastq" ]; then
		echo "File type must be .fq or .fastq"
		exit
	fi
	shift
	;;
*)
	echo "Input" $tag "unrecognized. Skipping."
	shift
	;;
esac
shift
done

# Define Output Directories
results_path=${input_path/data/results}
results_dir=${input_dir/data/results}
if [ ! -d $results_dir ]; then
mkdir $results_dir
fi

figs_dir=${input_dir/data/figs}
if [ ! -d $figs_dir ]; then
mkdir $figs_dir
fi

log_path=${results_path/$input_filetype/"log"}
log_file=${log_path##*/}

# Log STDOUT
exec > >(tee $log_path) 2>&1

# Begin Pipeline
echo "MMP P. falciparum Amplicon Analysis Pipeline"
echo "================================================================================"
echo "Started at:" `date`
echo "By user:" `whoami`
echo ""
echo "Initializing"
echo "--------------------------------------------------------------------------------"

# Assign Reference Genomes
pf_ref_path="data/resources/plasmodb/39/PlasmoDB-39_Pfalciparum3D7_Genome.fasta"
pf_ref_file=${pf_ref_path##*/}
hs_ref_path="data/resources/1000g/GRCh38_full_analysis_set_plus_decoy_hla.fa"
hs_ref_file=${hs_ref_path##*/}

# Print Settings
echo "Input Data" 
echo "  Directory:" $input_dir
echo "  File:" $input_file
echo "  ID:" $input_id
echo "Reference Data"
echo "  Pf:" $pf_ref_path
echo "  Hs:" $hs_ref_path
echo "Output"
echo "  Results Directory:" $results_dir
echo "  Log File:"  $log_file
echo "  Figure Directory:" $figs_dir

# Adapter Trimming and De-multiplexing
echo ""
echo "Trimming And De-multiplexing Nanopore Adapter Sequences with porechop"
echo "--------------------------------------------------------------------------------"
echo "Input:" $input_path
echo "Output:" $results_dir
echo "Running..."
porechop -i $input_path -b $results_dir
echo "Done."

# Sample Summary
all_samples=`ls $results_dir/*.fastq`
n_samples=`ls $all_samples | wc -w`
i=1
echo "De-multiplexing discovered"$n_samples" samples:"
echo "...................."
for current_sample in $all_samples;
do
	sample_file=${current_sample##*/}
	echo $i": "$sample_file
	i=$((i+1))
done

# By Barcode Analysis
echo ""
echo "Beginning By Barcode Analysis"
echo "--------------------------------------------------------------------------------"

i=1
for current_sample in $all_samples;
do
	current_sample_file=${current_sample##*/}
	current_sample_id=${current_sample_file%%.*}
	echo "-------------------------------------------------------------------------------" 
	echo "BARCODE" $i":" $current_sample_file
	echo "-------------------------------------------------------------------------------"
	echo "  Path:" $current_sample
	echo "  File:" $current_sample_file
	echo "  ID:" $current_sample_id

	# Mapping to Hs
	echo ""
	echo "  Mapping to Human Reference Genome with Minimap2"
	echo "  --------------------------------------------------------------------------------"
	hs_mapped_suffix="hs.sam"
	hs_mapped_path=${current_sample/"fastq"/$hs_mapped_suffix}
	echo "  Input:" $current_sample
	echo "  Output:" $hs_mapped_path
	echo "  Running..."
	minimap2 -ax map-ont \
				$hs_ref_path \
				$current_sample \
				> $hs_mapped_path	
	echo "  Done."

	# Converting 
	echo ""
	echo "  Converting to BAM"
	echo "  ...................."
	hs_bam_suffix="hs.bam"
	hs_bam_path=${hs_mapped_path/$hs_mapped_suffix/$hs_bam_suffix}
	echo "  Input:"  $hs_mapped_path
	echo "  Output:" $hs_bam_path
	echo "  Running..."
	samtools view -S -b $hs_mapped_path > $hs_bam_path
	echo "  Done."

	# Sorting
	echo ""
	echo "  Sorting BAM"
	echo "  ...................."
	hs_sorted_suffix="hs.sorted.bam"
	hs_sorted_path=${hs_bam_path/$hs_bam_suffix/$hs_sorted_suffix}
	echo "  Input:" $hs_bam_path
	echo "  Output:" $hs_sorted_path
	echo "  Running..."
	samtools sort $hs_bam_path -o $hs_sorted_path
	echo "  Done."

	# Indexing
	echo ""
	echo "  Indexing BAM"
	echo "  ...................."
	echo "  Running..."
	samtools index $hs_sorted_path
	echo "  Done."

	# Mapping to Pf
	echo ""
	echo "  Mapping to P.falciparum Reference Genome with Minimap2"
	echo "  --------------------------------------------------------------------------------"
	pf_mapped_suffix="pf.sam"
	pf_mapped_path=${current_sample/"fastq"/$pf_mapped_suffix}
	echo "  Input:" $current_sample
	echo "  Output:" $pf_mapped_path
	echo "  Running..."
	minimap2 -ax map-ont \
			$pf_ref_path \
			$current_sample \
			> $pf_mapped_path
	echo "  Done."

	# Converting 
    echo ""
	echo "  Converting to BAM"
	echo "  ...................."
	pf_bam_suffix="pf.bam"
	pf_bam_path=${pf_mapped_path/$pf_mapped_suffix/$pf_bam_suffix}
	echo "  Input:"  $pf_mapped_path
	echo "  Output:" $pf_bam_path
	echo "  Running..."
	samtools view -S -b $pf_mapped_path > $pf_bam_path
	echo "  Done."
	
	# Sorting
	echo ""
	echo "  Sorting BAM"
	echo "  ...................."
	pf_sorted_suffix="pf.sorted.bam"
	pf_sorted_path=${pf_bam_path/$pf_bam_suffix/$pf_sorted_suffix}
	echo "  Input:" $pf_bam_path
	echo "  Output:" $pf_sorted_path
	echo "  Running..."
	samtools sort $pf_bam_path -o $pf_sorted_path
	echo "  Done."

	# Indexing
	echo ""
	echo "  Indexing BAM"
	echo "  ...................."
	echo "  Running..."
	samtools index $pf_sorted_path
	echo "  Done."
	((++i))  # increment i
done


echo "================================================================================"
echo "Finished at: "`date`
echo "Time elapsed: " $(($SECONDS / 60))m  $(($SECONDS % 60))s
echo "================================================================================"
