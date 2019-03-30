# mmp-pipelines
A collection of pipelines for nanopore data analysis, built by the MMP.


## How to run the pipeline
The current pipeline allows for real-time adapter trimming and read mapping. Prepare your ONT experiment and enable real-time basecalling. Once your sequencing run has begun, type:

```
source activate mmp
python run_pf-mapping.py --source <fastq_pass> --target <data/ont/uk/2019-01-01> --wait <240>
```

The first line activates your conda environment. The second line initiates real-time adapter trimming (with `porechop`) and mapping (with `minimap2`). Base-called data from MinKNOW will be deposited in a `/fastq_pass` directory; set this as your `--source` directory. Each `.fastq` file will be moved from this directory to your `--target` before being trimmed and mapped. Your target directory should conform to `data/ont/<metadata, e.g. location>/<date>`.

Once you have collected enough data, stop MinKNOW and terminate `run_pf-mapping.py` by typing `crtl+c` in your terminal window. The final output of `run_pf-mapping.py` are `.sam` files; these can be converted to `.bam`, sorted and indexed by running:

```
python run_sam-to-bam.py --target <results/ont/uk/2019-01-01>
```

All files ending with `.sam` in the directory will be processed.

Next, you can choose to either search for a specific set of mutations (`run_mutation-search.py`) or scan for any non-synonymous mutations above a specified frequency (`run_mutation-scan.py`) in a gene of interest. E.g.:

```
python run_mutation-search.py --bam results/ont/uk/2019-01-01/BC01.sorted.bam \
                              --ini data/resource/pf-regions/kelch13.ini \
```

The above command will search the reads of `BC01.sorted.bam` for any of the mutations listed in `kelch13.ini` and deposit a table of result in `analysis/ont/uk/2019-01-01/BC01.reverse.KELCH13.search.csv`. Mutation scanning can be run in an analogous fashion.


## How to set up the pipeline
Firstly, clone the repository to your local machine:
```
git clone https://github.com/JasonAHendry/mmp-pipelines
```

Secondly, install the requisite software. With conda installed, run...
```
conda update conda
conda env create
```
Check installation by typing `porechop`, `minimap2`, `samtools`, and `bedtools`.

Thirdly, download reference genomes and deposit into the `data/resources` directory. Run:
```
./download_resources.sh
```
Check the data resources directory -- the human and *p.falciparum* reference genomes should be available.



### Tools
- **Porechop**. Oxford Nanopore Adapter Trimming. https://github.com/rrwick/Porechop
- **Minimap2**. Oxford Nanopore Read Mapping. https://github.com/lh3/minimap2
- **SAMtools**. Manipulation of SAM/BAM files. http://samtools.sourceforge.net/
- **bedtools**. Genome Arithmetic. https://bedtools.readthedocs.io/en/latest/
- **Integrative Genomics Viewer**. Visualization. http://software.broadinstitute.org/software/igv/home

### Resources
- **plasmodb**. *Plasmodium* reference genomes. http://plasmodb.org/plasmo/
- **1000g**. Human reference genome. http://www.internationalgenome.org/
- **pf-crosses**. *P. falciparum* genomic region annotations. https://www.malariagen.net/data/pf-crosses-1.0
