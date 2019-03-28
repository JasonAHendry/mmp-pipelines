# mmp-pipelines
A collection of pipelines for nanopore data analysis, built by the MMP.




## Setup
With conda installed, run:
```
conda env create
```

To download reference genomes, run:
```
./download_resources.sh
```

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
