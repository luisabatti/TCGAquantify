# TCGAquantify
An R script that analyzes TCGA RNA-seq datasets and exports differentially expressed genes from specific cancer types. At the end, it plots separate graphs showing which cancer types overexpress each gene of interest the most compared to normal tissue.

## Usage
1) Change CancerProject variable to include or remove cancer types to be analyzed;
2) Modify MainDirectory variable to change the folder in which all files will be saved;
3) Change genes_of_interest variable to modify which genes will be analyzed.

## Output
This script will output the following files:
1) .CSV table containing DESeq2 differential expresison results from each cancer type, will contain all genes and their log2 fold-change;
2) .RDS file with log2-normalized counts (from DESeq2) for all genes;
3) .RDS file with vst-transformed counts (from DESeq2) for all genes;
4) .CSV table with log2-normalized counts (from DESeq2) for each gene of interest;
5) .CSV table with differential expression analysis (from DESeq2) for each gene of interest (tumor vs. normal);
6) .CSV table with patient info for all cancer and normal samples;
7) .PDF figure with differential expression analysis (from DESeq2) for each gene of interest (tumor vs. normal).

## Requirements
This script requires the following packages:
- tidyverse
- TCGAbiolinks
- SummarizedExperiment
- DESeq2
- RColorBrewer
- rstatix
- EDASeq

Also requires a table with gene_name, gene_id, and gc_content that match the current annotation used by GDC. Both files are included in this repository, but if needed they can be downloaded from GENCODE (gene_name, gene_id) and biomart (gc_content).

## Credits
Written by Luis Abatti
