# TCGAquantify
Analyze TCGA datasets to export differentially expressed genes from specific cancer types. At the end, plot a graph showing which cancer types overexpress each gene of interest the most compared to normal tissue.

## Usage
1) Change CancerProject variable to include or remove cancer types to be analyzed;
2) Modify MainDirectory variable to change the folder in which all files will be saved;
3) Change genes_of_interest variable to modify which genes will be analyzed.

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

## Requirements
Written by Luis Abatti
