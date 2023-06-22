# TCGAquantify
A package with functions to analyze TCGA RNA-seq datasets and export differentially expressed genes from specific cancer types. 

## Requirements
This script requires the following packages:
- tidyverse
- TCGAbiolinks
- SummarizedExperiment
- DESeq2
- RColorBrewer
- EDASeq

## Usage
These functions will return the following data types:
1) tcga_deseq2(): a dds object from DESeq2 with differential expression analysis between primary tumor and normal tissue for each TCGA project ID;
2) tcga_dea(): a table with differential expression analysis (from DESeq2) for each gene of interest (tumor vs. normal) across all TCGA project IDs;
3) tcga_norm_counts(): a table with log2-normalized counts (from DESeq2) for each gene of interest across all TCGA project IDs;

## Settings
1) tcga_projects: variable to specify cancer types to be analyzed. Must be a list of these TCGA project IDs: 

"TCGA-BLCA", "TCGA-BRCA", "TCGA-CHOL", "TCGA-COAD", "TCGA-ESCA", "TCGA-GBM", "TCGA-HNSC", "TCGA-KICH", "TCGA-KIRC", "TCGA-KIRP", "TCGA-LIHC", "TCGA-LUAD", "TCGA-LUSC", "TCGA-PAAD", "TCGA-PCPG", "TCGA-PRAD", "TCGA-SARC", "TCGA-STAD", "TCGA-THCA", "TCGA-THYM", "TCGA-UCEC"

2) genes_of_interest: variable to modify which genes will be analyzed. Must be a list of ENSEMBL gene names.

3) output_directory: Folder where to save RNA-seq files downloaded from TCGA. The default is "./RNAseq_files/"

4) gene_annotation: must be a table with ENSEMBL gene_id, gene_name, and gc_content in percentage (0-1). An example table for the hg38 genome is included in the ./data folder. A similar table can be easily generated using the biomart tool from ENSEMBL. 

## Credits
Written by Luis Abatti
