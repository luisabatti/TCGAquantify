#Script by Luis E. Abatti. This script will analyze TCGA datasets to export differential expression analysis (DEA) from specific genes. At the end, plot a graph showing which cancer types overexpress each gene the most compared to normal tissue.
###### Required libraries
library(tidyverse)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(DESeq2)
library(RColorBrewer)
library(rstatix)
library(EDASeq)
######

CancerProject <- c("TCGA-BLCA", "TCGA-BRCA", "TCGA-CHOL", "TCGA-COAD", "TCGA-ESCA", "TCGA-GBM", "TCGA-HNSC", "TCGA-KICH", "TCGA-KIRC", "TCGA-KIRP", "TCGA-LIHC", "TCGA-LUAD", "TCGA-LUSC", "TCGA-PAAD", "TCGA-PCPG", "TCGA-PRAD", "TCGA-SARC", "TCGA-STAD", "TCGA-THCA", "TCGA-THYM", "TCGA-UCEC")
#These are the TCGA cancer projects that have solid tumor and normal tissue samples

MainDirectory <- paste0(getwd(),"/")
#Indicate working directory

RNAseqDirectory <- paste0(MainDirectory,"RNAseq_files","/")
#Indicate directory to save RNAseq files

ifelse(!dir.exists(file.path(RNAseqDirectory)), dir.create(file.path(RNAseqDirectory)), FALSE)
#Creates folder to export RNAseq files

ifelse(!dir.exists(file.path(MainDirectory,"Tables")), dir.create(file.path(MainDirectory,"Tables")), FALSE)
#Creates folder to export Tables 

ifelse(!dir.exists(file.path(MainDirectory,"Figures")), dir.create(file.path(MainDirectory,"Figures")), FALSE)
#Creates folder to export Figures

genes_of_interest <- c("SOX2", "PUM1", "ACTB", "SOX9")
#These are the genes which will be used in the analysis

gtf_hg38_annotation <- read_csv(paste0(MainDirectory,"gencode.v36.genes.csv"))
#Import hg38 gene annotation, currently TCGA uses gencode V36, can be downloaded from here: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/gencode.v42.chr_patch_hapl_scaff.annotation.gtf.gz
  
gtf_hg38_gc <- read_csv(paste0(MainDirectory,"gencode.v36.gc_content.csv"))
#Imports table with GC content for every gene, exported from ENSEMBL biomart

gtf_hg38 <- inner_join(gtf_hg38_annotation, gtf_hg38_gc, by = c("gene_id", "gene_name"))
#Joins both tables together to create a database of genes, their ids, names, and gc content

############################### PROCESS DIFFERENTIAL EXPRESSION ANALYSIS USING DESEQ2
TCGA_DEA <- data.frame(gene_id = character(),
                       baseMean = numeric(),
                       log2FoldChange = numeric(),
                       lfcSE = numeric(),
                       pvalue = numeric(),
                       padj = numeric(),
                       gene_name = character(),
                       project_id = character())
#Creates dataframe that will add the results from DESeq2 analysis bellow

####FUNCTION TO IMPORT EACH PROJECT SEPARATELY
TCGA_deseq2 <- function(importData) {
  
  message("Querying ",importData)
  
  query <- GDCquery(project = importData,
                    data.category = "Transcriptome Profiling",
                    data.type = "Gene Expression Quantification",
                    workflow.type = "STAR - Counts",
                    sample.type = c("Primary Tumor", "Solid Tissue Normal"))
  #Queries GDC for RNA-seq data
  
  message(importData,": Downloading dataset")
  
  GDCdownload(query = query,
              directory = RNAseqDirectory)
  #Download queried files
  
  message(importData,": Preparing dataset")
  
  dataDown <- GDCprepare(query = query,
                         save = F,
                         directory =  RNAseqDirectory,
                         summarizedExperiment = T)
  #Import files from each dataset into a table
  
  message(importData,": Extracting clinical info")
  
  TCGA_info <- as.data.frame(SummarizedExperiment::colData(dataDown)) %>%
    #Extract colData from SummarizedExperiment
    remove_rownames %>%
    dplyr::select(project_id, patient, sample, barcode, shortLetterCode)
  #Creates table with patient information
  
  TCGA_info %>%
    write_csv(paste0(RNAseqDirectory,importData,"/",importData,"_TCGAtools_TCGA_info.csv"))
  #Save TCGA info from each dataset
  
  message(importData,": Extracting normal samples")
  
  NT <- assay(dataDown)[,(TCGA_info %>%
                            filter(shortLetterCode == "NT")) %>%
                          .$barcode]
  #Creates NT variable with expression from normal tissue
  
  message(importData,": Extracting tumour samples")
  
  TP <- assay(dataDown)[,(TCGA_info %>%
                            filter(shortLetterCode == "TP")) %>%
                          .$barcode]
  #Creates TP variable with expression from tumor tissue
  
  message(importData,": Creating DESeq2 matrix")
  
  deseq2_table <- cbind(NT, TP) %>% .[, order(colnames(.))]
  #Creates DESeq2 matrix and sort by column names
  
  message(importData,": Creating DESeq2 phenodata")
  
  phenodata <- TCGA_info %>% 
    dplyr::select(barcode, shortLetterCode) %>% 
    arrange(barcode) %>% 
    mutate(shortLetterCode = as.factor(shortLetterCode)) %>% 
    column_to_rownames(., var = "barcode")
  #Creates phenodata matrix with patient, barcode and condition (NT or TP), sorted by sample barcode to match deseq2_table
  
  message(importData,": Normalizing reads")
  
  EDASeq_normalize <- function(Counts, gcContentFile) {
    gtf_hg38_gc <- gcContentFile %>%
      dplyr::select(gene_id, gc_content)
    #Import GC data from Ensembl bioMart
    
    rawCounts <- Counts
    #Import counts
    
    commonGenes <-  intersect(gtf_hg38_gc$gene_id, rownames(rawCounts))
    #Intersect genes that are in both GC data and in the rawCounts
    
    gcInfo <- gtf_hg38_gc %>% filter(gene_id %in% commonGenes) %>%
      column_to_rownames("gene_id")
    #Selects genes from GC table that are common to both
    
    rawCounts <- rawCounts[commonGenes, ]
    #Selects genes from counts table that are common to both
    
    set <- EDASeq::newSeqExpressionSet(as.matrix(rawCounts), featureData = gcInfo)
    #Import raw counts to EDAseq
    
    message(importData,": Within Lane Normalization")
    
    dataOffset <- EDASeq::withinLaneNormalization(set,"gc_content",
                                                  which="full",offset=TRUE)
    #Run Within Lane Normalization
    
    message(importData,": Between Lane Normalization")
    
    dataOffset <- EDASeq::betweenLaneNormalization(dataOffset,
                                                   which="full",offset=TRUE)
    #Run Between Lane Normalization
    
    return(dataOffset)
  }
  
  dataNorm <- EDASeq_normalize(deseq2_table, gtf_hg38_gc) 
  
  if(all(rownames(phenodata) %in% colnames(dataNorm))){
    print(paste(importData,":","All","rownames","in","colnames"))
  } else {
    print(paste(importData,":","Rownames","are","not","all","in","colnames"))
    break
  }
  #Check whether rows in phenodata are also in columns in the normalized dataset
  
  if(all(rownames(phenodata) == colnames(dataNorm))){
    print(paste(importData,":","All","rownames","equal","colnames"))
  } else {
    print(paste(importData,":","Rownames","are","not","qual","colnames"))
    break
  }
  #Check whether the rownames in the phenodata correspond to the colnames in the normalized dataset
  
  message(importData,": Running DESeq2")
  
  dds <- DESeqDataSetFromMatrix(EDASeq::counts(dataNorm), phenodata,  design = ~ shortLetterCode)
  #Load normalized counts and phenodata into dds, set design to shortLetterCode (TP vs NT; tumor vs. normal)
  
  dds$shortLetterCode <- relevel(dds$shortLetterCode, ref = "NT")
  #Establishes normal samples as the default for comparison
  
  message(importData,": Estimating Offset Normalization Factors")
  
  normFactors <- exp(-1 * EDASeq::offst(dataNorm))
  #Uses EDASeq offets to calculate normalization factors
  normFactors <- normFactors / exp(rowMeans(log(normFactors)))
  normalizationFactors(dds) <- normFactors
  
  message(importData,": Removing low expressing genes with normalized sum < 10 across samples")
  
  keep <- apply(counts(dds, normalized=TRUE), 1, sum) >= 10
  #Only keep genes that have a sum of >= 10 reads across samples
  
  dds <- dds[keep,]
  #Filters genes before analysis
  
  dds <- DESeq(dds, parallel=F)
  
  res <- lfcShrink(dds, coef=2, type="apeglm", parallel=F)
  #Uses DESeq2 lfcShrink algorithm to shrink log2FC
  
  message(importData,": Generating DE table")
  
  result_deseq2 <- as.data.frame(res) %>% 
    rownames_to_column(var = "gene_id") %>%
    inner_join(., gtf_hg38, by = "gene_id") %>%
    mutate(project_id = importData)
  #Inserts gene name and project_ids columns
  
  result_deseq2 %>%
    write_csv(paste0(RNAseqDirectory,importData,"/",importData,"_TCGAtools_DESeq2_results.csv"))
  #Exports a table containing the DESeq2 results from each project ID, will contain all genes and their log2fc
  
  result_deseq2_GOI <- result_deseq2 %>%
    filter(gene_name %in% genes_of_interest)
  #Extracts genes of interest from DESeq2 results
  
  ########### Extract log2-normalized reads from DESeq2 -> use this for any analysis involving ANOVA
  message(importData,": Extracting log2-normalized reads")
  
  result_deseq2_log2 <- normTransform(dds)
  #Calculates log2(counts + 1)
  
  normalized_res <- as.data.frame(assay(result_deseq2_log2)) %>%
    rownames_to_column(var = "gene_id") %>%
    inner_join(., gtf_hg38, by = "gene_id") %>%
    mutate(project_id = importData) %>%
    dplyr::select(project_id, gene_id, gene_name, everything())
  
  message(importData,": Saving log2-normalized reads RDS file")
  
  saveRDS(normalized_res, paste0(RNAseqDirectory,importData,"/",importData,"_TCGAtools_DESeq2_log2_normalized_counts.rds"))
  
  
  ########## Extract log2 CPM -> use this only if asked to show CPM, not a great normal distribution compared to log2(Counts + 1) which is also normalized to library size
  message(importData,": Extracting log2-normalized CPM/FPM")
  
  result_deseq2_log2CPM <- log2(fpm(dds, robust = TRUE) + 1)
  #Calculates log2CPM(x + 1)
  
  normalized_log2CPM <- as.data.frame(result_deseq2_log2CPM) %>%
    rownames_to_column(var = "gene_id") %>%
    inner_join(., gtf_hg38, by = "gene_id") %>%
    mutate(project_id = importData) %>%
    dplyr::select(project_id, gene_id, gene_name, everything())
  
  message(importData,": Saving log2-normalized CPM RDS file")
  
  saveRDS(normalized_log2CPM, paste0(RNAseqDirectory,importData,"/",importData,"_TCGAtools_DESeq2_log2_normalized_CPM.rds"))
  
  
  ########## Extract vst-transformed reads from DESeq2 -> use this for plotting or grouping
  message(importData,": Calculating vst")
  
  vsd <- vst(dds, blind=FALSE)
  #Calculates variance stabilizing transformation - useful for visualization or clustering
  
  normalized_vsd <- as.data.frame(assay(vsd)) %>%
    rownames_to_column(var = "gene_id") %>%
    inner_join(., gtf_hg38, by = "gene_id") %>%
    mutate(project_id = importData) %>%
    dplyr::select(project_id, gene_id, gene_name, everything())
  
  message(importData,": Saving vst-normalized counts RDSfile")
  
  saveRDS(normalized_vsd, paste0(RNAseqDirectory,importData,"/",importData,"_TCGAtools_DESeq2_vst_normalized_counts.rds"))
  
  
  message(importData,": Extracting normalized log2Counts")
  
  for (gene in genes_of_interest)
  {
    message(importData,paste(": Extracting normalized", gene, "log2Counts"))
    
    log2Counts <- normalized_res %>%
      filter(gene_name == gene) %>%
      #Filters gene of interest
      gather(4:ncol(.), key = "barcode", value = "log2Counts") %>%
      inner_join(., TCGA_info, by = c("project_id", "barcode")) %>%
      #Joins TCGA information
      mutate(sample = str_sub(sample, start = 1L, end = 15L)) %>%
      #Remove the last character of the samples (vial) so that we average expression coming from the same tissue, even from different vials. This is so we can average expression in all vials from the same sample (example: some patients had 3 samples from the same tissue)
      dplyr::select(project_id, patient, sample, shortLetterCode, gene_id, log2Counts) %>%
      group_by(project_id, patient, sample, shortLetterCode) %>%
      summarise(log2Counts = median(log2Counts)) %>%
      #Takes median of multiple RNA-seq counts from the same patient
      ungroup()
    
    log2Counts %>%
      write_csv(paste0(RNAseqDirectory,importData,"/",importData,"_TCGAtools_DESeq2_",gene,"_log2_normalized_counts.csv"))
    
  }
  
  return(result_deseq2_GOI)
  #Returns DESeq2 DE analysis between NT and TP for genes of interest
}

for (ids in CancerProject) 
{
  TCGA_DEA <- TCGA_DEA %>% full_join(., TCGA_deseq2(ids), by = c("gene_id", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj", "gene_name", "project_id"))
  #Loops through the project ids to gather differential expression data from genes of interest into a table
}

TCGA_DEA %>% 
  mutate(project_id = str_remove(project_id, "TCGA-")) %>%
  dplyr::select(project_id, gene_name, baseMean, log2FoldChange, lfcSE, pvalue, padj) %>%
  dplyr::rename(cancer_type = project_id, base_mean = baseMean, log2FC = log2FoldChange, p = pvalue) %>%
  arrange(desc(log2FC)) %>%
  write_csv(paste0(MainDirectory,"Tables/","GDC_TCGA_GOI_DEA.csv"))
#Exports DEA table


################################# IMPORT TCGA_INFO

TCGA_processInfo <- function(processInfo) {
  message(processInfo,": Reading TCGA Info")
  suppressMessages(TCGA_info_csv <- read_csv(paste0(RNAseqDirectory,processInfo,"/",processInfo,"_TCGAtools_TCGA_info.csv")))
  #Imports each file with TCGA info from each dataset
  return(TCGA_info_csv)
}

TCGA_info <- TCGA_processInfo(CancerProject[1])

for (ids in CancerProject[-1])
{
  message(ids,": Importing patient info")
  suppressMessages(TCGA_info <- full_join(TCGA_info, TCGA_processInfo(ids)))
  #Loops through the project id entries from CancerProject to process them into a single table
  message(ids,": Finished merging patient info")
}

TCGA_info %>% write_csv(paste0(RNAseqDirectory,"GDC_TCGAtools_TCGA_info.csv"))
#Export TCGA info table containing patient information

##### PLOT DEA

TCGA_colours <- c("BLCA" = "#994C00",
                  "BRCA" = "#005F00",
                  "CHOL" = "#65bb3e",
                  "COAD" = "#c74eaf",
                  "ESCA" = "#53c478",
                  "GBM" = "#db487e",
                  "HNSC" = "#676cc5",
                  "KICH" = "#d6a347",
                  "KIRC" = "#dc7333",
                  "KIRP" = "#4ec2b9",
                  "LGG" = "#dc996a",
                  "LIHC" = "#c9413a",
                  "LUAD" = "#004C99",
                  "LUSC" = "#4C0099",
                  "PAAD" = "#3c8963",
                  "PCPC" = "#cf8bc7",
                  "PRAD" = "#8f872d",
                  "SARC" = "#9cb46e",
                  "STAD" = "#660033",
                  "TGCT" = "#d69d36",
                  "THCA" = "#9f6433",
                  "THYM" = "#676b2a",
                  "UCEC" = "#606060")
#Specify colors to each cancer type

plot_TCGA_DEA <- function(plot_gene) {
  
  TCGA_DEA_plot <- TCGA_DEA %>%
    filter(gene_name == plot_gene) %>%
    mutate(project_id = str_remove(project_id, "TCGA-")) %>%
    ggplot(aes(x=log2FoldChange, y=-log10(padj), colour = as.factor(project_id), shape = as.factor(gene_name))) +
    #geom_rect(aes(xmin=1,xmax=Inf,ymin=2,ymax=Inf), fill="#F9F7F7", size = 0) +
    geom_point(size = 4, alpha = 0.8) +
    geom_errorbarh(aes(xmax = log2FoldChange + lfcSE, xmin = log2FoldChange - lfcSE), size = 0.8, height = 0.05) +
    ggrepel::geom_text_repel(aes(label = project_id), size = 5, fontface = 'bold', point.padding = unit(5, 'lines'), max.overlaps = Inf) +
    geom_hline(yintercept=(-log10(0.01)), linetype="dashed") +
    geom_vline(xintercept=1, linetype="dashed") +
    xlab(bquote(~Log[2]~ "FC over normal tissue")) +
    ylab(bquote(~-Log[10]~ "adjusted P")) +
    scale_x_continuous(limits = c(-1.5, 6), breaks=seq(-1,6,1)) +
    scale_y_continuous(trans='log10', limits = c(NA, 100), breaks=c(0.01,0.1,1,2,10,100)) +
    theme_classic() +
    scale_colour_manual(values = TCGA_colours) +
    theme(plot.title = element_text(size=10, hjust = 0.5),
          legend.position="none",
          axis.text=element_text(size=12, face="bold"),
          axis.title=element_text(size=18, face="bold"),
          strip.text = element_text(size=14, face="bold"))
  
  return(TCGA_DEA_plot)
}

for (gene in genes_of_interest) {
  ggsave(file=paste0(MainDirectory,"Figures/GDC_TCGA_DEA_",gene,"_plot.pdf"), plot=plot_TCGA_DEA(gene), width=200, height=150, dpi = 300, units = "mm")
}
