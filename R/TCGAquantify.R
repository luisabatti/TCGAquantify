tcga_deseq2 <- function(tcga_project, output_directory, gene_annotation) {
  
  ############
  
  ### This script will analyze TCGA RNA-seq projects to export differential expression analysis (DEA) from specific genes.
  
  ### tcga_projects: should be a vector one or more of the following: c("TCGA-BLCA", "TCGA-BRCA", "TCGA-CHOL", "TCGA-COAD", "TCGA-ESCA", "TCGA-GBM", "TCGA-HNSC", "TCGA-KICH", "TCGA-KIRC", "TCGA-KIRP", "TCGA-LIHC", "TCGA-LUAD", "TCGA-LUSC", "TCGA-PAAD", "TCGA-PCPG", "TCGA-PRAD", "TCGA-SARC", "TCGA-STAD", "TCGA-THCA", "TCGA-THYM", "TCGA-UCEC")
  
  ### genes_of_interest: should be a vector of one or more ENSEMBL gene names, for example "ESR1"
  
  ### gene_annotation: should be a table with gene_id (ENSEMBL), gene_name, and gc_content
  
  ############
  
  load_packages <- c("tidyverse", "TCGAbiolinks", "SummarizedExperiment", "DESeq2", "RColorBrewer", "EDASeq")
  
    require_packages <- function(x) {
      if (!require(x, character.only = TRUE)) {
        stop(paste("Please install the following package:", x))
      }
    }
  
    for (package in load_packages) {
      require_packages(package)
    }
  
  output_directory <- file.path(output_directory)
  
  ifelse(!dir.exists(output_directory), dir.create(output_directory), FALSE)
  #Creates folder to export RNAseq files
  
  message("Querying ",tcga_project)
  
  query <- GDCquery(project = tcga_project,
                    data.category = "Transcriptome Profiling",
                    data.type = "Gene Expression Quantification",
                    workflow.type = "STAR - Counts",
                    sample.type = c("Primary Tumor", "Solid Tissue Normal"))
  #Queries GDC for RNA-seq data
  
  message(tcga_project,": Downloading dataset")
  
  GDCdownload(query = query,
              directory = output_directory)
  #Download queried files
  
  message(tcga_project,": Preparing dataset")
  
  dataDown <- GDCprepare(query = query,
                         save = F,
                         directory =  output_directory,
                         summarizedExperiment = T)
  #Import files from each dataset into a table
  
  message(tcga_project,": Extracting clinical info")
  
  TCGA_info <- as.data.frame(SummarizedExperiment::colData(dataDown)) %>%
    #Extract colData from SummarizedExperiment
    remove_rownames %>%
    dplyr::select(project_id, patient, sample, barcode, shortLetterCode)
  #Creates table with patient information
  
  TCGA_info %>%
    write_csv(paste0(output_directory,"/",tcga_project,"/",tcga_project,"_TCGAtools_TCGA_info.csv"))
  #Save TCGA info from each dataset
  
  message(tcga_project,": Extracting normal samples")
  
  NT <- assay(dataDown)[,(TCGA_info %>%
                            filter(shortLetterCode == "NT")) %>%
                          .$barcode]
  #Creates NT variable with expression from normal tissue
  
  message(tcga_project,": Extracting tumour samples")
  
  TP <- assay(dataDown)[,(TCGA_info %>%
                            filter(shortLetterCode == "TP")) %>%
                          .$barcode]
  #Creates TP variable with expression from tumor tissue
  
  message(tcga_project,": Creating DESeq2 matrix")
  
  deseq2_table <- cbind(NT, TP) %>% .[, order(colnames(.))]
  #Creates DESeq2 matrix and sort by column names
  
  message(tcga_project,": Creating DESeq2 phenodata")
  
  phenodata <- TCGA_info %>% 
    dplyr::select(barcode, shortLetterCode) %>% 
    arrange(barcode) %>% 
    mutate(shortLetterCode = as.factor(shortLetterCode)) %>% 
    column_to_rownames(., var = "barcode")
  #Creates phenodata matrix with patient, barcode and condition (NT or TP), sorted by sample barcode to match deseq2_table
  
  message(tcga_project,": Normalizing reads")
  
  EDASeq_normalize <- function(Counts, gcContentFile) {
    gc_content <- gcContentFile %>%
      dplyr::select(gene_id, gc_content)
    #Import GC data from Ensembl bioMart
    
    rawCounts <- Counts
    #Import counts
    
    commonGenes <-  intersect(gc_content$gene_id, rownames(rawCounts))
    #Intersect genes that are in both GC data and in the rawCounts
    
    gcInfo <- gc_content %>% filter(gene_id %in% commonGenes) %>%
      column_to_rownames("gene_id")
    #Selects genes from GC table that are common to both
    
    rawCounts <- rawCounts[commonGenes, ]
    #Selects genes from counts table that are common to both
    
    set <- EDASeq::newSeqExpressionSet(as.matrix(rawCounts), featureData = gcInfo)
    #Import raw counts to EDAseq
    
    message(tcga_project,": Within Lane Normalization")
    
    dataOffset <- EDASeq::withinLaneNormalization(set,"gc_content",
                                                  which="full",offset=TRUE)
    #Run Within Lane Normalization
    
    message(tcga_project,": Between Lane Normalization")
    
    dataOffset <- EDASeq::betweenLaneNormalization(dataOffset,
                                                   which="full",offset=TRUE)
    #Run Between Lane Normalization
    
    return(dataOffset)
  }
  
  dataNorm <- EDASeq_normalize(deseq2_table, gcContentFile = gene_annotation) 
  
  if(!all(rownames(phenodata) %in% colnames(dataNorm))){
    print(paste(tcga_project,":","Rownames","are","not","all","in","colnames"))
    break
  }
  #Check whether rows in phenodata are also in columns in the normalized dataset
  
  if(!all(rownames(phenodata) == colnames(dataNorm))){
    print(paste(tcga_project,":","Rownames","are","not","qual","colnames"))
    break
  }
  #Check whether the rownames in the phenodata correspond to the colnames in the normalized dataset
  
  message(tcga_project,": Running DESeq2")
  
  dds <- DESeqDataSetFromMatrix(EDASeq::counts(dataNorm), phenodata,  design = ~ shortLetterCode)
  #Load normalized counts and phenodata into dds, set design to shortLetterCode (TP vs NT; tumor vs. normal)
  
  dds$shortLetterCode <- relevel(dds$shortLetterCode, ref = "NT")
  #Establishes normal samples as the default for comparison
  
  message(tcga_project,": Estimating Offset Normalization Factors")
  
  normFactors <- exp(-1 * EDASeq::offst(dataNorm))
  #Uses EDASeq offets to calculate normalization factors
  normFactors <- normFactors / exp(rowMeans(log(normFactors)))
  normalizationFactors(dds) <- normFactors
  
  message(tcga_project,": Removing low expressing genes with normalized sum < 10 across samples")
  
  keep <- apply(counts(dds, normalized=TRUE), 1, sum) >= 10
  #Only keep genes that have a sum of >= 10 reads across samples
  
  dds <- dds[keep,]
  #Filters genes before analysis
  
  dds <- DESeq(dds, parallel=F)
  
  return(dds)
  
}

tcga_dea <- function(tcga_projects, genes_of_interest, output_directory = "./RNAseq_files/", gene_annotation) {

  ###### This function will output a table with the DESeq2 differential expression analysis for each gene and TCGA project ID
  
  results_df <- data.frame()
  
  for (tcga_project in tcga_projects) {
    dds <- tcga_deseq2(tcga_project, output_directory, gene_annotation)
    
    res <- lfcShrink(dds, coef=2, type="apeglm", parallel=F)
    #Uses DESeq2 lfcShrink algorithm to shrink log2FC
    
    message(tcga_project,": Generating DE table")
    
    result_deseq2 <- as.data.frame(res) %>% 
      rownames_to_column(var = "gene_id") %>%
      inner_join(., gene_annotation, by = "gene_id") %>%
      mutate(project_id = tcga_project)
    #Inserts gene name and project_ids columns
    
    result_filtered <- result_deseq2 %>%
      filter(gene_name %in% genes_of_interest)
    #Filters for gene of interest list
    
    
    results_df <- bind_rows(results_df, result_filtered)
  }
  
  return(results_df)
}

tcga_norm_counts <- function(tcga_projects, genes_of_interest, output_directory = "./RNAseq_files/", gene_annotation) {
  ########### This function will export log2-normalized counts from DESeq2 for each gene and project ID
  
  results_df <- data.frame()
  
  for (tcga_project in tcga_projects) {
    
    message(tcga_project,": Extracting log2-normalized counts")
    
    dds <- tcga_deseq2(tcga_project, output_directory, gene_annotation)
    
    result_deseq2_log2 <- normTransform(dds)
    #Calculates log2(counts + 1)
    
    normalized_res <- as.data.frame(assay(result_deseq2_log2)) %>%
      rownames_to_column(var = "gene_id") %>%
      inner_join(., gene_annotation, by = "gene_id") %>%
      mutate(project_id = tcga_project) %>%
      dplyr::select(project_id, gene_id, gene_name, everything())
    
    TCGA_info <- read_csv(paste0(output_directory,"/",tcga_project,"/",tcga_project,"_TCGAtools_TCGA_info.csv"))
    
    log2Counts <- normalized_res %>%
      filter(gene_name %in% genes_of_interest) %>%
      #Filters gene of interest
      gather(4:ncol(.), key = "barcode", value = "log2Counts") %>%
      inner_join(., TCGA_info, by = c("project_id", "barcode")) %>%
      #Joins TCGA information
      mutate(sample = str_sub(sample, start = 1L, end = 15L)) %>%
      #Remove the last character of the samples (vial) so that we average expression coming from the same tissue, even from different vials. This is so we can average expression in all vials from the same sample (example: some patients had 3 samples from the same tissue)
      dplyr::select(project_id, patient, sample, shortLetterCode, gene_id, gene_name, log2Counts) %>%
      group_by(project_id, patient, sample, shortLetterCode, gene_id, gene_name) %>%
      summarise(log2Counts = median(log2Counts)) %>%
      #Takes median of multiple RNA-seq counts from the same patient
      ungroup()
    
    results_df <- bind_rows(results_df, log2Counts)
    
  }
  
  return(results_df)
}