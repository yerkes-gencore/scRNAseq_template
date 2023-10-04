# Getting cell cycle genes for Mus musculus
# Based on https://github.com/hbctraining/scRNA-seq_online/blob/master/lessons/cell_cycle_scoring.md
library(RCurl)
library(AnnotationHub)
library(dplyr)

# Inferred by orthology searches against Human list published in Tirosh, I, et al.:
cc_file <- getURL("https://raw.githubusercontent.com/hbc/tinyatlas/master/cell_cycle/Mus_musculus.csv")
cc_cycle_genes <- read.csv(text = cc_file)

## These annotations are Ensemble IDs but we need gene names.
# Connect to AnnotationHub
ah <- AnnotationHub()
# Access the Ensembl database for organism
ahDb <- query(ah, pattern = c("Mus musculus", "EnsDb"), ignore.case = TRUE)
# Acquire the latest annotation files
id <- ahDb %>% mcols() %>% rownames() %>% tail(n=1)
# Download the appropriate Ensembldb databases
edb <- ah[[id]]
# Extract gene-level information from database
annotations <- genes(edb, return.type = "data.frame")
# Select annotaitons of interest
annotations <- annotations %>% dplyr::select(gene_id, gene_name, seq_name, gene_biotype, description)

## Now we can use these annotaitons to get teh corresponding gene names for the Ensembl IDs of the cell cycle genes.
# Get gene names for Ensemble IDs for each gene
cell_cycle_markers <- cc_cycle_genes %>% dplyr::left_join(annotations, by =c("geneID" = "gene_id"))
# Acquire the S phase
s_genes <- cell_cycle_markers %>% dplyr::filter(phase == "S") %>% pull("gene_name")
# Acquire the G2M phase genes
g2m_genes <- cell_cycle_markers %>% dplyr::filter(phase == "G2/M") %>% pull("gene_name")
