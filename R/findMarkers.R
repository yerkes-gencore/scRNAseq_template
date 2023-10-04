#!/usr/bin/env Rscript
## This script is intended to be ran as a background job using a processed
## Seurat object with clusters set

library(Seurat)
library(here)

## You may not want to run Seurat's implementation of MAST, as it deviates
## from the original publication's workflow, leaving it as default for now...
MAST_installed <- require('MAST')
if (MAST_installed){
  DGE_method <- 'MAST'
} else {
  DGE_method <- 'wilcox'
}

here::i_am('R/findMarkers.R')

# Load data ----
seurat_obj <- readRDS(here::here('saved_rds/obj_post-reclustering.Rds'))

# Finding markers ----
## Run find markers with desired thresholds
cluster_markers <- Seurat::FindAllMarkers(seurat_obj, assay = 'RNA', return.thresh = 0.01, min.pct = 0.25, densify = TRUE, test.use = DGE_method)

# Save results ----
saveRDS(cluster_markers, here::here('saved_rds/cluster_markers_test.Rds'), compress = FALSE)
