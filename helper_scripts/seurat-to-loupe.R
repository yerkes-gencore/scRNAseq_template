if (!require(loupeR)) {
  devtools::install_github('10XGenomics/loupeR')
}

library(loupeR)
library(Seurat)
library(tidyr)
## Interactive prompt to agree to license
loupeR::setup()

obj <- readRDS('/Volumes/yerkes/genomelab/illumina/runs/Analysis/2023_Analyses/p23120_Amanda/Analysis/saved_rds/LNMC/t_cells/CD8/CD8s_reprocessed.Rds')
metadata <- obj@meta.data %>%
  as.data.frame() %>%
  dplyr::select(cluster_labels, capID, HTO_classification)
metadata <- setNames(lapply(names(metadata), function(i) as.factor(metadata[[i]])), colnames(metadata))
# metadata <- split(metadata, colnames(metadata))
loupe_obj <- create_loupe(
  obj@assays$RNA@counts,
  clusters = metadata,
  output_dir = '/Users/dgratz/Downloads',
  projections = list(umap = obj@reductions$umap@cell.embeddings),
  feature_ids = rownames(obj@assays$RNA@counts),
  output_name = 'loupeR_test',
  force = TRUE
)
clusters <- readRDS('/Users/dgratz/Downloads/clusters.RDS')
