## For alternatives to what's in the template scripts.
## These are unlikely to be used frequently, if at all, but
## we can preserve the code here in case we wish to reuse them

## Doublet finder ----

#' A wrapper for DoubletFinder functions
#'
#' Takes a pre-processes Seurat object for a single capture and performs
#' DoubletFinder functions. The object is returned with an additional metadata
#' column of DoubletFinder classification calls. This can optionally be used with
#' Antibody capture data for hashed studies to provide 'ground truth' calls
#' and improve the doublet detection algorithm.
#'
#' @param obj A Seurat object for a single capture
#' @param PCs Number of PCs to use
#' @param sct Is data SCTransformed?
#' @param cores Number of cores to use for parallelization
#' @param ground_truth Optional: An nCell-length character vector of ground-truth doublet classifications (e.g., "Singlet" or "Doublet") used to gauge performance of logistic regression models trained using pANN vectors during ROC analysis.
#' @param doublet_rate Estimated doublet rate, based on sequencing device specifications, library prep, etc.
#' @param assay Assay to perform DoubletFinder on, default 'RNA'
#'
#' @returns A Seurat object with additional metadata columns, including doublet predictions
#'  and additional data used to rerun DoubletFinder with annotation data if needed
#'
#' @export
#'
#' @import Seurat
#' @import DoubletFinder
#' @import dplyr
#'
#' @examples
#' \dontrun{
#' objs <- lapply(objs, runDoubletFinder,
#'   PCs = 10,
#'   sct = FALSE,
#'   cores = 1)
#' }
# runDoubletFinder <- function(obj,
#                              PCs,
#                              sct,
#                              cores,
#                              ground_truth = NULL,
#                              doublet_rate = 0.075,
#                              assay = 'RNA'
# ){
#   Seurat::DefaultAssay(obj) <- assay
#   sweep.res.list <- DoubletFinder::paramSweep(obj,
#                                                  PCs = 1:PCs,
#                                                  sct = sct,
#                                                  num.cores = cores)
#   ## DF internally subsets to 10k cells
#   ground_truth <- ground_truth[rownames(sweep.res.list[[1]])]
#   if (!is.null(ground_truth)){
#     sweep.stats <- DoubletFinder::summarizeSweep(sweep.res.list,
#                                                  GT = TRUE,
#                                                  GT.calls = ground_truth)
#   } else {
#     sweep.stats <- DoubletFinder::summarizeSweep(sweep.res.list, GT = FALSE)
#   }
#
#   bcmvn <- DoubletFinder::find.pK(sweep.stats)
#   pK <- bcmvn[bcmvn$BCmetric==max(bcmvn$BCmetric),'pK']
#   ## Assuming 7.5% doublet formation rate - tailor for your dataset
#   nExp_poi <- round(doublet_rate*nrow(obj@meta.data))
#   obj <- DoubletFinder::doubletFinder(obj,
#                                          PCs = 1:PCs,
#                                          pK = as.numeric(as.character(pK)),
#                                          nExp = nExp_poi,
#                                          sct = sct)
#   obj@meta.data <- obj@meta.data %>%
#     dplyr::rename(DF_pANN = starts_with('pANN')) %>%
#     dplyr::rename(DF_classifications = starts_with('DF.classifications'))
#   return(obj)
# }


# objs <- lapply(objs, NormFindVarFeatScaleData, norm_method = 'logNorm')
# objs <- lapply(objs, RunPCA)
# objs_filt <- lapply(objs, runDoubletFinder,
#                     PCs = 30,
#                     sct = FALSE,
#                     cores = 4) ## adjust for your machine
#
# doubletfinder_results <- data.table::rbindlist(
#   lapply(objs_filt, function(x) {
#     as.list(table(x$DF_classifications))
#   }), idcol = 'Sample')
# doubletfinder_results %>%
#   knitr::kable()

## Azimuth ----

# objs_filt_mapped <- lapply(objs_filt,
#                            RunAzimuth,
#                            reference = 'bonemarrowref', ## A seurat/Azimuth reference dataset
#                            umap.name = 'refUMAP.bm')    ## name for reference UMAP added to query data
# for (obj in names(objs_filt_mapped)){
#   DefaultAssay(objs_filt_mapped[[obj]]) <- 'RNA'
# }
#
# #### Plotting mapping results
# bm.ref <- LoadData('bonemarrowref.SeuratData', type = 'azimuth')
#
# p1 <- DimPlot(bm.ref$map,
#               reduction = "refUMAP",    ## This will be specific to the reference
#               group.by = "celltype.l2", ## This will be specific to the reference
#               label = TRUE,
#               label.size = 3,
#               raster = FALSE) +
#   NoLegend() + ggtitle("Reference annotations")
# p2 <- DimPlot(objs_filt_mapped$R1,
#               reduction = "refUMAP.bm",           ## umap.name used in RunAzimuth
#               group.by = "predicted.celltype.l2", ## This will be specific to the reference
#               label = TRUE,
#               label.size = 3,
#               repel = TRUE) +
#   NoLegend() + ggtitle("Query transferred labels")
# p1 + p2

#### Create azimuth reference ----

#' Create an Azimuth reference and save to file
#'
#' Create an Azimuth compatible reference from an existing Seurat object
#'  with a counts matrix and cell-type annotations in a metadata column. This
#'  wrapper will process the reference in a standard SCTransform workflow.
#'  The column of metadata with annotations should be specified to the `metadata`
#'  parameter. The reference will be saved to the specified output folder to
#'  be referenced in calls to `Azimuth::RunAzimuth`.
#'
#' @param ref A Seurat object with a counts matrix and cell-type annotations
#' @param metadata_column Columns of metadata to transfer onto query datasets
#' @param output_folder Where to save the Azimuth reference files for future calls
#'  of `RunAzimuth`
#' @param ndims Number of dimensions to use for PCA and UMAP
#' @param \dots Arguments passed to Azimuth::AzimuthReference()
#'
#' @returns A Seurat object with AzimuthData stored in the tools slot for use with Azimuth
#' @export
#'
#' @note Code based on scripts from azimuth-references github, such as
#'  https://github.com/satijalab/azimuth-references/blob/master/human_motorcortex/scripts/export.R
#'
#' @import Seurat
#' @import dplyr
#'
#' @examples
#' \dontrun{
#'   ## Load reference data as Seurat object
#'   ref <- readRDS(here('saved_rds/tabula_sapiens_lymph_subset.Rds'))
#'   ## No processing needed, plug in
#'   azimuth_ref <- createAzimuthReference(ref,
#'                                         'cell_type',
#'                                         output_folder = here('reference/lymph_node_TS'))
#'
#'   ## Don't need the reference object for this, just the files
#'   combined.obj <- RunAzimuth(combined.obj,
#'                              reference = here('reference/lymph_node_TS'),
#'                              umap.name = 'refUMAP.lnTS',
#'                              assay = 'SCT')
#' }
# createAzimuthReference <- function(ref,
#                                    metadata_column,
#                                    output_folder,
#                                    ndims = 50,
#                                    ...){
#   if (!dir.exists(output_folder)){
#     message('Creating output dir')
#     dir.create(output_folder, recursive = TRUE)
#   } else {
#     files <- dir(output_folder)
#     if ('idx.annoy' %in% files | 'ref.Rds' %in% files){
#       stop('Reference files already exist at the destination provided by "output_folder".\n
#             Please provide a different path')
#     }
#   }
#   if (ndims < 50){
#     stop('Azimuth requires references to be built with at least 50 dimensions')
#   }
#
#   ## in case the object is a subset, remove unwanted label factors
#   ref@meta.data <- ref@meta.data %>% droplevels()
#
#   ref <- Seurat::SCTransform(ref,
#                              #method = "glmGamPoi",
#                              new.assay.name = "SCT")
#   ref <- Seurat::RunPCA(ref,
#                         assay='SCT',
#                         npcs = ndims)
#   ## Return model for later map projections
#   ref <- Seurat::RunUMAP(ref,
#                          dims = 1:ndims,
#                          reduction = 'pca',
#                          return.model = TRUE,
#                          umap.method = "uwot")
#   ## Create azimuth compatible reference
#   ## Can't specify namespace here due to weird stuff with a tools() call
#   ## https://github.com/satijalab/azimuth/issues/155
#   ## And can't do an import here cause I want to leave azimuth
#   ## as a suggested package, so settling for 2 notes
#   ref <- AzimuthReference(object = ref,
#                           refUMAP = "umap",
#                           refDR = "pca",
#                           refAssay = "SCT",
#                           dims = 1:ndims,
#                           metadata = metadata_column,
#                           reference.version = "1.0.0",
#                           k.param = 31,
#                           ...)
#
#   ## Save azimuth compatible ref
#   Seurat::SaveAnnoyIndex(object = ref[["refdr.annoy.neighbors"]],
#                          file = file.path(output_folder,
#                                           "idx.annoy"))
#   saveRDS(object = ref,
#           file = file.path(output_folder, "ref.Rds"))
#   return(ref)
# }
#

## 3d UMAP ----

# library(threejs)
# obj.bcells <- RunUMAP(obj.bcells,
#                       reduction = 'harmony',
#                       reduction.name = 'test_3d',
#                       n.components = 3,
#                       dims = 1:20)
# obj.bcells <- obj.bcells %>%
#   FindNeighbors(reduction = 'harmony', dims = 1:20) %>%
#   FindClusters()
# scatterplot3js(obj.bcells@reductions$test_3d@cell.embeddings[,1],
#                obj.bcells@reductions$test_3d@cell.embeddings[,2],
#                obj.bcells@reductions$test_3d@cell.embeddings[,3],
#                color=plyr::mapvalues(obj.bcells$seurat_clusters,
#                                      from = levels(obj.bcells$seurat_clusters),
#                                      to = rainbow(12)),
#                size = 0.1,
#                labels = pbmc.obj$cluster_labels)
# ```

## WNN Processing via harmony

# ```{r}
# DefaultAssay(obj.tcells) <- 'RNA'
# obj.tcells <- obj.tcells %>%
#   NormalizeData() %>%
#   FindVariableFeatures(nfeatures = 2000 + length(blacklist_genes))
#
# VariableFeatures(obj.tcells) <- setdiff(VariableFeatures(obj.tcells),
#                                         blacklist_genes)
# VariableFeatures(obj.tcells) <- Seurat::VariableFeatures(obj.tcells)[1:2000]
#
#
# obj.tcells <- obj.tcells %>%
#   ScaleData() %>%
#   RunPCA(reduction.name = 'rna_pca')
# ElbowPlot(obj.tcells, reduction = 'rna_pca' , ndims = 50)
# ```
#
#
# ```{r}
# DefaultAssay(obj.tcells) <- 'ADT'
# obj.tcells <- obj.tcells %>%
#   NormalizeData(normalization.method = 'CLR', margin = 2)
# # old_varfeat <- VariableFeatures(obj.tcells)
# # VariableFeatures(obj.tcells) <- rownames(obj.tcells[['ADT']])
# obj.tcells <- obj.tcells %>%
#   ScaleData() %>%
#   RunPCA(reduction.name = 'adt_pca', features = rownames(obj.tcells[['ADT']]))
# ElbowPlot(obj.tcells, reduction = 'adt_pca')
# ```
#
# ```{r}
# DefaultAssay(obj.tcells) <- 'RNA'
# obj.tcells <- harmony::RunHarmony(obj.tcells,
#                                   group.by.vars = c('capID'),
#                                   reduction = 'rna_pca',
#                                   reduction.save = 'harmony_RNA',
#                                   dims = 1:28,
#                                   max.iter.harmony = 20,
#                                   plot_convergance = TRUE)
# ```
#
# ```{r}
# DefaultAssay(obj.tcells) <- 'ADT'
# obj.tcells <- harmony::RunHarmony(obj.tcells,
#                                   group.by.vars = c('capID'),
#                                   reduction = 'adt_pca',
#                                   reduction.save = 'harmony_ADT',
#                                   dims = 1:12)
# ```
#
# ```{r}
# ElbowPlot(obj.tcells, reduction = 'harmony_RNA', ndims = 50)
# ElbowPlot(obj.tcells, reduction = 'harmony_ADT', ndims = 50)
# ```
#
#
# ```{r}
# obj.tcells <- FindMultiModalNeighbors(obj.tcells,
#                                       reduction.list = list('harmony_RNA', 'harmony_ADT'),
#                                       dims.list = list(1:30, 1:11),
#                                       modality.weight.name = 'RNA.weight',
#                                       prune.SNN = 1/20)
# ```
#
#
# ```{r}
# obj.tcells <- RunUMAP(obj.tcells,
#                       nn.name = 'weighted.nn',
#                       reduction.name = 'wnn.umap',
#                       reduction.key = 'wnnUMAP_')
# ```
#
# ```{r}
# DimPlot(obj.tcells, reduction = 'wnn.umap', raster = FALSE, label = TRUE)
# ```
