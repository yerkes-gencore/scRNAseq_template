## For alternatives to what's in the template scripts.
## These are unlikely to be used frequently, if at all, but
## we can preserve the code here in case we wish to reuse them

## Doublet finder ----

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
