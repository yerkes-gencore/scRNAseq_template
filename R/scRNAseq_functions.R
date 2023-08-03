# ## This file holds some small functions that may be useful, but haven't been
# ## robustly tested or added to packages, workflows, etc.
# ## if you find yourself frequently using one of these functions, we should
# ## try to improve the documentation of it and consider adding it
# ## to a package, or otherwise separating it from these other informal,
# ## possibly unused functions
#
## Derrik functions ----

### Cell migration table ----
# ## Example code to make table highlighting migration of cells from old
# ## clusters to new clusters

# cluster_migration <- table(bcell_obj$RNA_snn_res.1.4,bcell_obj$original_clustering_labels) %>%
#   as.data.frame.matrix()
# cluster_migration <- round(100* cluster_migration / rowSums(cluster_migration) ,2)
# cluster_migration <- cluster_migration[gtools::mixedsort(rownames(cluster_migration)),]
# cluster_migration <- cluster_migration[,gtools::mixedsort(colnames(cluster_migration))]
# cluster_migration_DT <-  datatable(
#   cluster_migration,
#   options = list(pageLength = 22,
#                  searching=FALSE,
#                  ordering=FALSE,
#                  lengthChange = FALSE)) %>%
#   DT::formatStyle(names(cluster_migration),
#                   background = styleColorBar(range(cluster_migration), 'lightblue'))
# cluster_migration_DT


### Extract values from GFF attributes ----
# extract_attributes <- function(gtf_attributes, att_of_interest, split_char = "; "){
#   att <- strsplit(gtf_attributes, split_char)
#   att <- gsub("\"","",unlist(att))
#   att <- str_trim(att, side="both")
#   if(!is.null(unlist(strsplit(att[grep(att_of_interest, att)], " ")))){
#     return( unlist(strsplit(att[grep(att_of_interest, att)], " "))[2])
#   }else{
#     return(NA)}
# }
#
### Bootstraping confidence intervals
# ## estimating cutting off extreme values
# ## 95% credible intervals for including 98% of data
# bootstrap_cutoffs <- function(capture, feature, k=1000, low_bound=0.01, high_bound=0.99){
#   boot_res <- matrix(data=NA, nrow=k, ncol=2)
#   for (i in 1:k){
#     x <- sample(capture@meta.data[[feature]], replace=TRUE, size=k)
#     boot_res[i,] <- quantile(x, probs=c(low_bound,high_bound))
#   }
#   return(list(low=quantile(boot_res[,1], probs=c(0.05)), high=quantile(boot_res[,2], probs=c(0.95))))
# }
#
#
### GO enrichment on cluster markers ----
# comp_clusters <- function(res, c1, c2){
#   tmp <- FindMarkers(combined.sct, group.by = paste0('integrated_snn_res.', res), ident.1=c1, ident.2=c2, logfc.threshold = 0.5, base=2) %>% mutate(TEST = p_val_adj<0.05)
#   tmp2 <- rownames(tmp %>% filter( p_val_adj<0.05))
#   go <- enrichGO(gene = tmp2,
#                  universe = rownames(combined.sct),
#                  keyType = "SYMBOL",
#                  OrgDb = org.Mm.eg.db,
#                  ont = "BP",
#                  pAdjustMethod = "BH",
#                  pvalueCutoff = 0.1,
#                  #qvalueCutoff = 1,
#                  readable=TRUE, pool=TRUE)
#   out <- go@result[,c("Description", "pvalue", "p.adjust", "geneID")]
#   return(list(table=out, genes=tmp))
# }
#
### histograms of most abundant genes across cells ----
# plot_top_genes <- function(data, n=20){
#   C <- data@assays$RNA@counts
#   C <- Matrix::t(Matrix::t(C)/Matrix::colSums(C)) * 100
#   most_expressed <- order(apply(C, 1, median), decreasing = T)[20:1]
#   boxplot(as.matrix(t(C[most_expressed, ])), cex = 0.1, las = 1, xlab = "% total count per cell", col = (scales::hue_pal())(20)[20:1], horizontal = TRUE, main=data@project.name)
# }
#
#
### median absolute deviation outlier detection ----
# is_outlier <- function(data, nmads){
#   data.mad <- mad(data)
#   outlier = (
#     (data < median(data) - nmads * data.mad) |
#       (median(data) + nmads * data.mad < data)
#   )
#   outlier
# }
#
# detect_outliers <- function(capture, nmads=4, hbc_format=TRUE){
#   ## Log transforming data to normalize it
#   if (hbc_format) {
#     outliers <- (is_outlier(log1p(capture$nUMI), nmads) |
#                    is_outlier(log1p(capture$nGene), nmads) |
#                    is_outlier(log1p(capture$log10GenesPerUMI), nmads))
#   } else {
#     outliers <- (is_outlier(log1p(capture$nCount_RNA), nmads) |
#                    is_outlier(log1p(capture$nFeature_RNA), nmads))
#   }
#   capture$outlier <- outliers
#   capture
# }
#
#
## Micah functions ----
#
# # This is a hacky way to assign more than one output from a function to different objects
# # This saves memory and cleans up code when returning a large seurat obj along with plots
# # Shamelessly copied verbatim from code by Gabor Grothendieck at: https://stat.ethz.ch/pipermail/r-help/2004-June/053343.html
# list <- structure(NA,class="result")
# "[<-.result" <- function(x,...,value) {
#   args <- as.list(match.call())
#   args <- args[-c(1:2,length(args))]
#   length(value) <- length(args)
#   for(i in seq(along=args)) {
#     a <- args[[i]]
#     if(!missing(a)) eval.parent(substitute(a <- v,list(a=a,v=value[[i]])))
#   }
#   x
# }
#
#
# quickSCTandPCA <- function(seurat_obj, regress.out = NULL, npcs = 20) {
#   seurat_obj <-
#     SCTransform(seurat_obj, vst.flavor = "v2", verbose = FALSE,
#                 vars.to.regress = regress.out)
#   seurat_obj <- RunPCA(seurat_obj, npcs = npcs, verbose = FALSE)
#   return(seurat_obj)
# }
#
# # should plot elbow plot between these two
# quickCluster <- function(seurat_obj, plot_title=NULL,
#                          regress.out = NULL, npcs = 20, res = 0.4) {
#   ## Basic normalization and clustering -> save obj to file
#   seurat_obj <- RunUMAP(seurat_obj, reduction = "pca", dims = 1:npcs, verbose=FALSE)
#   seurat_obj <- FindNeighbors(seurat_obj, reduction = "pca", dims = 1:npcs, verbose=FALSE)
#   seurat_obj <- FindClusters(seurat_obj, resolution = res, verbose=FALSE)
#
#   DimPlotTitle <- paste(as.character(unique(obj_capture$capID)), plot_title)
#   p1 <- DimPlot(seurat_obj, label = T, repel = T) + ggtitle(DimPlotTitle)
#   #ggsave(here("plots", "UMAP.png"), p1, height=8, width=10)
#
#   return(list(data = seurat_obj, plots = p1))
# }
#
