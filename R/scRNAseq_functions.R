#### Extract values from GFF attributes
extract_attributes <- function(gtf_attributes, att_of_interest, split_char = "; "){
  att <- strsplit(gtf_attributes, split_char)
  att <- gsub("\"","",unlist(att))
  att <- str_trim(att, side="both")
  if(!is.null(unlist(strsplit(att[grep(att_of_interest, att)], " ")))){
    return( unlist(strsplit(att[grep(att_of_interest, att)], " "))[2])
  }else{
    return(NA)}
}

## estimating cutting off extreme values 
## 95% credible intervals for including 98% of data
bootstrap_cutoffs <- function(capture, feature, k=1000, low_bound=0.01, high_bound=0.99){
  boot_res <- matrix(data=NA, nrow=k, ncol=2)
  for (i in 1:k){
    x <- sample(capture@meta.data[[feature]], replace=TRUE, size=k)
    boot_res[i,] <- quantile(x, probs=c(low_bound,high_bound))
  }
  return(list(low=quantile(boot_res[,1], probs=c(0.05)), high=quantile(boot_res[,2], probs=c(0.95))))
}

## Grabs QC files for samples and combines them in datatable
capture_QC_metrics <- function(samples){
  files <- lapply(samples$FileID, function(x) {read.csv(Sys.glob(file.path(config$rootDir, config$alignmentDir, x,"outs/metrics_summary.csv")))})
  qc <- t(data.table::rbindlist(files))
  colnames(qc) <- samples$Label
  rownames(qc) <- unlist(lapply(rownames(qc), gsub, pattern=".", replacement=" ", fixed=TRUE))
  qc
}

## grabs QC for hashed data
capture_QC_metrics_hashed <- function(samples){
  files <- lapply(samples$FileID, function(x) {read.csv(Sys.glob(file.path(config$rootDir, config$alignmentDir, x,"outs/per_sample_outs/",x,"metrics_summary.csv")))})
  qc <- data.table::rbindlist(files, use.names=TRUE, fill=TRUE,idcol="ID")
  qc$ID <- plyr::mapvalues(qc$ID, from=1:length(samples$Label), to=samples$Label)
  #rownames(qc) <- unlist(lapply(rownames(qc), gsub, pattern=".", replacement=" ", fixed=TRUE))
  qc
}

comp_clusters <- function(res, c1, c2){
  tmp <- FindMarkers(combined.sct, group.by = paste0('integrated_snn_res.', res), ident.1=c1, ident.2=c2, logfc.threshold = 0.5, base=2) %>% mutate(TEST = p_val_adj<0.05)
  tmp2 <- rownames(tmp %>% filter( p_val_adj<0.05))
  go <- enrichGO(gene = tmp2,
                 universe = rownames(combined.sct),
                 keyType = "SYMBOL",
                 OrgDb = org.Mm.eg.db,
                 ont = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.1,
                 #qvalueCutoff = 1,
                 readable=TRUE, pool=TRUE)
  out <- go@result[,c("Description", "pvalue", "p.adjust", "geneID")]
  return(list(table=out, genes=tmp))
}

## histograms of most abundant genes across cells
plot_top_genes <- function(data, n=20){
  C <- data@assays$RNA@counts
  C <- Matrix::t(Matrix::t(C)/Matrix::colSums(C)) * 100
  most_expressed <- order(apply(C, 1, median), decreasing = T)[20:1]
  boxplot(as.matrix(t(C[most_expressed, ])), cex = 0.1, las = 1, xlab = "% total count per cell", col = (scales::hue_pal())(20)[20:1], horizontal = TRUE, main=data@project.name) 
}

qc_plots_hash <- function(cap){
  plot1 <- FeatureScatter(cap, feature1 = "nCount_RNA", feature2 = "percent.mito", group.by = "HTO_classification.global")  +
    theme(legend.position="none") + ggtitle(cap@project.name)
  plot2 <- FeatureScatter(cap, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "HTO_classification.global") +   ggtitle(cap@project.name)
  plot1 + plot2 
}

### Outlier detection for a vector of values, based on median absolute deviation
is_outlier <- function(data, nmads){
  data.mad <- mad(data)
  outlier = (
    (data < median(data) - nmads * data.mad) |
      (median(data) + nmads * data.mad < data)
  )
  outlier
}
detect_outliers <- function(capture, nmads=4){
  ## Log transforming data to normalize it
  outliers <- (is_outlier(log1p(capture$nCount_RNA), nmads) |
                 is_outlier(log1p(capture$nFeature_RNA), nmads))
  capture$outlier <- outliers
  capture
}

## Returns top labels from map query result for each cluster 
## including  median score and proportion of cluster
mapquery_top_results <- function(obj, top_n, label_column, score_column){
  obj@meta.data %>%
    dplyr::select(seurat_clusters, {{ label_column }}, {{ score_column }}) %>%
    group_by(seurat_clusters, pick({{ label_column }})) %>% 
    summarise(median_score = median(.data[[ score_column ]]),
              n = n(), .groups = 'drop_last') %>% 
    mutate(freq = n / sum(n)) %>%
    top_n(2) %>% 
    dplyr::select(-n) %>% 
    arrange(seurat_clusters, desc(freq))
}
