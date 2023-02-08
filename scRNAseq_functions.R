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

### Make captures from samplesheet
generate_captures <- function(samplesheet, min.cells=10, min.features=200, mt_pattern="^mt-", rb_pattern= "^Rp[sl]"){
  files <- unlist(lapply(samplesheet$FileID, function(x) {Sys.glob(file.path(config$rootDir, config$alignmentDir, x,"outs/filtered_feature_bc_matrix"))}))
  counts <- lapply(files, Read10X)
  a <- mapply(create_from_10X, counts, project=samplesheet$Label, min.cells=10, min.features=200, mt_pattern="^mt-", rb_pattern= "^Rp[sl]")
  names(a) <- samplesheet$Label
  a
}

## wrapper for create seurat object, collects some extra metadata
create_from_10X <- function(data, project, min.cells=10, min.features=200, mt_pattern="^mt-", rb_pattern= "^Rp[sl]"){
  obj <- CreateSeuratObject(counts=data, project=project, min.cells = min.cells, min.features=min.features)
  
  # Get Metadata
  DefaultAssay(obj) <- "RNA"
  # obj[["percent.mito"]] <- PercentageFeatureSet(obj, features = mt.genes[mt.genes %in% row.names(obj)])
  obj[["percent.mito"]] <- PercentageFeatureSet(obj, pattern = mt_pattern)
  obj[["percent_ribo"]] = PercentageFeatureSet(obj, pattern = rb_pattern)
  obj$log10GenesPerUMI <- log10(obj$nFeature_RNA) / log10(obj$nCount_RNA)
  obj
}

## create hashed 10x captures
create_hashed_objects <- function(file, label, hash_table, mito.pattern = "^MT", ribo.pattern="^RP[SL]") {
  data <- Read10X_h5(Sys.glob(file.path(config$rootDir, config$alignmentDir, file,"outs/per_sample_outs/",file,"count/sample_filtered_feature_bc_matrix.h5")))
  unhashed <- CreateSeuratObject(counts=data$`Gene Expression`, project = label)
  unhashed[["percent.mito"]] <- PercentageFeatureSet(unhashed, pattern=mito.pattern)
  unhashed[["percent.ribo"]] <- PercentageFeatureSet(unhashed, pattern=ribo.pattern)
  unhashed <- NormalizeData(unhashed)
  unhashed <- FindVariableFeatures(unhashed, selection.method = "mean.var.plot")
  unhashed <- ScaleData(unhashed, features=VariableFeatures(unhashed))
  adt_data <- data$`Antibody Capture`[rownames(data$`Antibody Capture`)[!grepl("Hash", rownames(data$`Antibody Capture`), ignore.case = TRUE)],]
  unhashed[['ADT']] <- CreateAssayObject(counts = adt_data)
  unhashed <- NormalizeData(unhashed, assay = "ADT", normalization.method = "CLR", margin = 2)
  hto_data <- data$`Antibody Capture`[rownames(data$`Antibody Capture`)[grepl("Hash", rownames(data$`Antibody Capture`), ignore.case = TRUE)],]
  rownames(hto_data) <- plyr::mapvalues(substr(rownames(hto_data), 1, 10), from=hash_table$HTO, to=hash_table$Label)
  #rownames(hto_data) <- paste("hash", 1:nrow(hto_data), sep="_")
  unhashed[['HTO']] <- CreateAssayObject(counts=hto_data)
  unhashed <- NormalizeData(unhashed, assay='HTO', normalization.method = 'CLR')
  HTODemux(unhashed, assay='HTO', positive.quantile = 0.99)
}

## plot ncountRNAxpercent.mito and ncountRNAxnFeatureRNA in the same plot with lines for possible cutoffs
qc_plots_cutoffs <- function(obj, cutoffs, show_outliers=FALSE) {
  project <- obj@project.name
  plot1 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.mito",
                          group.by = (if (show_outliers) 'outlier' else NULL)) +
    theme(legend.position=(if (show_outliers) 'right' else 'none')) +
    labs(color = 'Outlier') +
    geom_hline(yintercept=cutoffs$percent.mito.max, lty='dashed') +
    geom_vline(xintercept=cutoffs$nCount_RNA.max, lty='dashed') + 
    geom_vline(xintercept=cutoffs$nCount_RNA.min, lty='dashed') +
    ggtitle(project)
  plot2 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
                          group.by = (if (show_outliers) 'outlier' else NULL)) +
    theme(legend.position=(if (show_outliers) 'right' else 'none')) +
    labs(color = 'Outlier') +
    geom_hline(yintercept=cutoffs$nFeature_RNA.max, lty='dashed') +
    geom_hline(yintercept=cutoffs$nFeature_RNA.min, lty='dashed') + 
    geom_vline(xintercept=cutoffs$nCount_RNA.max, lty='dashed') + 
    geom_vline(xintercept=cutoffs$nCount_RNA.min, lty='dashed') +
    ggtitle(project)
  print(plot1 + plot2)
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

qc_plots <- function(cap, show_outliers=FALSE){
  plot1 <- FeatureScatter(cap, feature1 = "nCount_RNA", feature2 = "percent.mito",
                          group.by = (if (show_outliers) 'outlier' else NULL))  +
    theme(legend.position=(if (show_outliers) 'right' else 'none')) +
    labs(color = 'Outlier') +
    ggtitle(cap@project.name)
  plot2 <- FeatureScatter(cap, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
                          group.by = (if (show_outliers) 'outlier' else NULL)) +
    theme(legend.position=(if (show_outliers) 'right' else 'none')) +
    labs(color = 'Outlier') +
    ggtitle(cap@project.name)
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

## Runs doublet finder, with or without hashing. Latent doublet rate
## may need to be changed depending on library prep
doublet_finder <- function(obj, pcs, hashed,
                           cores = 1, doublet_rate = 0.075){
  sweep.res.list <- paramSweep_v3(obj, PCs = 1:pcs, sct = TRUE,
                                  num.cores = cores)
  if (hashed) {
    sweep.stats <- summarizeSweep(sweep.res.list, GT = TRUE, 
                                  GT.calls = obj$HTO_classification.global)
  } else {
    sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  }
  bcmvn <- find.pK(sweep.stats)
  pK <- bcmvn[bcmvn$BCmetric==max(bcmvn$BCmetric),'pK']
  annotations <- obj@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations) ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
  nExp_poi <- round(doublet_rate*nrow(obj@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  results <- doubletFinder_v3(obj, PCs = 1:pcs, pK = as.numeric(as.character(pK)), nExp = nExp_poi.adj, sct = TRUE)
  return(results)
}


## Read GEX data and correct ambient RNA with soupX
read_counts_soupX <- function(label, file, mito.pattern="^MT-") {
  filtered_data <- Read10X(
    dir(file.path(config$rootDir, config$alignmentDir, file, 'outs'),
        recursive = TRUE, pattern = 'filtered_feature_bc_matrix$',
        include.dirs = TRUE, full.names = TRUE))
  if (!is.list(filtered_data)) {
    filtered_data <- list(`Gene Expression` = filtered_data)
  }
  
  clus <- read.csv(file.path(dir(
    file.path(config$rootDir, config$alignmentDir, file, 'outs'),
    recursive = TRUE, pattern = 'clustering', include.dirs = TRUE, full.names = TRUE),
    'gene_expression_graphclust/clusters.csv'))
  clusters <- clus$Cluster
  names(clusters) <- clus$Barcode
  
  unfiltered_data <- Read10X(dir(file.path(
    config$rootDir, config$alignmentDir, file, 'outs'),
    recursive = TRUE, pattern = 'raw_feature_bc_matrix$',
    include.dirs = TRUE, full.names = TRUE))
  if (!is.list(unfiltered_data)) {
    unfiltered_data <- list(`Gene Expression` = unfiltered_data)
  }
  
  sc <- SoupChannel(tod = unfiltered_data$`Gene Expression`,
                    toc = filtered_data$`Gene Expression`)
  sc <- setClusters(sc, clusters)
  sc <- autoEstCont(sc)
  fixed_counts <- adjustCounts(sc, roundToInt = TRUE)
  compare <- 1- (rowSums(fixed_counts) / rowSums(filtered_data$`Gene Expression`)) 
  count_diff <- data.frame(diff=compare, total=rowSums(filtered_data$`Gene Expression`)) %>%
    rownames_to_column('Gene')
  soupX_qc_plot_1 <- 
    ggplot(count_diff,
           aes(x=diff, y=total, label=Gene, color = grepl(mito.pattern, Gene))) +
    scale_color_manual(name = 'Mitochondrial genes',
                       values = setNames(c('green', 'grey'), c(TRUE,FALSE))) + 
    geom_point(alpha=0.5) + 
    geom_text_repel(data=subset(count_diff, diff > 0.25 | total > 1000000),
                    position=position_jitter()) + 
    scale_y_log10() + 
    labs(x='Portion of reads removed', y='Original count total', title = label) +
    theme_bw()
  message('SoupX diagnostic plots saved in "plots/"')
  ggsave(paste0('plots/', label, '_QC1_soupX.png'), plot = soupX_qc_plot_1)
  
  soupX_qc_plot_2 <- 
    ggplot(count_diff, aes(x=diff)) + 
    geom_histogram() +
    theme_bw() + 
    labs(x='Proportion of reads removed', y='Frequency',
         title = 'Portion of reads removed per sample')
  ggsave(paste0('plots/', label, '_QC2_soupX.png'), plot = soupX_qc_plot_2)
  filtered_data$`Gene Expression` <- fixed_counts
  filtered_data
}

## Create seurat object from counts (ideally the corrected ones from SoupX)
create_and_preprocess_Seurat_object_from_counts <- function(
    counts, project_label,
    mito.pattern = "^MT", ribo.pattern="^RP[SL]",
    min.cells=10, min.features=0,
    s.features, g2m.features) {
  
  obj <- CreateSeuratObject(counts=counts$`Gene Expression`, project = project_label,
                            min.cells = min.cells, min.features=min.features)
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj, selection.method = "mean.var.plot")
  obj <- ScaleData(obj, features=VariableFeatures(obj))
  obj <- CellCycleScoring(obj, s.features=s.features, g2m.features = g2m.features)
  obj[["percent.mito"]] <- PercentageFeatureSet(obj, pattern=mito.pattern)
  obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern=ribo.pattern)
  obj
}

## Add hash data for Antibody capture assay, using counts object from read_counts_soupX
add_hash_data <- function(obj, counts, hash_table){
  hto_data <- counts$`Antibody Capture`[
    rownames(counts$`Antibody Capture`)[
      grepl("Hash", rownames(counts$`Antibody Capture`), ignore.case = TRUE)
    ],
  ]
  rownames(hto_data) <- plyr::mapvalues(substr(rownames(hto_data), 1, 10),
                                        from=hash_table$HTO, to=hash_table$Label)
  obj[['HTO']] <- CreateAssayObject(counts=hto_data)
  obj <- NormalizeData(obj, assay='HTO', normalization.method = 'CLR')
  HTODemux(obj, assay='HTO', positive.quantile = 0.99)
}