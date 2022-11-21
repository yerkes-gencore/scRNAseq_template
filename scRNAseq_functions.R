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
qc_vln_plots <- function(obj, mt.cutoff = 5, nfeature.cutoff.low = 400, nfeature.cutoff.high = 6000,
                         ncount.low=1000, ncount.high=20000, bootstrap=TRUE) {
  if (bootstrap) {
    nfeature.cutoff <- bootstrap_cutoffs(obj, feature='nFeature_RNA')
    nfeature.cutoff.low <- nfeature.cutoff$low
    nfeature.cutoff.high <- nfeature.cutoff$high
    mt.cutoff <- bootstrap_cutoffs(obj, feature='percent.mito')$high
    ncount <- bootstrap_cutoffs(obj, feature='nCount_RNA')
    ncount.low <- ncount$low
    ncount.high <- ncount$high
  }
  plot1 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.mito") +
    geom_hline(yintercept = mt.cutoff, lty='dashed') +
    theme(legend.position="none") + 
    geom_vline(xintercept=ncount.high, lty='dashed') + 
    geom_vline(xintercept=ncount.low, lty='dashed')
  plot2 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
    geom_hline(yintercept=nfeature.cutoff.low, lty='dashed') +
    geom_hline(yintercept=nfeature.cutoff.high, lty='dashed') + 
    geom_vline(xintercept=ncount.high, lty='dashed') + 
    geom_vline(xintercept=ncount.low, lty='dashed')
  print(plot1 + plot2)
  return(data.frame(nfeature=c(nfeature.cutoff.low,nfeature.cutoff.high), ncount=c(ncount.low,ncount.high), mito=c(0, mt.cutoff), row.names = c('low', 'high')))
  
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


