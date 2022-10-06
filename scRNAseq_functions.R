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

## plot ncountRNAxpercent.mito and ncountRNAxnFeatureRNA in the same plot with lines for possible cutoffs
qc_vln_plots <- function(obj, mt.cutoff = 5, nfeature.cutoff.low = 400, nfeature.cutoff.high = 6000 ) {
  plot1 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.mito") +
    geom_hline(yintercept = mt.cutoff, lty='dashed') + aes(alpha=0.02) +
    theme(legend.position="none")
  plot2 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
    geom_hline(yintercept=nfeature.cutoff.low, lty='dashed') +
    geom_hline(yintercept=nfeature.cutoff.high, lty='dashed')
  plot1 + plot2 
}

## Grabs QC files for samples and combines them in datatable
capture_QC_metrics <- function(samples){
  files <- lapply(samples$FileID, function(x) {read.csv(Sys.glob(file.path(config$rootDir, config$alignmentDir, x,"outs/metrics_summary.csv")))})
  qc <- t(data.table::rbindlist(files))
  colnames(qc) <- samples$Label
  rownames(qc) <- unlist(lapply(rownames(qc), gsub, pattern=".", replacement=" ", fixed=TRUE))
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
