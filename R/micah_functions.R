## Functions for basic pre-processing of scRNA
## Using this to compare different cell calling and ambient reads handling methods
## on Abbie_10X data (ignoring ADT data for now)

# This is a hacky way to assign more than one output from a function to different objects
# This saves memory and cleans up code when returning a large seurat obj along with plots
# Shamelessly copied verbatim from code by Gabor Grothendieck at: https://stat.ethz.ch/pipermail/r-help/2004-June/053343.html
list <- structure(NA,class="result")
"[<-.result" <- function(x,...,value) {
  args <- as.list(match.call())
  args <- args[-c(1:2,length(args))]
  length(value) <- length(args)
  for(i in seq(along=args)) {
    a <- args[[i]]
    if(!missing(a)) eval.parent(substitute(a <- v,list(a=a,v=value[[i]])))
  }
  x
}

## For each capture: 
## read in, calculate basic stats and run basic clustering
read_counts <- function(capID, file, mito.pattern = "^MT", ribo.pattern = "^RP[SL]") {
  ## read in counts
  seurat_data <- Seurat::Read10X(data.dir = file)
  ## create Seurat object
  seurat_obj <- CreateSeuratObject(counts = seurat_data$`Gene Expression`,
                                   min.features = 100, # necessary for sctransform::vst but may not be desirable for some initial diagnostics
                                   min.cells = 1, # necessary for cell cycle scoring (and possibly other things)
                                   project = capID)
  show(seurat_obj)
  
  ## calculate qc metrics -> save metadata in memory
  # Add number of genes per UMI for each cell to metadata
  seurat_obj$log10GenesPerUMI <- log10(seurat_obj$nFeature_RNA) / log10(seurat_obj$nCount_RNA)
  # Mitochondrial ratio
  seurat_obj$mitoRatio <- PercentageFeatureSet(object = seurat_obj, pattern = mito.pattern)/100
  # Ribosomal ratio
  seurat_obj$riboRatio <- PercentageFeatureSet(object = seurat_obj, pattern = ribo.pattern)/100
  # Create metadata dataframe
  metadata <- seurat_obj@meta.data
  # Add cell IDs to metadata
  metadata$cells <- rownames(metadata)
  # Rename columns to be more readable/intuitive
  metadata <- metadata %>%
    dplyr::rename(capID = orig.ident,
                  nUMI = nCount_RNA,
                  nGene = nFeature_RNA)
  seurat_obj@meta.data <- metadata
  return(seurat_obj)
}

## Check for effect of cell cycle and mitoRatio on PCA
checkCC <- function(seurat_obj, plot_title=NULL) {
  # A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
  # segregate this list into markers of G2/M phase and markers of S phase
  s_genes <- cc.genes$s.genes
  g2m_genes <- cc.genes$g2m.genes
  # Score cells for cell cycle
  seurat_phase <- NormalizeData(seurat_obj, verbose=FALSE)
  seurat_phase <- CellCycleScoring(seurat_phase, verbose=FALSE,
                                   g2m.features = g2m_genes, 
                                   s.features = s_genes)
  seurat_phase <- FindVariableFeatures(seurat_phase, verbose=FALSE,
                                       selection.method = "vst",
                                       nfeatures = 2000)
  seurat_phase <- ScaleData(seurat_phase, verbose=FALSE)
  seurat_phase <- RunPCA(seurat_phase, npcs=20, verbose=FALSE)
  
  DimPlotTitle <- paste(as.character(unique(obj_capture$capID)), plot_title)
  
  # Plot the PCA colored by cell cycle phase
  (p_cc <- DimPlot(seurat_phase,
                   reduction = "pca",
                   group.by= "Phase",
                   split.by = "Phase") +
      ggtitle(DimPlotTitle))
  #ggsave(here("plots/cc", "cc.png"), p_cc, height=5, width=10)
  # To Do: save plot as png with reasonable filename
  
  # Turn mitoRatio into categorical factor vector based on quartile values
  mito_qs <- summary(seurat_phase@meta.data$mitoRatio)
  seurat_phase@meta.data$mitoFr <- cut(seurat_phase@meta.data$mitoRatio, 
                                       breaks=c(-Inf, mito_qs[2], mito_qs[3], mito_qs[4], Inf), 
                                       labels=c("Low","Medium","Medium high", "High"))
  # Plot the PCA colored by mito cat
  (p_mt <- DimPlot(seurat_phase,
                   reduction = "pca",
                   group.by= "mitoFr",
                   split.by = "mitoFr") +
      ggtitle(DimPlotTitle))
  #ggsave(here("plots", "mt.png"), p_mt, height=5, width=10)
  
  # Turn mitoRatio into categorical factor vector based on quartile values
  ribo_qs <- summary(seurat_phase@meta.data$riboRatio)
  seurat_phase@meta.data$riboFr <- cut(seurat_phase@meta.data$riboRatio, 
                                       breaks=c(-Inf, ribo_qs[2], ribo_qs[3], ribo_qs[4], Inf), 
                                       labels=c("Low","Medium","Medium high", "High"))
  # Plot the PCA colored by mito cat
  (p_rb <- DimPlot(seurat_phase,
                   reduction = "pca",
                   group.by= "riboFr",
                   split.by = "riboFr") +
      ggtitle(DimPlotTitle))
  #ggsave(here("plots", "rb.png"), p_rb, height=5, width=10)
  
  plots <- list(pca_cc=p_cc, pca_mt=p_mt, pca_rb=p_rb)
  return(list(data = seurat_phase,
              plots = plots))
}

quickSCTandPCA <- function(seurat_obj, regress.out = NULL, npcs = 20) {
  seurat_obj <- 
    SCTransform(seurat_obj, vst.flavor = "v2", verbose = FALSE,
                vars.to.regress = regress.out)
  seurat_obj <- RunPCA(seurat_obj, npcs = npcs, verbose = FALSE)
  return(seurat_obj)
}

# should plot elbow plot between these two
quickCluster <- function(seurat_obj, plot_title=NULL,
                         regress.out = NULL, npcs = 20, res = 0.4) {
  ## Basic normalization and clustering -> save obj to file
  seurat_obj <- RunUMAP(seurat_obj, reduction = "pca", dims = 1:npcs, verbose=FALSE)
  seurat_obj <- FindNeighbors(seurat_obj, reduction = "pca", dims = 1:npcs, verbose=FALSE)
  seurat_obj <- FindClusters(seurat_obj, resolution = res, verbose=FALSE)
  
  DimPlotTitle <- paste(as.character(unique(obj_capture$capID)), plot_title)
  p1 <- DimPlot(seurat_obj, label = T, repel = T) + ggtitle(DimPlotTitle)
  #ggsave(here("plots", "UMAP.png"), p1, height=8, width=10)
  
  return(list(data = seurat_obj, plots = p1))
}

# Cell-level metrics
defaultCutoffs <- list(mitoRatio.max = Inf,
                       mitoRatio.min = 0.2,
                       nUMI.max = Inf,
                       nUMI.min = 500,
                       nGene.max = Inf,
                       nGene.min = 250,
                       log10GenesPerUMI.max = Inf,
                       log10GenesPerUMI.min = 0.80)

### Outlier detection for a vector of values, based on median absolute deviation
is_outlier <- function(data, nmads){
  data.mad <- mad(data)
  outlier = (
    (data < median(data) - nmads * data.mad) |
      (median(data) + nmads * data.mad < data)
  )
  outlier
}

detect_outliers <- function(capture, nmads=4, hbc_format=TRUE){
  ## Log transforming data to normalize it
  if (hbc_format) {
    outliers <- (is_outlier(log1p(capture$nUMI), nmads) |
                   is_outlier(log1p(capture$nGene), nmads) |
                   is_outlier(log1p(capture$log10GenesPerUMI), nmads))
  } else {
    outliers <- (is_outlier(log1p(capture$nCount_RNA), nmads) |
                   is_outlier(log1p(capture$nFeature_RNA), nmads))
  }
  capture$outlier <- outliers
  capture
}

## Filter based on arbitrary cutoffs and/or outliers (outliers may not be working!)
filterThresholdsAndOutliers <- 
  function(capture, cutoffs = defaultCutoffs, keep_all = FALSE,
           by_outlier = FALSE, nmads = 4, filterName = "filtThreshOnly",
           by_threshold = TRUE, hbc_format = TRUE, gene_level = TRUE) {
    filt <- capture
    if (by_outlier & !by_threshold) {
      print("Ignoring 'filterName' argument because only filtering by outliers; logical stored in 'outliers' metadata column.")
    } else {
      print(paste0("Creating new metadata column to store filter logical: ", filterName))
    }
    
    ## Cell-level filtering
    # Filter by outlier detection
    if (!by_outlier) {
      filt$outlier <- FALSE
    } else {
      filt <- detect_outliers(filt, hbc_format=hbc_format, nmads)
      filt_subsetted <- subset(filt, subset = (outlier == FALSE)) # in case you only want to filter by outliers; will be overwritten otherwise
    }
    
    # Filter by manually specified cutoffs
    if(by_threshold) {
      if (hbc_format) {
        #print(cutoffs)
        metadata <- filt@meta.data
        # Create new column with cell-level filtering logical
        metadata <- metadata %>% mutate(
          !!filterName :=
            nUMI >= cutoffs$nUMI.min & nUMI <= cutoffs$nUMI.max &
            nGene >= cutoffs$nGene.min & nGene <= cutoffs$nGene.max &
            log10GenesPerUMI >= cutoffs$log10GenesPerUMI.min & log10GenesPerUMI <= cutoffs$log10GenesPerUMI.max &
            mitoRatio >= cutoffs$mitoRatio.min & mitoRatio <= cutoffs$mitoRatio.max &
            outlier == FALSE)
        filt@meta.data <- metadata
        # Create a new subsetted object
        filt_subsetted <- subset(filt, subset =
                                   nUMI >= cutoffs$nUMI.min & nUMI <= cutoffs$nUMI.max &
                                   nGene >= cutoffs$nGene.min & nGene <= cutoffs$nGene.max &
                                   log10GenesPerUMI >= cutoffs$log10GenesPerUMI.min & log10GenesPerUMI <= cutoffs$log10GenesPerUMI.max &
                                   mitoRatio >= cutoffs$mitoRatio.min & mitoRatio <= cutoffs$mitoRatio.max &
                                   outlier == FALSE)
        # Print the counts of filtered and unfiltered cells
        print(filt@meta.data %>% dplyr::count(get(filterName)))
        
      } else {
        print(cutoffs)
        # Create new column with cell-level filtering logical
        metadata <- metadata %>% mutate(
          !!filterName :=
            nFeature_RNA > cutoffs$nFeature_RNA.min & 
            nFeature_RNA < cutoffs$nFeature_RNA.max &
            percent.mito < cutoffs$percent.mito.max & 
            nCount_RNA   > cutoffs$nCount_RNA.min & 
            nCount_RNA   < cutoffs$nCount_RNA.max &
            outlier == FALSE)
        filt@meta.data <- metadata
        # Create a new subsetted object
        filt_subsetted <- subset(filt, subset =
                                   nUMI >= cutoffs$nUMI.min & nUMI <= cutoffs$nUMI.max &
                                   nGene >= cutoffs$nGene.min & nGene <= cutoffs$nGene.max &
                                   log10GenesPerUMI >= cutoffs$log10GenesPerUMI.min & log10GenesPerUMI <= cutoffs$log10GenesPerUMI.max &
                                   mitoRatio >= cutoffs$mitoRatio.min & mitoRatio <= cutoffs$mitoRatio.max &
                                   outlier == FALSE)
        # Print the counts of filtered and unfiltered cells
        print(filt@meta.data %>% dplyr::count(get(filterName)))
      }
    }
    
    ## Gene-level filtering
    if (gene_level) {
      print("Only keeping genes which are expressed in 10 or more cells")
      counts <- GetAssayData(object = filt_subsetted, slot = "counts")
      nonzero <- counts > 0
      # Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
      keep_genes <- Matrix::rowSums(nonzero) >= 10
      filtered_counts <- counts[keep_genes, ]
      filt_subsetted <- CreateSeuratObject(filtered_counts, meta.data = filt@meta.data)
    }
    print("Note: Gene-level filtering not retained in unsubsetted metadata. Use subsetted if gene-level filtering is needed.")
    
    # Save summaries of filtered genes in a table and return
    filtSummary <- list(counts = list(nCellsIn = ncol(capture),
                                      nCellsRemoved = ncol(capture) - ncol(filt_subsetted),
                                      nCellsOut = ncol(filt_subsetted),
                                      nGenesIn = nrow(counts),
                                      nGenesRemoved = sum(Matrix::rowSums(nonzero) < 10),
                                      nGenesOut = sum(keep_genes)),
                        names = list(CellsIn = colnames(capture),
                                     CellsRemoved = colnames(capture)[!(colnames(capture)%in%colnames(filt_subsetted))],
                                     CellsOut = colnames(filt_subsetted),
                                     GenesIn = rownames(counts),
                                     GenesRemoved = colnames(counts)[!(colnames(counts)%in%colnames(keep_genes))],
                                     GenesOut = colnames(keep_genes)))
    
    print(paste0('Removing ', filtSummary$counts$nCellsRemoved, ' cells of ', filtSummary$counts$nCellsIn, ' from capture ', capture@project.name, ' due to quality filtering cutoffs'))
    print(paste0(filtSummary$counts$nGenesOut, " genes remain after removing ", filtSummary$counts$nGenesRemoved, " genes with nonzero counts in less than 10 cells."))
    return(list(all_data = filt, subset = filt_subsetted, filtSummary = filtSummary))
  }

## QC plotting functions (for multiple captures at once)
plotQC_ridges <- function(metadata, cutoffs = defaultCutoffs) {
  # Visualize the number of cell counts per sample
  qc_p1 <- metadata %>% 
    ggplot(aes(x=capID, group=capID, fill=capID)) + 
    geom_bar() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    theme(plot.title = element_text(hjust=0.5, face="bold")) +
    ggtitle("N Cells") +
    theme(legend.position="none")
  
  # Visualize the distribution of genes detected per cell via histogram
  qc_p2 <- metadata %>% 
    ggplot(aes(x=nGene, y=capID, fill=capID, color=capID)) + 
    geom_density_ridges() + 
    scale_x_log10() + 
    theme_classic() +
    geom_vline(xintercept = cutoffs$nGene.min) +
    theme(legend.position="none") +
    ggtitle("N Genes")
  
  # Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI (novelty score)
  qc_p3 <- metadata %>% 
    ggplot(aes(x=log10GenesPerUMI, y=capID, fill=capID, color=capID)) + 
    geom_density_ridges() + 
    theme_classic() +
    geom_vline(xintercept = cutoffs$log10GenesPerUMI.min) +
    theme(legend.position="none") +
    ggtitle("Complexity Score")
  
  # Visualize the distribution of mitochondrial gene expression detected per cell
  qc_p4 <- metadata %>% 
    ggplot(aes(x=mitoRatio, y=capID, fill=capID, color=capID)) + 
    geom_density_ridges() + 
    theme_classic() +
    geom_vline(xintercept = cutoffs$mitoRatio.min) +
    theme(legend.position="none") +
    ggtitle("Mitochondrial ratio")
  
  allCaps_QC_ridges <- ggarrange(qc_p1, qc_p2, qc_p3, qc_p4, ncol=2, nrow=2)
  return(allCaps_QC_ridges)
  # To do: distribution of reads per cell
  # Need to read in cellranger summary output for this
}

## Joint filtering effects
plotQC_joint <- function(metadata, cutoffs,
                         nUMI_cutoff=500, nGene_cutoff=250, 
                         mitoRatio_cutoff=0.2,
                         color) {
  # Visualize the correlation between genes detected and number of UMIs and 
  # determine whether strong presence of cells with low numbers of genes/UMIs
  if (color == "mitoRatio") {
    allCaps_QC_joint <- metadata %>% 
      ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
      geom_point(alpha=0.2) + 
      scale_colour_gradient(low = "goldenrod1", high = "red", limits=c(0,1)) +
      stat_smooth(method=lm,) +
      scale_x_log10() + 
      scale_y_log10() + 
      theme_classic() +
      geom_vline(xintercept = cutoffs$nUMI.min) +
      geom_hline(yintercept = cutoffs$nGene.min) +
      facet_wrap(~capID)
  } else if (color == "riboRatio") {
    allCaps_QC_joint <- metadata %>% 
      ggplot(aes(x=nUMI, y=nGene, color=riboRatio)) + 
      geom_point(alpha=0.2) + 
      scale_colour_gradient(low = "goldenrod1", high = "red", limits=c(0,1)) +
      stat_smooth(method=lm,) +
      scale_x_log10() + 
      scale_y_log10() + 
      theme_classic() +
      geom_vline(xintercept = cutoffs$nUMI.min) +
      geom_hline(yintercept = cutoffs$nGene.min) +
      facet_wrap(~capID)
  }
  return(allCaps_QC_joint)
}

## Utility functions for combining and plotting qc and clustering metrics from all caps

merge_metadata <- function(metadata_allCaps) {
  merged_metadata <- list()
  for (slot in names(metadata_allCaps)) {
    merged_metadata[[slot]] <- 
      lapply(metadata_allCaps[[slot]], function(x) tibble(x, cellID = row.names(x))) %>%
      bind_rows() %>% tibble() %>%
      select(capID, cellID, everything())
  }
  return(merged_metadata)
}

create_arrangedPlot <- 
  function(qc_plots, data="filtered_counts", stage="input", plot="pca_cc") {
    arranged_titles <- list(pca_cc = "Cell cycle",
                            pca_mt = "Mito ratio",
                            pca_rb = "Ribo ratio")
    
    if (grepl("pca_*", plot)) {
      arrangedPlot <- ggarrange(plotlist=qc_plots[[stage]][[plot]], legend="none") %>%
        annotate_figure(top=text_grob(arranged_titles[[plot]],
                                      face="bold", size=20))
      return(arrangedPlot)
    } else if (plot=="qc_ridges") {
      arrangedPlot <- plotQC_ridges(merged_metadata[[stage]], cutoffs=Abbie10XCutoffs)
      return(arrangedPlot)
    } else if (plot=="qc_joint") {
      arrangedPlot <- plotQC_joint(merged_metadata[[stage]], cutoffs=Abbie10XCutoffs,
                                   color="mitoRatio")
      return(arrangedPlot)
    } else if (plot=="umap") {
      input_umaps_noLeg <- lapply(qc_plots[[stage]][["umap"]], 
                                  function(x) x+theme(legend.position="none"))
      arrangedPlot <- ggarrange(plotlist=input_umaps_noLeg)
      return(arrangedPlot)
    } else if (plot=="elbow") {
      arrangedPlot <- ggarrange(plotlist=qc_plots[[stage]][[plot]], legend="none")
      return(arrangedPlot)
    }
  }

create_arrangedPlots <- 
  function(qc_plots, data="filtered_counts") {
    arrangedPlots <- list()
    for (stage in names(qc_plots)) {
      for (plot in names(qc_plots[[stage]])) {
        arrangedPlots[[data]][[stage]][[plot]] <- create_arrangedPlot(qc_plots, stage=stage, plot=plot)
      }
    }
    return(arrangedPlots)
  }

save_arrangedPlot <- 
  function(arrangedPlots, data="filtered_counts", stage="input", plot="pca_cc") {
    path <- here("plots", data, stage, paste0(plot, ".png"))
    print(paste0("Saving arrangedPlots$",data,"$",stage,"$",plot," ",
                 "here: ", path))
    ggsave(path,
           arrangedPlots[[data]][[stage]][[plot]],
           height = 8, width = 8, units = "in", bg="white")
  }

save_arrangedPlots <- 
  function(arrangedPlots, data="filtered_counts") {
    for (stage in names(arrangedPlots[[data]])) {
      for (plot in names(arrangedPlots[[data]][[stage]])) {
        if (!dir.exists(here("plots", data, stage))) {
          dir.create(here("plots", data, stage, plot))
        }
        save_arrangedPlot(arrangedPlots, data=data, stage=stage, plot=plot)
      }
    }
  }