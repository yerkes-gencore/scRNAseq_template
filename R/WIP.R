# Trajectory inference

```{r}
library(slingshot)
sce <- as.SingleCellExperiment(combined.obj)
sce <- slingshot(sce, clusterLabels = 'seurat_clusters', reducedDim = 'UMAP')
test <- as.Seurat(sce)

dimred <- combined.obj@reductions$umap
clustering <- combined.obj$seurat_clusters
counts <- as.matrix(combined.obj@assays$RNA@counts[combined.obj@assays$RNA@var.features,])
lineages <- getLineages(data = dimred, clusterLabels = clustering)
```


```{r}
DimPlot(combined.obj, reduction = 'umap')
ssds <- SlingshotDataSet(sce)
pto <- as.PseudotimeOrdering(sce)
lines(SlingshotDataSet(sce))
```


## Most common features

Plotting highly expressed genes. We expect to see constitutively expressed genes (e.g. mitochondrial, ribosomal, etc).

```{r most_abundant_genes, fig.width=10, fig.height=4}
plot_top_genes <- function(data, n=20){
  C <- data@assays$RNA@counts
  C <- Matrix::t(Matrix::t(C)/Matrix::colSums(C)) * 100
  most_expressed <- order(apply(C, 1, median), decreasing = T)[20:1]
  par(mar=c(5,12,4,1)+.1)
  boxplot(as.matrix(t(C[most_expressed, ])), las=1, cex=0.1, xlab = "% total count per cell", col = (scales::hue_pal())(20)[20:1], horizontal = TRUE, main=data@project.name) 
}
tmp <- lapply(captures, plot_top_genes)
```


##How much of each capture is in each cluster

```{r}
if (hashed){
  tmp <- as.data.frame.matrix(table(combined.obj$seurat_clusters,
                                    combined.obj$ID))
  knitr::kable(tmp, caption = 'Number of cells')
  tmp <- tmp %>% mutate(across(.cols= everything(), ~ 100*.x / sum(.x))) #, .names = "{.col}_percent"
  knitr::kable(tmp, caption = 'Percentage of cells from capture', digits = 1)
} else {
  tmp <- as.data.frame.matrix(table(combined.obj$seurat_clusters,
                                    combined.obj$orig.ident))
  knitr::kable(tmp, caption = 'Number of cells')
  tmp <- tmp %>% mutate(across(.cols= everything(), ~ 100*.x / sum(.x))) #, .names = "{.col}_percent"
  knitr::kable(tmp, caption = 'Percentage of cells from capture', digits = 1)
}

```