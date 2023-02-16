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