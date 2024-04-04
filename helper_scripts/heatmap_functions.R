# *** Avoid issues with mismatched SubjectID and column params ***
## Problem: I think when DGEList is called, the columns (SubjectIDs) are put in alphabetical order
## Solution (safest for whole project): When creating bulk DGEList, modify column order to match the metadata order desired for all downstream plotting.
## Solution (most expedient for plots): Input all column metadata as a dataframe and re-order final column for `hmap` and column metadata tbl.
## Solution: Do both of the above and check the final `hmap` object before plotting heatmap.

# *** Avoid issues with mismatched genes and row params ***
## Input gene list as a vector if no row-splitting or row label coloring is desired
## Otherwise, input gene list as tibble `gene_tbl` with a `gene` column.
##  If coloring row labels, include a `gene_label_color` column with the desired colors for each gene.
##  If splitting rows, include a `row_split_label` column with the desired row splitting labels.

# *** Avoid issues with one heatmap wrapper interfering with another ***
## Re-write a base ComplexHeatmap::Heatmap wrapper to handle gene list or dataframe inputs, column/row splitting and coloring

## *** Heatmap ***

plotHeatmap <- function(exprs_mat = NULL, gene_tbl = NULL, sample_tbl = NULL,
                        sampleID = "sample",
                        cluster_genes = FALSE, cluster_samples = FALSE,
                        box_width = unit(3.5, "mm"), width_buffer = unit(5, "mm"),
                        box_height = unit(3.5, "mm"), height_buffer = unit(10, "mm"),
                        column_title = " ", column_gap = unit(2, "mm"),
                        colors = c("blue", "white", "red"),
                        scale_min = -2, scale_max = 2,
                        heatmap_legend_param = list(at = c(scale_min, 0, scale_max),
                                                    labels = c(scale_min, 0, scale_max),
                                                    title = "log2 fold\ndifference\nfrom\nmedian\nexpression"),
                        slice_labels_rot = 90) {
  ## Check for missing/misordered sample names
  # If the colnames(exprs_mat) do not match sample_tbl, warn and rearrange columns to match sample_tbl
  if (!setequal(sample_tbl[[sampleID]], colnames(exprs_mat))) {
    message(paste0("The samples and/or order in `sample_tbl` do not match `colnames(exprs_mat)`. Subsetting `exprs_mat` to match `sample_tbl`."))
    missing_from_exprs_mat <- sample_tbl[[sampleID]][!sample_tbl[[sampleID]] %in% colnames(exprs_mat)]
    if (length(missing_from_exprs_mat) > 0) {
      message(paste0("The following genes in `sample_tbl[[", sampleID, "]]` are missing in `colnames(exprs_mat)`: ", paste0(missing_from_exprs_mat, collapse = ", ")))
    }
    missing_from_sample_tbl <- colnames(exprs_mat)[!colnames(exprs_mat) %in% sample_tbl[[sampleID]]]
    if (length(missing_from_sample_tbl) > 0) {
      message(paste0("The following genes in `colnames(exprs_mat)` are missing in `sample_tbl[[", sampleID, "]]`: ", paste0(missing_from_sample_tbl, collapse = ", ")))
    }
    exprs_mat <- exprs_mat[,sample_tbl[[sampleID]]]
  }

  ## Check for missing genes
  # Genes in gene_tbl missing in rownames(exprs_mat)
  missing_genes <- gene_tbl$gene[!gene_tbl$gene %in% rownames(exprs_mat)]
  if (length(missing_genes) > 0) {
    # If there are any missing_genes, filter them out of the gene_tbl and report them in a message.
    gene_tbl <- gene_tbl[gene_tbl$gene %in% rownames(exprs_mat), ]
    if (length(gene_tbl$gene) == 0) {
      stop("No data for requested genes")
    }
    if (length(gene_tbl$gene) > 0) {
      message(paste0("Genes ", paste0(missing_genes, collapse = ", "),
                     " not found in `exprs_mat`."))
    }
  }

  # Add color column if missing
  if (is.null(suppressWarnings(sample_tbl$color))) {
    sample_tbl$color <- "black"
  }

  ## Create final matrix
  hm_mat <- exprs_mat[gene_tbl$gene, ]
  baseline <- matrixStats::rowMedians(hm_mat)
  hm_mat <- hm_mat - baseline

  ## Plot heatmap
  ComplexHeatmap::Heatmap(
    hm_mat,
    heatmap_legend_param = heatmap_legend_param,
    width = ncol(hm_mat) * box_width + width_buffer,
    height = nrow(hm_mat) * box_height + height_buffer,
    column_title = column_title,
    cluster_rows = cluster_genes,
    cluster_columns = cluster_samples,
    column_split = if (!is.null(suppressWarnings(sample_tbl$group))) {
      sample_tbl$group
    } else {
      NULL
    },
    column_names_gp = grid::gpar(col = sample_tbl$color),
    row_split = if (!is.null(suppressWarnings(gene_tbl$split_label))) {
      gene_tbl$split_label
    } else {
      NULL
    },
    row_names_gp = grid::gpar(col = gene_tbl$color),
    top_annotation = (if (!is.null(suppressWarnings(sample_tbl$group))) {
      ComplexHeatmap::HeatmapAnnotation(
        foo = ComplexHeatmap::anno_block(
          gp = grid::gpar(col = NA),
          labels = distinct(sample_tbl, group, color)$group,
          labels_gp = grid::gpar(col = distinct(sample_tbl, group, color)$color, fontsize = 10),
          labels_rot = slice_labels_rot,
          height = unit(2, "cm")))
    }
    else {
      NULL
    }),
    column_gap = column_gap,
    col = circlize::colorRamp2(c(scale_min, 0, scale_max), colors)
  )
}

## Function to create top genes table

getTopDegTbl <- function(
    resultsTable, contrast_id,
    groupID = "grp_chl", sampleID = "sample", arrange_by = "pvalue", direction = "unequal",
    padj_cutoff = 0.05, slice_n = 50, filter_pattern = NULL) {

  if (! arrange_by %in% c("pvalue", "log2FoldChange")) {
    stop("'arrange_by' parameter must be set to 'pvalue' (default) or 'log2FoldChange'")
  }

  if (arrange_by == "log2FoldChange" | !is.null(padj_cutoff)) {
    if (nrow(filter(resultsTable,
                    contrast == contrast_id,
                    padj < padj_cutoff)) == 0) {
      print("No significant DEGs.")
      return(NULL)
    } else if (nrow(filter(resultsTable,
                           contrast == contrast_id,
                           padj < padj_cutoff)) == 1) {
      print("Only one significant DEG.")
      return(NULL)
    }
  }

  deg_tbl <- resultsTable %>% filter(contrast == contrast_id)
  if (!is.null(filter_pattern)) {
    deg_tbl <- deg_tbl %>% filter(!grepl(filter_pattern, gene))
  }
  if (direction == "equal") {
    deg_tbl <- deg_tbl %>%
      group_by(-sign(log2FoldChange))
    if (arrange_by == "pvalue") {
      deg_tbl <- deg_tbl %>% arrange(pvalue) %>% slice_head(n = round(slice_n/2)) %>%
        ungroup()
    } else if (arrange_by == "log2FoldChange") {
      deg_tbl <- deg_tbl %>% filter(padj < padj_cutoff) %>%
        arrange(desc(abs(log2FoldChange))) %>%
        slice_head(n = round(slice_n/2)) %>%
        ungroup() %>%
        arrange(desc(log2FoldChange))
    }
  } else if (direction == "up") {
    deg_tbl <- deg_tbl %>% filter(log2FoldChange > 0)
    if (arrange_by == "pvalue") {
      deg_tbl <- deg_tbl %>% arrange(pvalue) %>%
        slice_head(n = slice_n)
    } else if (arrange_by == "log2FoldChange") {
      deg_tbl <- deg_tbl %>% arrange(desc(log2FoldChange)) %>%
        filter(padj < padj_cutoff) %>% slice_head(n = slice_n)
    }
  } else if (direction == "down") {
    deg_tbl <- deg_tbl %>% filter(log2FoldChange < 0)
    if (arrange_by == "pvalue") {
      deg_tbl <- deg_tbl %>% arrange(pvalue) %>% slice_head(n = slice_n)
    } else if (arrange_by == "log2FoldChange") {
      deg_tbl <- deg_tbl %>% arrange(log2FoldChange) %>%
        filter(padj < padj_cutoff) %>% slice_head(n = slice_n)
    }
  } else if (direction == "unequal") {
    if (arrange_by == "pvalue") {
      deg_tbl <- deg_tbl %>% arrange(pvalue) %>% slice_head(n = slice_n)
    } else if (arrange_by == "log2FoldChange") {
      deg_tbl <- deg_tbl %>% arrange(desc(abs(log2FoldChange))) %>%
        filter(padj < padj_cutoff) %>% slice_head(n = slice_n) %>%
        arrange(desc(log2FoldChange))
    }
  } else {
    stop("'direction' parameter must be set to 'uneven' (default), 'equal', 'up', or 'down'")
  }
  return(deg_tbl)
}
