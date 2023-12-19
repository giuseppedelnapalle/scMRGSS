#!/usr/bin/env Rscript

library(Seurat)
library(ggplot2)

# plot average expression of genes from two groups
scatter_plot_genes_avg_exp <- function(object, assay="RNA", group_by="orig.ident",
                               title="", path_file){
  avg <- as.data.frame(log1p(AverageExpression(object, assays=assay, 
                                               group.by=group_by, verbose = FALSE)[[assay]]))
  p <- ggplot(avg, aes_string(colnames(avg)[1], colnames(avg)[2])) + geom_point() + ggtitle(title)
  ggsave(path_file, plot = p, width = 18, height = 12, units = "cm")
}

log1p_avg <- function(mat){
  log1p(apply(mat, 1, mean))
}
