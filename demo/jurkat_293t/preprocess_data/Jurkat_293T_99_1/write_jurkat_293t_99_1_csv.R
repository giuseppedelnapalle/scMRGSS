#!/usr/bin/env Rscript
# write Jurkat_293T_50_50 dataset to CSV files

options(stringsAsFactors = F)
library(Seurat)
library(ggplot2)

# working directory
wd <- "/home/nikola/Project_Data/R_data/tests/Jurkat_293T/Jurkat_293T_99_1"
setwd(wd)

# dataset names
ds_n <- "jurkat_293t_99_1"

# data directories
d_rt <- "/home/nikola/Documents/Rsch/resources/scRNA_seq/Jurkat_293T"
dir_dt <- paste(d_rt, ds_n, sep = "/")
dir_analysis <- "/home/nikola/Documents/Rsch/resources/scRNA_seq/Jurkat_293T/jurkat_293t_99_1_analysis/analysis_csv"
dir_tsne <- paste(dir_analysis, "tsne", sep = "/")
dir_pca <- paste(dir_analysis, "pca", sep = "/")

# output directory
dir_o <- paste(wd, "output", sep = "/")
# dir.create(dir_o, recursive = T)

# 1 load data -------------------------------------------------------------

# read 10x files
mat <- Read10X(data.dir = dir_dt)

head(mat)[,1:6]
dim(mat)
# [1] 32738  4185

t_sne <- read.csv(paste(dir_tsne, "projection.csv", sep = "/"), header = T)
head(t_sne)
dim(t_sne)

pca <- read.csv(paste(dir_pca, "projection.csv", sep = "/"), header = T)
head(pca)
dim(pca)

# 2 annotate cells --------------------------------------------------------

# initialize Seurat object
seurat <- CreateSeuratObject(counts = mat, project = ds_n)

# store mitochondrial percentage in object meta data
seurat <- PercentageFeatureSet(seurat, pattern = "^MT-", col.name = "percent.mt")

# run sctransform
# remove confounding sources of variation
seurat <- SCTransform(seurat, vars.to.regress = "percent.mt", verbose = FALSE)

# These are now standard steps in the Seurat workflow for visualization and clustering
seurat <- RunPCA(seurat, verbose = FALSE)
seurat <- RunUMAP(seurat, dims = 1:30, verbose = FALSE)

seurat <- FindNeighbors(seurat, dims = 1:30, verbose = FALSE)
seurat <- FindClusters(seurat, resolution = .2, verbose = FALSE)

DimPlot(seurat, reduction = "umap", label = TRUE) + NoLegend()
# DimPlot(seurat, reduction = "pca", label = TRUE) + NoLegend()
fig <- DimPlot(seurat, reduction = "umap", label = TRUE) + NoLegend()
ggsave(paste(dir_o, paste0("umap_", ds_n, ".pdf"), sep = "/"),
       width = 16, height = 16,units = "cm")

# clustering result
table(seurat$seurat_clusters)
# 0    1    2    3    4    5 
# 1783 1148  597  355  272   30
head(seurat$seurat_clusters)

# create meta data
c_ann <- ifelse(seurat$seurat_clusters != 5, "jurkat", "293t")
head(c_ann)

# 3 examine embedding results ---------------------------------------------

# check if barcodes of t_sne & mat are in the same order
sum(match(t_sne$Barcode, colnames(mat)) == c(1:nrow(t_sne))) == nrow(t_sne)

df <- cbind(t_sne, data.frame(label=c_ann))
head(df)

# scatterplots
# t-SNE
fig <- ggplot(df, aes(x=TSNE.1, y=TSNE.2, color=label)) + geom_point(size=2)
ggsave(paste(dir_o, paste0("scatterplot_tsne_", ds_n, ".pdf"), sep = "/"), plot = fig, 
       width = 12, height = 8, units = "in", dpi = 300)

# PCA
sum(match(pca$Barcode, colnames(mat)) == c(1:nrow(pca))) == nrow(pca)
df_2 <- cbind(pca, data.frame(label=c_ann))
fig <- ggplot(df_2, aes(x=PC.1, y=PC.2, color=label)) + geom_point(size=2)
ggsave(paste(dir_o, paste0("scatterplot_pca_", ds_n, ".pdf"), sep = "/"), plot = fig, 
       width = 12, height = 8, units = "in", dpi = 300)

# 4 write data to CSV files -----------------------------------------------

write.csv(mat, file = paste(dir_o, "mat_jurkat_293t_99_1.csv", sep = "/"),
          quote = F, row.names = T)

# meta data
c_ann <- ifelse(seurat$seurat_clusters != 5, "jurkat", "293t")
head(c_ann)

saveRDS(c_ann, file = paste(dir_o, "c_ann_jurkat_293t_99_1.rds", sep = "/"))
write.table(c_ann, file = paste(dir_o, "c_ann_jurkat_293t_99_1.csv", sep = "/"),
          quote = F, sep = ",", row.names = F, col.names = F)
