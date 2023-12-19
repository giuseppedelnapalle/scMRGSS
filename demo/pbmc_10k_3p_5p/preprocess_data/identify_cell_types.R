#!/usr/bin/env Rscript
# assign cell type identities to clusters
# reference https://satijalab.org/seurat/articles/pbmc3k_tutorial.html

options(stringsAsFactors = F)
library(Seurat)
library(SingleR)
# library(celldex)
library(ggplot2)
library(dplyr)

# working directory
wd <- "/home/nikola/Project_Data/R_data/pbmc_10k_integration"
setwd(wd)

# input directory
# dir_i <- paste(wd, "input", sep = "/")
dir_i <- "/home/nikola/Documents/Rsch/resources/scRNA_seq/blood/PBMC_10k"
# dir.create(dir_i, recursive = T)

# output directory
dir_o <- paste(wd, "output", sep = "/")
# dir.create(dir_o, recursive = T)

# load functions
# source("custom_seurat_functions.R")
# source("helper_functions.R")

# dataset names
ds_n <- "pbmc10k_3p"
ds_n_2 <- "pbmc10k_5p"

# 1 dataset 1 -------------------------------------------------------------

# set up Seurat object
matrix_3p <- Read10X_h5(paste(dir_i, "3p_pbmc10k_filt.h5", sep = "/"),use.names = T)
srat_3p   <- CreateSeuratObject(matrix_3p,project = "pbmc10k_3p")

rm(matrix_3p)

# QC and selecting cells for further analysis
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
srat_3p[["percent.mt"]] <- PercentageFeatureSet(srat_3p, pattern = "^MT-")
srat_3p[["percent.rb"]] <- PercentageFeatureSet(srat_3p, pattern = "^RP[SL]")

# Visualize QC metrics as a violin plot
VlnPlot(srat_3p, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(srat_3p, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(srat_3p, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# plot1 + plot2

srat_3p <- subset(srat_3p, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 12)

# Normalizing the data

# global-scaling normalization --------------------------------------------

# global-scaling normalization
# srat_3p <- NormalizeData(srat_3p, normalization.method = "LogNormalize", scale.factor = 10000)

# # Identification of highly variable features (feature selection)
# srat_3p <- FindVariableFeatures(srat_3p, selection.method = "vst", nfeatures = 2000)
# 
# # Identify the 10 most highly variable genes
# top10 <- head(VariableFeatures(srat_3p), 10)
# 
# # plot variable features with and without labels
# plot1 <- VariableFeaturePlot(srat_3p)
# plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
# plot1 + plot2
# 
# # Scaling the data
# all.genes <- rownames(srat_3p)
# system.time({
#   srat_3p <- ScaleData(srat_3p, features = all.genes)
# })

# SCTransform -------------------------------------------------------------

# Single SCTransform command replaces NormalizeData(), ScaleData(),
# and FindVariableFeatures()
srat_3p <- SCTransform(srat_3p, method = "glmGamPoi", 
                       vars.to.regress = "percent.mt", verbose = FALSE)

# Perform linear dimensional reduction ------------------------------------

# Perform linear dimensional reduction
srat_3p <- RunPCA(srat_3p, features = VariableFeatures(object = srat_3p))

# Examine and visualize PCA results a few different ways
print(srat_3p[["pca"]], dims = 1:5, nfeatures = 5)

DimPlot(srat_3p, reduction = "pca")

DimHeatmap(srat_3p, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(srat_3p, dims = 1:15, cells = 500, balanced = TRUE)

# Determine the ‘dimensionality’ of the dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time

# srat_3p <- JackStraw(srat_3p, num.replicate = 100)
# # Error in JackStraw(srat_3p, num.replicate = 100) : 
# #   JackStraw cannot be run on SCTransform-normalized data.
# # Please supply a non-SCT assay.
# 
# srat_3p <- ScoreJackStraw(srat_3p, dims = 1:20)
# 
# JackStrawPlot(srat_3p, dims = 1:15)

ElbowPlot(srat_3p)

# Cluster the cells
srat_3p <- FindNeighbors(srat_3p, dims = 1:10)
srat_3p <- FindClusters(srat_3p, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(srat_3p), 5)

# Run non-linear dimensional reduction (UMAP/tSNE)
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
srat_3p <- RunUMAP(srat_3p, dims = 1:10)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
# DimPlot(srat_3p, reduction = "umap", label = T)
fig <- DimPlot(srat_3p, reduction = "umap", label = T)
# ggsave(paste(dir_o, paste0(paste("dim_red_plot",ds_n,"SCTransform",sep = "_"),".png"),sep = "/"), 
#        plot = fig, width = 18, height = 12, units = "cm")
ggsave(paste(dir_o, paste0(paste("dim_red_plot",ds_n,"SCTransform",sep = "_"),".pdf"),sep = "/"), 
       plot = fig, width = 18, height = 12, units = "cm")

# Finding differentially expressed features (cluster biomarkers)
# # find all markers of cluster 2
# cluster2.markers <- FindMarkers(srat_3p, ident.1 = 2, min.pct = 0.25)
# head(cluster2.markers, n = 5)
# 
# # find all markers distinguishing cluster 5 from clusters 0 and 3
# cluster5.markers <- FindMarkers(srat_3p, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
# head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
# It is recommended to do differential expression on RNA assay, and not the SCTransform
DefaultAssay(srat_3p) <- "RNA"
srat_3p <- NormalizeData(srat_3p)
srat_3p <- FindVariableFeatures(srat_3p, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(srat_3p)
srat_3p <- ScaleData(srat_3p, features = all.genes)

system.time({
  srat_3p.markers <- FindAllMarkers(srat_3p, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
})
srat_3p.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

dim(srat_3p.markers)
head(srat_3p.markers)

# the ROC test returns the ‘classification power’ for any individual marker 
# (ranging from 0 - random, to 1 - perfect)
# cluster0.markers <- FindMarkers(srat_3p, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

VlnPlot(srat_3p, features = c("MS4A1", "CD79A"))

# you can plot raw counts as well
# VlnPlot(srat_3p, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)

FeaturePlot(srat_3p, features = c("MS4A1", "GNLY", 
                                  "CD3E", # T-cell, NK-cell
                                  "CD14", 
                                  "FCER1A", # Langerhans cell, DC, basophil
                                  "FCGR3A", # mono, NK, MP
                                  "LYZ", "PPBP", "CD8A"))
# Naive CD4+ T
FeaturePlot(srat_3p, features = c("IL7R", "CCR7"))
# CD14+ Mono
FeaturePlot(srat_3p, features = c("CD14", "LYZ"))
# Memory CD4+
FeaturePlot(srat_3p, features = c("IL7R", "S100A4"))
# CD8+ T
FeaturePlot(srat_3p, features = c("CD8B"))
# B
FeaturePlot(srat_3p, features = c("MS4A1"))
# CD16+ Mono
FeaturePlot(srat_3p, features = c("FCGR3A", "MS4A7"))
# NK
FeaturePlot(srat_3p, features = c("GNLY", "NKG7"))
# DC
FeaturePlot(srat_3p, features = c("LILRA4", "TPM2"))
# Platelet
FeaturePlot(srat_3p, features = c("PPBP", "GP1BB"))
# MAIT
FeaturePlot(srat_3p, features = c("KLRB1", "CXCR6"))
# Plasmablast
FeaturePlot(srat_3p, features = c("CD38", "CD59"))

srat_3p.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(srat_3p, features = top10$gene) + NoLegend()

# Assigning cell type identity to clusters
# manually assign cell types

# manually assign cell types ----------------------------------------------

n_c <- length(levels(as.factor(srat_3p.markers$cluster)))
srat_3p_mk_top <- sapply(0:(n_c-1), function(i) {
  print(paste0("cluster: ", i))
  print(head(srat_3p.markers[srat_3p.markers$cluster==i,]$gene))
  head(srat_3p.markers[srat_3p.markers$cluster==i,]$gene, n=30)
})

# new.cluster.ids <- c("CD14+ Mono", "Naive CD8+ T", "Helper CD4+ T", "Mono", 
#                      "Naive CD4+ T", "B", "Naive CD4+ T", "CD8+ T", 
#                      "CD16+ Mono","NK","DC", "Platelet",
#                      "B_subtype")
new.cluster.ids <- c("cd14_mono", "naive_cd8_t", "helper_cd4_t", "mono", 
                     "naive_cd4_t", "b", "naive_cd4_t_2", "effector_t", 
                     "cd16_mono","nk","dc", "plt",
                     "b_sub")

names(new.cluster.ids) <- levels(srat_3p)
srat_3p <- RenameIdents(srat_3p, new.cluster.ids)
fig <- DimPlot(srat_3p, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
# ggsave(paste(dir_o, paste0(paste("dim_red_plot",ds_n,"curated",sep = "_"),".png"),sep = "/"), 
#        plot = fig, width = 18, height = 12, units = "cm")
ggsave(paste(dir_o, paste0(paste("dim_red_plot",ds_n,"curated",sep = "_"),".pdf"),sep = "/"), 
       plot = fig, width = 18, height = 12, units = "cm")

# export cell annotation data
c_ann_3p <- srat_3p@active.ident
write.table(c_ann_3p, file = paste(dir_o, paste0("c_ann_", ds_n,".csv"), sep = "/"), sep = ",", 
            quote = FALSE, col.names = FALSE)
saveRDS(c_ann_3p, file = paste(dir_o, paste0("c_ann_", ds_n,".rds"), sep = "/"))

# Cell type annotation using SingleR --------------------------------------

# get reference datasets from celldex package
hpca.ref <- celldex::HumanPrimaryCellAtlasData()
dice.ref <- celldex::DatabaseImmuneCellExpressionData()
monaco.ref <- celldex::MonacoImmuneData()

# convert our Seurat object to single cell experiment (SCE)
sce_3p <- as.SingleCellExperiment(srat_3p)
sce_3p

# hpca
hpca.main <- SingleR(test = sce_3p,assay.type.test = 1,ref = hpca.ref,labels = hpca.ref$label.main)
hpca.fine <- SingleR(test = sce_3p,assay.type.test = 1,ref = hpca.ref,labels = hpca.ref$label.fine)
# dice
dice.main <- SingleR(test = sce_3p,assay.type.test = 1,ref = dice.ref,labels = dice.ref$label.main)
dice.fine <- SingleR(test = sce_3p,assay.type.test = 1,ref = dice.ref,labels = dice.ref$label.fine)
# monaco
monaco.main <- SingleR(test = sce_3p,assay.type.test = 1,ref = monaco.ref,labels = monaco.ref$label.main)
monaco.fine <- SingleR(test = sce_3p,assay.type.test = 1,ref = monaco.ref,labels = monaco.ref$label.fine)

table(hpca.main$pruned.labels)
table(dice.main$pruned.labels)
table(monaco.main$pruned.labels)

table(hpca.fine$pruned.labels)
table(dice.fine$pruned.labels)
table(monaco.fine$pruned.labels)

# add the annotations to the Seurat object metadata
srat_3p@meta.data$hpca.main   <- hpca.main$pruned.labels
srat_3p@meta.data$dice.main   <- dice.main$pruned.labels
srat_3p@meta.data$monaco.main <- monaco.main$pruned.labels
srat_3p@meta.data$hpca.fine   <- hpca.fine$pruned.labels
srat_3p@meta.data$dice.fine   <- dice.fine$pruned.labels
srat_3p@meta.data$monaco.fine <- monaco.fine$pruned.labels

# visualize the fine-grained annotations
srat_3p <- SetIdent(srat_3p, value = "hpca.fine")
fig <- DimPlot(srat_3p, label = T , repel = T, label.size = 3) + NoLegend()
# ggsave(paste(dir_o, paste0(paste("dim_red_plot",ds_n,"hpca.fine",sep = "_"),".png"),sep = "/"), 
#        plot = fig, width = 18, height = 12, units = "cm")
ggsave(paste(dir_o, paste0(paste("dim_red_plot",ds_n,"hpca.fine",sep = "_"),".pdf"),sep = "/"), 
       plot = fig, width = 18, height = 12, units = "cm")

srat_3p <- SetIdent(srat_3p, value = "dice.fine")
fig <- DimPlot(srat_3p, label = T , repel = T, label.size = 3) + NoLegend()
# ggsave(paste(dir_o, paste0(paste("dim_red_plot",ds_n,"dice.fine",sep = "_"),".png"),sep = "/"), 
#        plot = fig, width = 18, height = 12, units = "cm")
ggsave(paste(dir_o, paste0(paste("dim_red_plot",ds_n,"dice.fine",sep = "_"),".pdf"),sep = "/"), 
       plot = fig, width = 18, height = 12, units = "cm")

srat_3p <- SetIdent(srat_3p, value = "monaco.fine")
fig <- DimPlot(srat_3p, label = T , repel = T, label.size = 3) + NoLegend()
# ggsave(paste(dir_o, paste0(paste("dim_red_plot",ds_n,"monaco.fine",sep = "_"),".png"),sep = "/"), 
#        plot = fig, width = 18, height = 12, units = "cm")
ggsave(paste(dir_o, paste0(paste("dim_red_plot",ds_n,"monaco.fine",sep = "_"),".pdf"),sep = "/"), 
       plot = fig, width = 18, height = 12, units = "cm")

rm(plot1, plot2)

# saveRDS(srat_3p, file = paste(dir_o, "srat_3p.rds", sep = "/"))
saveRDS(srat_3p.markers, file = paste(dir_o, "srat_3p.markers.rds", sep = "/"))

# load data
# load("identify_cell_types.RData")
# srat_3p <- readRDS(paste(dir_o, "srat_3p.rds", sep = "/"))
# srat_3p.markers <- readRDS(paste(dir_o, "srat_3p.markers.rds", sep = "/"))

# 2 dataset 2 -------------------------------------------------------------

# set up Seurat object
matrix_5p <- Read10X_h5(paste(dir_i, "5p_pbmc10k_filt.h5", sep = "/"),use.names = T)$`Gene Expression`
srat_5p   <- CreateSeuratObject(matrix_5p,project = "pbmc10k_5p")

rm(matrix_5p)

# QC and selecting cells for further analysis
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
srat_5p[["percent.mt"]] <- PercentageFeatureSet(srat_5p, pattern = "^MT-")
srat_5p[["percent.rb"]] <- PercentageFeatureSet(srat_5p, pattern = "^RP[SL]")

# Visualize QC metrics as a violin plot
VlnPlot(srat_5p, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(srat_5p, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(srat_5p, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# plot1 + plot2

srat_5p <- subset(srat_5p, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 8)

# Normalizing the data

# global-scaling normalization --------------------------------------------

# # global-scaling normalization
# srat_5p <- NormalizeData(srat_5p, normalization.method = "LogNormalize", scale.factor = 10000)
# 
# # Identification of highly variable features (feature selection)
# srat_5p <- FindVariableFeatures(srat_5p, selection.method = "vst", nfeatures = 2000)
# 
# # Identify the 10 most highly variable genes
# top10 <- head(VariableFeatures(srat_5p), 10)
# 
# # plot variable features with and without labels
# plot1 <- VariableFeaturePlot(srat_5p)
# plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
# # plot1 + plot2
# 
# # Scaling the data
# all.genes <- rownames(srat_5p)
# system.time({
#   srat_5p <- ScaleData(srat_5p, features = all.genes)
# })

# SCTransform -------------------------------------------------------------

# Single SCTransform command replaces NormalizeData(), ScaleData(),
# and FindVariableFeatures()
srat_5p <- SCTransform(srat_5p, method = "glmGamPoi",
                       vars.to.regress = "percent.mt", verbose = FALSE)

# Perform linear dimensional reduction ------------------------------------

# Perform linear dimensional reduction
srat_5p <- RunPCA(srat_5p, features = VariableFeatures(object = srat_5p))

# Examine and visualize PCA results a few different ways
print(srat_5p[["pca"]], dims = 1:5, nfeatures = 5)

DimPlot(srat_5p, reduction = "pca")

DimHeatmap(srat_5p, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(srat_5p, dims = 1:15, cells = 500, balanced = TRUE)

# Determine the ‘dimensionality’ of the dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time

# srat_5p <- JackStraw(srat_5p, num.replicate = 100)
# # Error in JackStraw(srat_5p, num.replicate = 100) : 
# #   JackStraw cannot be run on SCTransform-normalized data.
# # Please supply a non-SCT assay.
# 
# srat_5p <- ScoreJackStraw(srat_5p, dims = 1:20)
# 
# JackStrawPlot(srat_5p, dims = 1:15)

ElbowPlot(srat_5p)

# Cluster the cells
srat_5p <- FindNeighbors(srat_5p, dims = 1:10)
srat_5p <- FindClusters(srat_5p, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(srat_5p), 5)

# Run non-linear dimensional reduction (UMAP/tSNE)
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
srat_5p <- RunUMAP(srat_5p, dims = 1:10)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
fig <- DimPlot(srat_5p, reduction = "umap", label = T)
# ggsave(paste(dir_o, paste0(paste("dim_red_plot",ds_n_2,"SCTransform",sep = "_"),".png"),sep = "/"),
#        plot = fig, width = 18, height = 12, units = "cm")
ggsave(paste(dir_o, paste0(paste("dim_red_plot",ds_n_2,"SCTransform",sep = "_"),".pdf"),sep = "/"),
       plot = fig, width = 18, height = 12, units = "cm")
# ggsave(paste(dir_o, paste0(paste("dim_red_plot",ds_n_2,sep = "_"),".png"),sep = "/"), 
#        plot = fig, width = 18, height = 12, units = "cm")

# Finding differentially expressed features (cluster biomarkers)
# # find all markers of cluster 2
# cluster2.markers <- FindMarkers(srat_5p, ident.1 = 2, min.pct = 0.25)
# head(cluster2.markers, n = 5)
# 
# # find all markers distinguishing cluster 5 from clusters 0 and 3
# cluster5.markers <- FindMarkers(srat_5p, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
# head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
# It is recommended to do differential expression on RNA assay, and not the SCTransform
DefaultAssay(srat_5p) <- "RNA"
srat_5p <- NormalizeData(srat_5p)
srat_5p <- FindVariableFeatures(srat_5p, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(srat_5p)
srat_5p <- ScaleData(srat_5p, features = all.genes)

srat_5p.markers <- FindAllMarkers(srat_5p, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
srat_5p.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

dim(srat_5p.markers)
head(srat_5p.markers)

# the ROC test returns the ‘classification power’ for any individual marker 
# (ranging from 0 - random, to 1 - perfect)
# cluster0.markers <- FindMarkers(srat_5p, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

VlnPlot(srat_5p, features = c("MS4A1", "CD79A"))

# you can plot raw counts as well
# VlnPlot(srat_5p, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)

FeaturePlot(srat_5p, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", 
                                  "FCGR3A", "LYZ", "PPBP", "CD8A"))
# Naive CD4+ T
FeaturePlot(srat_5p, features = c("IL7R", "CCR7"))
# CD14+ Mono
FeaturePlot(srat_5p, features = c("CD14", "LYZ"))
# Memory CD4+
FeaturePlot(srat_5p, features = c("IL7R", "S100A4"))
# CD8+ T
FeaturePlot(srat_5p, features = c("CD8B"))
# B
FeaturePlot(srat_5p, features = c("MS4A1"))
# CD16+ Mono
FeaturePlot(srat_5p, features = c("FCGR3A", "MS4A7"))
# NK
FeaturePlot(srat_5p, features = c("GNLY", "NKG7"))
# DC
FeaturePlot(srat_5p, features = c("LILRA4", "TPM2"))
# Platelet
FeaturePlot(srat_5p, features = c("PPBP", "GP1BB"))
# MAIT
FeaturePlot(srat_5p, features = c("KLRB1", "CXCR6"))
# Plasmablast
FeaturePlot(srat_5p, features = c("CD38", "CD59"))

srat_5p.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(srat_5p, features = top10$gene) + NoLegend()

# Assigning cell type identity to clusters
# manually assign cell types

# manually assign cell types ----------------------------------------------

n_c <- length(levels(as.factor(srat_5p.markers$cluster)))
srat_5p_mk_top <- sapply(0:(n_c-1), function(i) {
  print(paste0("cluster: ", i))
  print(head(srat_5p.markers[srat_5p.markers$cluster==i,]$gene))
  head(srat_5p.markers[srat_5p.markers$cluster==i,]$gene, n=30)
})

new.cluster.ids <- c("cd14_mono", "naive_cd4_t", "helper_cd4_t", "naive_cd4_t_2", 
                     "naive_cd8_t", "b", "mono", "effector_t", 
                     "cd14_mono","nk","mdc", "cd16_mono",
                     "dc", "plt")

names(new.cluster.ids) <- levels(srat_5p)
srat_5p <- RenameIdents(srat_5p, new.cluster.ids)
fig <- DimPlot(srat_5p, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
# ggsave(paste(dir_o, paste0(paste("dim_red_plot",ds_n_2,"curated",sep = "_"),".png"),sep = "/"), 
#        plot = fig, width = 18, height = 12, units = "cm")
ggsave(paste(dir_o, paste0(paste("dim_red_plot",ds_n_2,"curated",sep = "_"),".pdf"),sep = "/"), 
       plot = fig, width = 18, height = 12, units = "cm")

# export cell annotation data
c_ann_5p <- srat_5p@active.ident
write.table(c_ann_5p, file = paste(dir_o, paste0("c_ann_", ds_n_2,".csv"), sep = "/"), sep = ",", 
            quote = FALSE, col.names = FALSE)
saveRDS(c_ann_5p, file = paste(dir_o, paste0("c_ann_", ds_n_2,".rds"), sep = "/"))

# Cell type annotation using SingleR --------------------------------------

# get reference datasets from celldex package
hpca.ref <- celldex::HumanPrimaryCellAtlasData()
dice.ref <- celldex::DatabaseImmuneCellExpressionData()
monaco.ref <- celldex::MonacoImmuneData()

# convert our Seurat object to single cell experiment (SCE)
sce_5p <- as.SingleCellExperiment(srat_5p)
sce_5p

# hpca
system.time({
  hpca.main <- SingleR(test = sce_5p,assay.type.test = 1,ref = hpca.ref,labels = hpca.ref$label.main)
  hpca.fine <- SingleR(test = sce_5p,assay.type.test = 1,ref = hpca.ref,labels = hpca.ref$label.fine)
})
# dice
system.time({
  dice.main <- SingleR(test = sce_5p,assay.type.test = 1,ref = dice.ref,labels = dice.ref$label.main)
  dice.fine <- SingleR(test = sce_5p,assay.type.test = 1,ref = dice.ref,labels = dice.ref$label.fine)
})
# monaco
system.time({
  monaco.main <- SingleR(test = sce_5p,assay.type.test = 1,ref = monaco.ref,labels = monaco.ref$label.main)
  monaco.fine <- SingleR(test = sce_5p,assay.type.test = 1,ref = monaco.ref,labels = monaco.ref$label.fine)
})

table(hpca.main$pruned.labels)
table(dice.main$pruned.labels)
table(monaco.main$pruned.labels)

table(hpca.fine$pruned.labels)
table(dice.fine$pruned.labels)
table(monaco.fine$pruned.labels)

# add the annotations to the Seurat object metadata
srat_cp <- srat_5p
srat_cp@meta.data$hpca.main   <- hpca.main$pruned.labels
srat_cp@meta.data$dice.main   <- dice.main$pruned.labels
srat_cp@meta.data$monaco.main <- monaco.main$pruned.labels
srat_cp@meta.data$hpca.fine   <- hpca.fine$pruned.labels
srat_cp@meta.data$dice.fine   <- dice.fine$pruned.labels
srat_cp@meta.data$monaco.fine <- monaco.fine$pruned.labels

# visualize the fine-grained annotations
srat_cp <- SetIdent(srat_cp, value = "hpca.fine")
fig <- DimPlot(srat_cp, label = T , repel = T, label.size = 3) + NoLegend()
# ggsave(paste(dir_o, paste0(paste("dim_red_plot",ds_n_2,"hpca.fine",sep = "_"),".png"),sep = "/"), 
#        plot = fig, width = 18, height = 12, units = "cm")
ggsave(paste(dir_o, paste0(paste("dim_red_plot",ds_n_2,"hpca.fine",sep = "_"),".pdf"),sep = "/"), 
       plot = fig, width = 18, height = 12, units = "cm")

srat_cp <- SetIdent(srat_cp, value = "dice.fine")
fig <- DimPlot(srat_cp, label = T , repel = T, label.size = 3) + NoLegend()
# ggsave(paste(dir_o, paste0(paste("dim_red_plot",ds_n_2,"dice.fine",sep = "_"),".png"),sep = "/"), 
#        plot = fig, width = 18, height = 12, units = "cm")
ggsave(paste(dir_o, paste0(paste("dim_red_plot",ds_n_2,"dice.fine",sep = "_"),".pdf"),sep = "/"), 
       plot = fig, width = 18, height = 12, units = "cm")


srat_cp <- SetIdent(srat_cp, value = "monaco.fine")
fig <- DimPlot(srat_cp, label = T , repel = T, label.size = 3) + NoLegend()
# ggsave(paste(dir_o, paste0(paste("dim_red_plot",ds_n_2,"monaco.fine",sep = "_"),".png"),sep = "/"), 
#        plot = fig, width = 18, height = 12, units = "cm")
ggsave(paste(dir_o, paste0(paste("dim_red_plot",ds_n_2,"monaco.fine",sep = "_"),".pdf"),sep = "/"), 
       plot = fig, width = 18, height = 12, units = "cm")

rm(plot1, plot2)

# saveRDS(srat_5p, file = paste(dir_o, "srat_5p.rds", sep = "/"))
saveRDS(srat_5p.markers, file = paste(dir_o, "srat_5p.markers.rds", sep = "/"))

# srat_5p <- readRDS(paste(dir_o, "srat_5p.rds", sep = "/"))

# save objects and image --------------------------------------------------

save.image(file = "identify_cell_types.RData")
