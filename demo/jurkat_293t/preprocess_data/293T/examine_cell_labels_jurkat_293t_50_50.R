#!/usr/bin/env Rscript
# examine cell labels of jurkat_293t_50_50 dataset

options(stringsAsFactors = F)
library(ggplot2)

# working directory
wd <- "/home/nikola/Project_Data/R_data/tests/Jurkat_293T/293T"
setwd(wd)

# dataset names
ds_n <- "jurkat_293t_50_50"

# data directories
d_rt <- "/home/nikola/Documents/Rsch/resources/scRNA_seq/Jurkat_293T/jurkat 293t_50 50_analysis/analysis_csv"
dir_dt <- paste(d_rt, "tsne", sep = "/")
dir_dt_2 <- paste(d_rt, "pca", sep = "/")

# output directory
dir_o <- paste(wd, "output", "plots", sep = "/")
dir.create(dir_o, recursive = T)

# 1 load data -------------------------------------------------------------

t_sne <- read.csv(paste(dir_dt, "projection.csv", sep = "/"), header = T)
head(t_sne)
dim(t_sne)

pca <- read.csv(paste(dir_dt_2, "projection.csv", sep = "/"), header = T)
head(pca)
dim(pca)

c_ann <- readRDS(paste(dir_o, "c_ann_jurkat_293t_50_50.rds", sep = "/"))
head(c_ann)
dim(c_ann)

# 2 visualisation ---------------------------------------------------------

# check if barcodes of t_sne & c_ann are in the same order
sum(match(t_sne$Barcode, c_ann$barcode) == c(1:nrow(t_sne))) == nrow(t_sne)

df <- cbind(t_sne, data.frame(label=c_ann$label))
head(df)

# scatterplots
# t-SNE
fig <- ggplot(df, aes(x=TSNE.1, y=TSNE.2, color=label)) + geom_point(size=2)
ggsave(paste(dir_o, paste0("scatterplot_tsne_", ds_n, ".pdf"), sep = "/"), plot = fig, 
       width = 12, height = 8, dpi = 300)

# PCA
sum(match(pca$Barcode, c_ann$barcode) == c(1:nrow(pca))) == nrow(pca)
df_2 <- cbind(pca, data.frame(label=c_ann$label))
fig <- ggplot(df_2, aes(x=PC.1, y=PC.2, color=label)) + geom_point(size=2)
ggsave(paste(dir_o, paste0("scatterplot_pca_", ds_n, ".pdf"), sep = "/"), plot = fig, 
       width = 12, height = 8, dpi = 300)
