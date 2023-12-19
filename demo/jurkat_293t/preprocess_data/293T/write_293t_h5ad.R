#!/usr/bin/env Rscript
# write 293T dataset to a H5AD file

options(stringsAsFactors = F)
library(Seurat)
library(anndata)
library(Matrix)

# working directory
wd <- "/home/nikola/Project_Data/R_data/tests/Jurkat_293T/293T"
setwd(wd)

# dataset names
ds_n <- "293t"

# data directories
d_rt <- "/home/nikola/Documents/Rsch/resources/scRNA_seq/Jurkat_293T"
dir_dt <- paste(d_rt, ds_n, sep = "/")

# output directory
dir_o <- paste(wd, "output", sep = "/")
# dir.create(dir_o, recursive = T)

# 1 load data -------------------------------------------------------------

# read 10x file
mat <- Read10X(data.dir = dir_dt)
head(mat)[,1:6]
dim(mat)

# 2 create AnnData object -------------------------------------------------

# meta data
c_ann <- rep("293t", ncol(mat))
mt_dt <- data.frame(row.names = colnames(mat), label=c_ann)
head(mt_dt)

# feature data
ft_dt <- data.frame(row.names = rownames(mat), type=rep("gene", nrow(mat)))
head(ft_dt)

# transpose the matrix
mat_t <- t(mat)
dim(mat_t)

# convert dgCMatrix to dgRMatrix
mat_r <- as(as(as(mat_t, "dMatrix"), "generalMatrix"), "RsparseMatrix")

adt <- AnnData(X = mat_r, obs = mt_dt, var = ft_dt)

# 3 write data to H5AD file -----------------------------------------------

write_h5ad(adt, filename = paste(dir_o, paste0("mat_", ds_n, ".h5ad"), sep = "/"))
