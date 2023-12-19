#!/usr/bin/env Rscript
# write Jurkat_293T_50_50 dataset to H5AD files

options(stringsAsFactors = F)
library(Seurat)
library(anndata)
library(Matrix)

# working directory
wd <- "/home/nikola/Project_Data/R_data/tests/Jurkat_293T/Jurkat_293T_50_50"
setwd(wd)

# dataset names
ds_n <- "293t"
ds_n_2 <- "jurkat_293t_50_50"
ds_n_3 <- "jurkat"
ds_n_4 <- "jurkat_293t_99_1"

# data directories
d_rt <- "/home/nikola/Documents/Rsch/resources/scRNA_seq/Jurkat_293T"
dir_dt <- paste(d_rt, ds_n, sep = "/")
dir_dt_2 <- paste(d_rt, ds_n_2, sep = "/")
dir_mt <- paste(d_rt, "cell_labels", sep = "/")
dir_dt_3 <- paste(d_rt, ds_n_3, sep = "/")
dir_dt_4 <- paste(d_rt, ds_n_4, sep = "/")

# output directory
dir_o <- paste(wd, "output", sep = "/")
# dir.create(dir_o, recursive = T)

# 1 load data -------------------------------------------------------------

# read 10x files
mat_293t <- Read10X(data.dir = dir_dt)
mat_50_50 <- Read10X(data.dir = dir_dt_2)

mat_jk <- Read10X(data.dir = dir_dt_3)
# mat_99_1 <- Read10X(data.dir = dir_dt_4)

head(mat_293t)[,1:6]
head(mat_50_50)[,1:6]
head(mat_jk)[,1:6]
# head(mat_99_1)[,1:6]

dim(mat_293t)
dim(mat_50_50)
dim(mat_jk)
# dim(mat_99_1)

# meta data
# 3258 jurkat cells
c_ann <- read.table(paste(dir_mt, "293t_jurkat_cluster_mdf.txt", sep = "/"), 
                    header = F, sep = "\t")
head(c_ann)
dim(c_ann)
# [1] 9531    1

ncol(mat_293t) + ncol(mat_jk) + ncol(mat_50_50)
# [1] 9531

c_ann[(ncol(mat_293t)-2):(ncol(mat_293t)+2),]
# [1] "293t"   "293t"   "293t"   "jurkat" "jurkat"

c_ann[(ncol(mat_293t)+ncol(mat_jk)-5):(ncol(mat_293t)+ncol(mat_jk)),]
# [1] "jurkat" "jurkat" "jurkat" "jurkat" "jurkat" "jurkat"

c_ann[(ncol(mat_293t)+ncol(mat_jk)+ncol(mat_50_50)-5):nrow(c_ann),]
# [1] "jurkat" "jurkat" "293t"   "293t"   "293t"   "jurkat"

length(c_ann[(nrow(c_ann)-ncol(mat_50_50)+1):nrow(c_ann),])
# [1] 3388

# c_ann
# ncol(mat_293t) (2885) 293_t cells
# ncol(mat_jk) (3258) jurkat cells
# ncol(mat_50_50) (3388) jurkat_293t_50_50 cells

# 2 create AnnData object -------------------------------------------------

# meta data
c_ann_50_50 <- c_ann[(nrow(c_ann)-ncol(mat_50_50)+1):nrow(c_ann),]
mt_dt <- data.frame(row.names = colnames(mat_50_50), label=c_ann_50_50)
head(mt_dt)

# feature data
ft_dt <- data.frame(row.names = rownames(mat_50_50), type=rep("gene", nrow(mat_50_50)))
head(ft_dt)

# transpose the matrix
mat_t <- t(mat_50_50)
dim(mat_t)

# convert dgCMatrix to dgRMatrix
mat_r <- as(as(as(mat_t, "dMatrix"), "generalMatrix"), "RsparseMatrix")

adt <- AnnData(X = mat_r, obs = mt_dt, var = ft_dt)

# 3 write data to H5AD file -----------------------------------------------

write_h5ad(adt, filename = paste(dir_o, paste0("mat_", ds_n_2, ".h5ad"), sep = "/"))
