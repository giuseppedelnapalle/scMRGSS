#!/usr/bin/env Rscript
# write Jurkat dataset to CSV files

options(stringsAsFactors = F)
library(Seurat)

# working directory
wd <- "/home/nikola/Project_Data/R_data/tests/Jurkat_293T/Jurkat"
setwd(wd)

# dataset names
ds_n <- "jurkat"

# data directories
d_rt <- "/home/nikola/Documents/Rsch/resources/scRNA_seq/Jurkat_293T"
dir_dt <- paste(d_rt, ds_n, sep = "/")

# output directory
dir_o <- paste(wd, "output", sep = "/")
# dir.create(dir_o, recursive = T)

# # function directory
# dir_f <- "/home/nikola/Project_Data/R_data/functions_scripts/h5"
# 
# # load functions
# source(paste(dir_f, "h5_utilities.R", sep = "/"))

# 1 load data -------------------------------------------------------------

# read 10x files
matrix_jk <- Read10X(data.dir = dir_dt)

head(matrix_jk)[,1:6]
dim(matrix_jk)

# 2 write data to CSV files -----------------------------------------------

write.csv(matrix_jk, file = paste(dir_o, "mat_jurkat.csv", sep = "/"),
          quote = F, row.names = T)

# meta data
c_ann <- rep("jurkat", ncol(matrix_jk))
head(c_ann)

write.table(c_ann, file = paste(dir_o, "c_ann_jurkat.csv", sep = "/"),
          quote = F, sep = ",", row.names = F, col.names = F)
