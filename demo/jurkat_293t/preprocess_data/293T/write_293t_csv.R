#!/usr/bin/env Rscript
# write 293T dataset to CSV files

options(stringsAsFactors = F)
library(Seurat)

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

# read 10x files
mat <- Read10X(data.dir = dir_dt)

head(mat)[,1:6]
dim(mat)

# 2 write data to CSV files -----------------------------------------------

write.csv(mat, file = paste(dir_o, "mat_293t.csv", sep = "/"),
          quote = F, row.names = T)

# meta data
c_ann <- rep("293t", ncol(mat))
head(c_ann)

write.table(c_ann, file = paste(dir_o, "c_ann_293t.csv", sep = "/"),
          quote = F, sep = ",", row.names = F, col.names = F)
