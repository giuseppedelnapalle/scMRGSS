#!/usr/bin/env Rscript
# extract and preprocess scRNA-seq data of 293T

options(stringsAsFactors = F)
library(Seurat)

# working directory
wd <- "/home/nikola/Project_Data/R_data/tests/Jurkat_293T/293T"
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

# function directory
dir_f <- "/home/nikola/Project_Data/R_data/functions_scripts/h5"

# load functions
source(paste(dir_f, "h5_utilities.R", sep = "/"))

# 1 load data -------------------------------------------------------------

# read 10x files
matrix_sg <- Read10X(data.dir = dir_dt)
matrix_cmb <- Read10X(data.dir = dir_dt_2)

matrix_jk <- Read10X(data.dir = dir_dt_3)
matrix_99_1 <- Read10X(data.dir = dir_dt_4)

head(matrix_sg)[,1:6]
head(matrix_cmb)[,1:6]
head(matrix_jk)[,1:6]
head(matrix_99_1)[,1:6]

dim(matrix_sg)
dim(matrix_cmb)
dim(matrix_jk)
dim(matrix_99_1)

# meta data
c_ann <- read.table(paste(dir_mt, "293t_jurkat_cluster.txt", sep = "/"), 
                    header = F, sep = "\t")
head(c_ann)
dim(c_ann)
# [1] 9530    1

ncol(matrix_sg) + ncol(matrix_jk) + ncol(matrix_cmb)
# [1] 9531

c_ann[(ncol(matrix_sg)-2):(ncol(matrix_sg)+2),]
# [1] "293t"   "293t"   "293t"   "jurkat" "jurkat"

c_ann[(ncol(matrix_sg)+ncol(matrix_jk)-5):(ncol(matrix_sg)+ncol(matrix_jk)),]
# [1] "jurkat" "jurkat" "jurkat" "jurkat" "jurkat" "293t"

c_ann[(ncol(matrix_sg)+ncol(matrix_jk)+ncol(matrix_cmb)-5):nrow(c_ann),]
# [1] "jurkat" "293t"   "293t"   "293t"   "jurkat"

length(c_ann[(nrow(c_ann)-ncol(matrix_cmb)+1):nrow(c_ann),])

# c_ann
# ncol(matrix_sg) (2885) 293_t cells
# ncol(matrix_jk)-1 (3257) jurkat cells?
# ncol(matrix_cmb) (3388) jurkat_293t_50_50 cells

# 2 preprocessing ---------------------------------------------------------

# calculate n_features
n_f <- calc_n_features(matrix_sg)
n_f_2 <- calc_n_features(matrix_cmb)

summary(n_f)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1463    3041    3393    3416    3792    5784
summary(n_f_2)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1454    2979    3409    3407    3832    6071

# calculate n_counts
n_c <- calc_n_counts(matrix_sg)
n_c_2 <- calc_n_counts(matrix_cmb)

summary(n_c)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 3890   11443   14338   15270   17942   52302
summary(n_c_2)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 3799   10703   13987   15041   17973   55987

# filter cells based on number of features
matrix_sg_f <- filter_cells_n_f(matrix_sg, 2000)
dim(matrix_sg_f)

matrix_cmb_f <- filter_cells_n_f(matrix_cmb, 2000)
dim(matrix_cmb_f)

# filter genes
mat_lst <- filter_match_genes(matrix_sg_f, matrix_cmb_f, 3)
matrix_sg_f <- mat_lst[[1]]
matrix_cmb_f <- mat_lst[[2]]
dim(matrix_sg_f)
dim(matrix_cmb_f)

# 3 extract 293T data -----------------------------------------------------

# cell labels for jurkat_293t_50_50 dataset
# length(c_ann[(nrow(c_ann)-ncol(matrix_cmb)+1):nrow(c_ann),])
c_ann_cmb <- data.frame(barcode=colnames(matrix_cmb), 
                        label=c_ann[(nrow(c_ann)-ncol(matrix_cmb)+1):nrow(c_ann),])
head(c_ann_cmb)
dim(c_ann_cmb)

saveRDS(c_ann_cmb, file = paste(dir_o, "c_ann_jurkat_293t_50_50.rds", sep = "/"))

# filtered
c_ann_cmb_f <- c_ann_cmb[c_ann_cmb$barcode %in% colnames(matrix_cmb_f),]
dim(c_ann_cmb_f)
saveRDS(c_ann_cmb_f, file = paste(dir_o, "c_ann_jurkat_293t_50_50_filtered.rds", sep = "/"))

select <- colnames(matrix_cmb_f) %in% c_ann_cmb_f$barcode[c_ann_cmb_f$label == "293t"]
matrix_sg_ext <- matrix_cmb_f[, select]
dim(matrix_sg_ext)

write.csv(matrix_sg_f, file = paste(dir_o, "mat_293t_f.csv", sep = "/"),
          quote = F, row.names = T)
write.csv(matrix_sg_ext, file = paste(dir_o, "mat_293t_50_50_mix_f.csv", sep = "/"),
          quote = F, row.names = T)

# save objects and image --------------------------------------------------

rm(matrix_sg, matrix_cmb, matrix_jk, matrix_99_1)
rm(matrix_sg_f, matrix_cmb_f, matrix_sg_ext)
rm(mat_lst)

save.image(file = "data_prep_293T.RData")
