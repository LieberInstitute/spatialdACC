setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

library("SingleCellExperiment")
library("here")
library("scater")
library("scDblFinder")
library("BiocParallel")
library("tidyverse")

load(here("processed-data", "snRNA-seq", "01_QC", "sce_qc.rda"))

sce <- scDblFinder(sce, samples="Sample", BPPARAM=MulticoreParam(8,RNGseed=1234))

#save doublet scores sce
save(sce, file=here("processed-data", "snRNA-seq", "01_QC", "sce_doublet.rda"))

summary(sce$scDblFinder.score)
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
# 0.0002785 0.0006567 0.0010397 0.0644836 0.0031831 0.9997794

table(sce$Sample, sce$scDblFinder.class)
#               singlet doublet
#  10c_dACC_SVB    5312     257
# 1c_dACC_MRV     4206     229
# 2c_dACC_MRV     4406     225
# 3c_dACC_MRV     4068     260
# 4c_dACC_MRV     4035     345
# 5c_dACC_SVB     2606     159
# 6c_dACC_SVB     2836     160
# 7c_dACC_SVB     4218     210
# 8c_dACC_SVB     3941     239
# 9c_dACC_SVB     4658     394

sum(sce$scDblFinder.score > 0.99)
# [1] 1745

# check if doublets are overlapping with low qc
#most of the high mito are singlets
#none of the low sum/detected are doublets
table(sce$high_mito, sce$scDblFinder.class)
#        singlet doublet
# FALSE   34831    2288
# TRUE     5455     190

table(sce$low_sum, sce$scDblFinder.class)
#         singlet doublet
# FALSE   40002    2478
# TRUE      284       0

table(sce$low_detected, sce$scDblFinder.class)
#         singlet doublet
# FALSE   39930    2478
# TRUE      356       0

#we just flagged doublet scores, did not remove from sce
