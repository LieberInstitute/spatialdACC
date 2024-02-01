setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

library("SingleCellExperiment")
library("here")
library("scater")
library("scDblFinder")
library("BiocParallel")
library("tidyverse")

load(here("processed-data", "snRNA-seq", "01_QC", "sce_qc.rda"))

sce <- scDblFinder(sce, samples="Sample", BPPARAM=MulticoreParam(2,RNGseed=1234))

#save doublet scores sce
save(sce, file=here("processed-data", "snRNA-seq", "01_QC", "sce_doublet.rda"))

summary(sce$scDblFinder.score)
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
# 0.0000657 0.0010397 0.0012999 0.0568570 0.0030581 0.9999392

table(sce$scDblFinder.class)
# singlet doublet
# 35161    1871

table(sce$Sample, sce$scDblFinder.class)
#               singlet doublet
# 10c_dACC_SVB    4514     199
# 1c_dACC_MRV     3654     168
# 2c_dACC_MRV     3781     164
# 3c_dACC_MRV     3629     203
# 4c_dACC_MRV     3507     222
# 5c_dACC_SVB     2267     124
# 6c_dACC_SVB     2399     133
# 7c_dACC_SVB     3857     170
# 8c_dACC_SVB     3465     187
# 9c_dACC_SVB     4088     301

sum(sce$scDblFinder.score > 0.99)
# [1] 1421

# check if doublets are overlapping with low qc
#most of the high mito are singlets
table(sce$high_mito, sce$scDblFinder.class)
#        singlet doublet
# FALSE   35161    1871

#we just flagged doublet scores, did not remove from sce
