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
# 0.0000381 0.0004016 0.0009069 0.0609606 0.0018350 0.9999698

table(sce$scDblFinder.class)
# singlet doublet
# 39849    2325

table(sce$Sample, sce$scDblFinder.class)
#               singlet doublet
#  10c_dACC_SVB    5109     246
# 1c_dACC_MRV     4186     241
# 2c_dACC_MRV     4329     238
# 3c_dACC_MRV     4084     229
# 4c_dACC_MRV     4022     269
# 5c_dACC_SVB     2631     131
# 6c_dACC_SVB     2680     158
# 7c_dACC_SVB     4216     178
# 8c_dACC_SVB     3946     232
# 9c_dACC_SVB     4646     403

sum(sce$scDblFinder.score > 0.99)
# [1] 1666

# check if doublets are overlapping with low qc
#most of the high mito are singlets
table(sce$high_mito, sce$scDblFinder.class)
#        singlet doublet
# FALSE   34876    2156
# TRUE     4973     169

#we just flagged doublet scores, did not remove from sce
