setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')
suppressPackageStartupMessages(library("SpatialExperiment"))
suppressPackageStartupMessages(library("scuttle"))
suppressPackageStartupMessages(library("scran"))
suppressPackageStartupMessages(library("scater"))
suppressPackageStartupMessages(library("jaffelab"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("here"))
suppressPackageStartupMessages(library("sessioninfo"))

load(here("processed-data", "VSPG","01_QC", "spe_basic.Rdata"), verbose = TRUE)

dim(spe)
# [1]  28096 17349
length(table(spe$sample_id))
# 4
length(table(spe$brnum))
# 4

spe <- scuttle::addPerCellQC(
    spe,
    subsets = list(Mito = which(seqnames(spe) == "chrM")),
    BPPARAM = BiocParallel::MulticoreParam(4)
)

#### Check for low quality spots ####

## High mito
spe$high_mito_br <- isOutlier(spe$subsets_Mito_percent, nmads = 3, type = "higher", batch = spe$brnum)

table(spe$high_mito_br)
# FALSE  TRUE
# 17028   321

table(spe$sample_id, spe$high_mito_br)
#   FALSE TRUE
#   V12N28-333_A1  4174   75
#   V12N28-333_B1  4274   39
#   V12N28-333_C1  4401   85
#   V12N28-333_D1  4179  122

## low library size
spe$low_sum_br <- isOutlier(spe$sum, log = TRUE, type = "lower", batch = spe$brnum)
table(spe$low_sum_br)
# FALSE   TRUE
# 17301    48

## low detected features
spe$low_detected_br <- isOutlier(spe$detected, log = TRUE, type = "lower", batch = spe$brnum)
table(spe$low_detected_br)
# FALSE   TRUE
# 17279    70

## are all low sum are also low detected? Mostly
table(spe$low_sum_br, spe$low_detected_br)
#      FALSE   TRUE
#FALSE 17278    23
#TRUE      1    47

## are all low sum are also high mito? No
table(spe$low_sum_br, spe$high_mito_br)
#      FALSE   TRUE
#FALSE 16985   316
#TRUE     43     5

## Annotate spots to drop
spe$discard_auto_br <- spe$high_mito_br | spe$low_sum_br | spe$low_detected_br

table(spe$discard_auto_br)
# FALSE   TRUE
# 16965   384

## discard 4.4% of spots by brain
100 * sum(spe$discard_auto_br) / ncol(spe)
# 2.213384

save(spe, file = here::here("processed-data", "VSPG","01_QC", "spe_QC.Rdata"))
