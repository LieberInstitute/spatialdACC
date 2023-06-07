setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')
suppressPackageStartupMessages(library("SpatialExperiment"))
suppressPackageStartupMessages(library("scuttle"))
suppressPackageStartupMessages(library("scran"))
suppressPackageStartupMessages(library("scater"))
suppressPackageStartupMessages(library("jaffelab"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("here"))
suppressPackageStartupMessages(library("sessioninfo"))

load(here("processed-data", "05_QC", "spe_basic.Rdata"), verbose = TRUE)

dim(spe)
# [1]  29720 77604
length(table(spe$sample_id))
# 17
length(table(spe$brnum))
# 10

spe <- scuttle::addPerCellQC(
    spe,
    subsets = list(Mito = which(seqnames(spe) == "chrM")),
    BPPARAM = BiocParallel::MulticoreParam(4)
)

#### Check for low quality spots ####

## High mito
spe$high_mito_id <- isOutlier(spe$subsets_Mito_percent, nmads = 3, type = "higher", batch = spe$sample_id)
spe$high_mito_br <- isOutlier(spe$subsets_Mito_percent, nmads = 3, type = "higher", batch = spe$brnum)
table(spe$high_mito_id)
# FALSE  TRUE
# 77203   401

table(spe$high_mito_br)
# FALSE  TRUE
# 77384   220

table(spe$sample_id, spe$high_mito_id)
# FALSE TRUE
# V12J03-002_A1  4972    6
# V12J03-002_B1  4699    2
# V12J03-002_C1  4584    2
# V12N28-331_A1  4870   83
# V12N28-331_B1  4652   35
# V12N28-331_C1  4433   69
# V12N28-331_D1  4913   61
# V12N28-332_A1  4991    1
# V12N28-332_B1  4413   17
# V12N28-332_C1  4564   12
# V12N28-332_D1  4940    0
# V12N28-334_A1  4561   21
# V12N28-334_B1  4108   43
# V12N28-334_C1  4632   13
# V12N28-334_D1  4447   24
# V12Y31-080_B1  3908    0
# V12Y31-080_C1  3516   12

table(spe$high_mito_br, spe$brnum)
# Br2720 Br2743 Br3942 Br6423 Br6432 Br6471 Br6522 Br8325 Br8492 Br8667
# FALSE  14909   4632   4561   4447  13730  13659   9914   3516   4108   3908
# TRUE      14     13     21     24     88      5      0     12     43      0

#are high_mito_br spots and high_mito_id spots same? No
table(spe$high_mito_br,spe$high_mito_id)
#       FALSE  TRUE
# FALSE 77122   262
# TRUE     81   139

## low library size
spe$low_sum_id <- isOutlier(spe$sum, log = TRUE, type = "lower", batch = spe$sample_id)
table(spe$low_sum_id)
# FALSE   TRUE
# 76540  1064
spe$low_sum_br <- isOutlier(spe$sum, log = TRUE, type = "lower", batch = spe$brnum)
table(spe$low_sum_br)
# FALSE   TRUE
# 75108  2496

#are low_sum_br spots and low_sum_id spots same? no
table(spe$low_sum_br,spe$low_sum_id)
#       FALSE  TRUE
# FALSE 74960   148
# TRUE   1580   916

## low detected features
spe$low_detected_id <- isOutlier(spe$detected, log = TRUE, type = "lower", batch = spe$sample_id)
table(spe$low_detected_id)
# FALSE   TRUE
# 75836  1768
spe$low_detected_br <- isOutlier(spe$detected, log = TRUE, type = "lower", batch = spe$brnum)
table(spe$low_detected_br)
# FALSE   TRUE
# 74345  3259

#are low_detected_br spots and low_detected_id spots same? no
table(spe$low_detected_br,spe$low_detected_id)
#      FALSE   TRUE
#FALSE 74113   232
#TRUE   1723  1536

## are all low sum are also low detected? Mostly
table(spe$low_sum_br, spe$low_detected_br)
#      FALSE   TRUE
#FALSE 74345   763
#TRUE      0  2496
table(spe$low_sum_id, spe$low_detected_id)
#      FALSE   TRUE
#FALSE 75829   711
#TRUE      7  1057

## are all low sum are also high mito? No
table(spe$low_sum_br, spe$high_mito_br)
#      FALSE   TRUE
#FALSE 74943   165
#TRUE   2441    55
table(spe$low_sum_id, spe$high_mito_id)
#      FALSE   TRUE
#FALSE 76163   377
#TRUE   1040    24

## Annotate spots to drop
spe$discard_auto_br <- spe$high_mito_br | spe$low_sum_br | spe$low_detected_br
spe$discard_auto_id <- spe$high_mito_id | spe$low_sum_id | spe$low_detected_id

table(spe$discard_auto_br)
# FALSE   TRUE
# 74182  3422
table(spe$discard_auto_id)
# FALSE   TRUE
# 75464  2140

## discard 4.4% of spots by brain
100 * sum(spe$discard_auto_br) / ncol(spe)
# 4.409567
## discard 2.8% of spots by capture area
100 * sum(spe$discard_auto_id) / ncol(spe)
# 2.75759

save(spe, file = here::here("processed-data", "05_QC", "spe_QC.Rdata"))

## Reproducibility information
Sys.time()
proc.time()
options(width = 120)
session_info()
