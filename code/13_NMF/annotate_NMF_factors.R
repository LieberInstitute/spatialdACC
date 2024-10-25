setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')
library(SpatialExperiment)
library(scater)
library(RcppML)

# get NMF results from single nucleus data
x <- readRDS(file = here("processed-data", "snRNA-seq", "06_NMF", "nmf_results.RDS"))

patterns <- x@w

# we want to subset the patterns to only include some patterns
# we also want to rename the patterns to be more informative, such as "Astro - NMF32"
# Oligo: 26, 23, 27, 13, 43, 40, 36, 28, 9, 33, 39
# L5_6_NP: 46
# L6_b: 35
# L5_ET: 61
# L6_CT: 15
# L6_IT_Car3: 68
# L2_3_IT: 3, 11
# L5_IT: 38
# L6_IT: 32
# Pvalb: 10, 63
# SST: 52, 56
# SST Chodl: 51
# LAMP5: 37, 60
# Sncg: 55, 58
# Vip: 44, 47
# Endo: 75, 49
# Astro: 14, 21, 53, 65
# OPC: 17, 24
# VLMC: 59, 17
# microPVM: 19, 54, 57
# misc: 59, 64, 61, 15

# subset the patterns
patterns <- patterns[, c(26, 23, 27, 13, 43, 40, 36, 28, 9, 33, 39, 46, 35, 61, 15, 68, 3, 11, 38, 32, 10, 63, 52, 56, 51, 37, 60, 55, 58, 44, 47, 75, 49, 14, 21, 53, 65, 17, 24, 59, 17, 19, 54, 57, 59, 64, 61, 15)]

# rename the patterns
colnames(patterns) <- c("Oligo-NMF26", "Oligo-NMF23", "Oligo-NMF27", "Oligo-NMF13", "Oligo-NMF43", "Oligo-NMF40", "Oligo-NMF36", "Oligo-NMF28", "Oligo-NMF9", "Oligo-NMF33", "Oligo-NMF39", "L5_6_NP-NMF46", "L6_b-NMF35", "L5_ET-NMF61", "L6_CT-NMF15", "L6_IT_Car3-NMF68", "L2_3_IT-NMF3", "L2_3_IT-NMF11", "L5_IT-NMF38", "L6_IT-NMF32", "Pvalb-NMF10", "Pvalb-NMF63", "SST-NMF52", "SST-NMF56", "SST Chodl-NMF51", "LAMP5-NMF37", "LAMP5-NMF60", "Sncg-NMF55", "Sncg-NMF58", "Vip-NMF44", "Vip-NMF47", "Endo-NMF75", "Endo-NMF49", "Astro-NMF14", "Astro-NMF21", "Astro-NMF53", "Astro-NMF65", "OPC-NMF17", "OPC-NMF24", "VLMC-NMF59", "VLMC-NMF17", "microPVM-NMF19", "microPVM-NMF54", "microPVM-NMF57", "misc-NMF59", "misc-NMF64", "misc-NMF61", "misc-NMF15")

# save the patterns
saveRDS(patterns, file = here("processed-data", "snRNA-seq", "06_NMF", "nmf_patterns_subset.RDS"))
