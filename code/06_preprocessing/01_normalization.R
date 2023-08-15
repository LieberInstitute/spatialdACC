setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')
suppressPackageStartupMessages(library("here"))
suppressPackageStartupMessages(library("scran"))

load(here("processed-data", "05_QC", "spe_QC.Rdata"), verbose = TRUE)

#remove 1 spot with 100% mito reads
index <- which(colData(spe)$subsets_Mito_percent > 99)
spe <- spe[,-index]

#discard extremely low detected & sum spots
spe$discard_extreme <- spe$extreme_low_sum | spe$extreme_low_detected
spe <- spe[,!colData(spe)$discard_extreme]

# calculate library size factors
spe <- computeLibraryFactors(spe)
summary(sizeFactors(spe))
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 0.00051  0.35785  0.76570  1.00000  1.31600 11.07102

pdf(here("plots", "06_preprocessing", "hist_size_factors.pdf"), width = 21, height = 10)
hist(sizeFactors(spe), breaks = 20)
dev.off()

spe <- logNormCounts(spe)

save(spe, file = here::here("processed-data", "06_preprocessing", "spe_norm.Rdata"))
