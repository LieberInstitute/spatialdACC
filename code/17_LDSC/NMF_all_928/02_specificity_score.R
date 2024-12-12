setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

# Load libraries
library(here)

topN.mat <- readRDS(here("processed-data", "17_LDSC", "top928_intermediate_mat.rds"))

write.csv(topN.mat,here::here("processed-data", "17_LDSC", "NMF_score_all_928.csv"))
