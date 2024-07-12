setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

# Load libraries
library(here)

is_top_10_percent <- function(column) {
    top_10_percent_threshold <- quantile(column, 0.9)
    as.numeric(column >= top_10_percent_threshold)
}

dat <- read.table(here::here("processed-data", "17-2_LDSC_enrichment", "snRNA-seq_aggregated_de.tsv"),header=T)

res <- apply(dat, 2, is_top_10_percent)
rownames(res) <- rownames(dat)
write.csv(res,here::here("processed-data", "17-2_LDSC_enrichment", "snRNA-seq_score.csv"))
