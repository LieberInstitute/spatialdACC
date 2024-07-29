setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

# Load libraries
library(here)

# there are a lot of NAs (esp WM-CC and L2) because there are some genes in those
# domains that are only downregulated and therefore their t stats go into the
# upregulated domain

is_top_10_percent <- function(column) {
    top_10_percent_threshold <- quantile(column, 0.9, na.rm = T)
    result <- as.numeric(column >= top_10_percent_threshold)
    result[is.na(result)] <- 0  # Replace NA values with 0
    return(result)
}

dat <- read.table(here::here("processed-data", "17-3_LDSC_pairwise", "visium_aggregated_de.tsv"),header=T)
# count NAs in each column
na_count <- colSums(is.na(dat))

dat.norm <- t(apply(dat, 1, function(x) {
    x / sum(x, na.rm = TRUE)
}))
res <- apply(dat.norm, 2, is_top_10_percent)
rownames(res) <- rownames(dat.norm)
write.csv(res,here::here("processed-data", "17-3_LDSC_pairwise", "visium_score.csv"))
