setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

# Load libraries
library(here)

topN.mat <- readRDS(here("processed-data", "17_LDSC", "top928_intermediate_mat.rds"))
cell_types <- sub("\\.NMF\\d+$", "", colnames(topN.mat))

result <- topN.mat %>%
    as.data.frame() %>%
    mutate(Gene = rownames(.)) %>%
    pivot_longer(-Gene, names_to = "Pattern", values_to = "Count") %>%
    mutate(CellType = sub("\\.NMF\\d+$", "", Pattern)) %>%
    group_by(Gene, CellType) %>%
    summarize(Sum = sum(Count), .groups = "drop") %>%
    pivot_wider(names_from = CellType, values_from = Sum, values_fill = 0)

result <- result %>%
    mutate(TotalOverlaps = rowSums(across(-Gene)),
           CategoryCount = rowSums(select(., -Gene) > 0))

result <- result %>%
    filter(TotalOverlaps > 1)

result <- result %>%
    mutate(Classification = ifelse(CategoryCount == 1, "Single Category", "Multiple Categories"))


dim(topN.mat)
remove_list <- which(result$Classification=="Multiple Categories")
topN.mat <- topN.mat[-remove_list,]
dim(topN.mat)

write.csv(topN.mat,here::here("processed-data", "17_LDSC", "NMF_score_928.csv"))
