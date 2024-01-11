setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

library("SingleCellExperiment")
library("here")
library("scater")
library("scran")
library("sessioninfo")
library("purrr")
library("tidyverse")

load(here("processed-data", "04_build_sce", "1c-10c_sce_raw.rda"))

##get droplet score filepaths
droplet_paths <- list.files(here("processed-data", "snRNA-seq", "01_QC"),
                            full.names = TRUE
)

names(droplet_paths) <- gsub("st", "s", gsub("droplet_scores_|.Rdata", "", basename(droplet_paths)))

e.out <- lapply(droplet_paths, function(x) get(load(x)))

#### Compile drop empty info ####
logs <- list.files(here("code", "snRNA-seq", "01_QC", "logs"), pattern = "drops_k.[0-9]", full.names = TRUE)
logs <- map(logs, readLines)

knee_lower <- map_dbl(logs, ~ parse_number(.x[grepl("knee_lower =", .x)]))
names(knee_lower) <- gsub("st", "s", map_chr(logs, ~str_sub(.x[grepl("Running Sample: ", .x)], " ", 3)))

FDR_cutoff <- 0.001

drop_summary <- stack(map_int(e.out, nrow)) %>%
    rename(total_n = values) %>%
    left_join(stack(map_int(e.out, ~ sum(.x$FDR < FDR_cutoff, na.rm = TRUE))) %>%
                  rename(non_empty = values)) %>%
    select(Sample = ind, total_n, non_empty) %>%
    left_join(stack(knee_lower) %>% rename(Sample = ind, lower_cutoff = values))

write_csv(drop_summary, file = here("processed-data", "snRNA-seq", "01_QC", "drop_summary.csv"))

drop_summary %>%
    arrange(non_empty)

summary(drop_summary$non_empty)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 2765    4217    4404    4276    4582    5569

##which sample is problematic?
drop_summary$Sample[which.max(drop_summary$non_empty)]
# [1] 10c_dACC_SVB

#make Louise-style barplot
drop_barplot <- drop_summary %>%
    mutate(empty = total_n - non_empty) %>%
    select(-total_n) %>%
    pivot_longer(!Sample, names_to = "drop_type", values_to = "n_drop") %>%
    ggplot(aes(x = Sample, y = n_drop, fill = drop_type)) +
    geom_col() +
    scale_y_continuous(trans = "log10") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(drop_barplot, filename = here("plots", "snRNA-seq", "01_QC", "drop_barplot.png"), width = 9)


#### Eliminate empty droplets ####
dim(sce)
#[1]    36601 19473661

e.out.all <- do.call("rbind", e.out)[colnames(sce), ]
sce <- sce[, which(e.out.all$FDR <= 0.001)]

dim(sce)
# [1]    36601 42764

#save drops removed sce
save(sce,file=here("processed-data", "snRNA-seq", "01_QC", "sce_drops_removed.rda"))


## Remove genes with no data
no_expr <- which(rowSums(counts(sce)) == 0)
length(no_expr)
# [1] 1606
length(no_expr) / nrow(sce) * 100
# [1] 4.387858
sce <- sce[-no_expr, ]
dim(sce)
# [1]    34995 19473661

## Remove spots without counts
if (any(colSums(counts(sce)) == 0)) {
    message("removing cells without counts for sce")
    sce <- sce[, -which(colSums(counts(sce)) == 0)]
    dim(sce)
}

#[1]    34995 14472681

#explore mitochondrial read percents
sce <- scuttle::addPerCellQC(
    sce,
    subsets = list(Mito = which(seqnames(sce) == "chrM")))

fivenum(colData(sce)$subsets_Mito_percent)
# [1]   0   0   0   0 100

#### Check for low quality spots ####

## High mito
sce$high_mito <- isOutlier(sce$subsets_Mito_percent, nmads = 3, type = "higher", batch = sce$Sample)
table(sce$high_mito)
# FALSE  TRUE
# 13417507  1055174

table(sce$Sample, sce$high_mito)
#FALSE    TRUE
#10c_dACC_SVB 1275052  108720
#1c_dACC_MRV  1365858   94307
#2c_dACC_MRV  1372737  115198
#3c_dACC_MRV  1293772  104122
#4c_dACC_MRV  1241509  123297
#5c_dACC_SVB  1198363  106234
#6c_dACC_SVB  1320320  120659
#7c_dACC_SVB  1603770   67412
#8c_dACC_SVB  1474837  110440
#9c_dACC_SVB  1271289  104785


## low library size
sce$low_sum <- isOutlier(sce$sum, log = TRUE, nmads = 3, type = "lower", batch = sce$Sample)
table(sce$low_sum)
# FALSE
# 14472681

## low detected features
sce$low_detected <- isOutlier(sce$detected, log = TRUE, nmads = 3, type = "lower", batch = sce$Sample)
table(sce$low_detected)
# FALSE
# 14472681

#### QC plots ####
pdf(here("plots", "snRNA-seq", "01_QC", "QC_violin_plots.pdf"), width = 21)

hist(colData(sce)$subsets_Mito_percent, breaks = 400, ylim = c(0,2000))

plotColData(sce, x = "Sample", y = "subsets_Mito_percent", colour_by = "high_mito") +
    ggtitle("Mito Percent")

plotColData(sce, x = "Sample", y = "sum", colour_by = "low_sum") +
    scale_y_log10() +
    ggtitle("Total UMIs")

plotColData(sce, x = "Sample", y = "detected", colour_by = "low_detected") +
    scale_y_log10() +
    ggtitle("Detected genes")

dev.off()
