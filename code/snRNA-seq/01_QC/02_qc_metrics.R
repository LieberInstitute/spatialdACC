setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

library("SingleCellExperiment")
library("here")
library("scater")
library("scran")
library("sessioninfo")
library("purrr")
library("tidyverse")
library("scDblFinder")
library("jaffelab")
library("batchelor")
library("patchwork")

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

#boxplot of distribution of empty droplets
drop_boxplot <- drop_summary %>%
    ggplot(aes(x="", y=non_empty)) +
    geom_boxplot()  +
    stat_summary(
        aes(label = round(stat(y), 1)),
        geom = "text",
        fun.y = function(y) { o <- boxplot.stats(y)$out; if(length(o) == 0) NA else o },
        hjust = -1
    ) +
    ggtitle("Distribution of Non-empty Droplets") +
    xlab("") +
    ylab("Number of Non-empty droplets")

ggsave(drop_boxplot, filename = here("plots", "snRNA-seq", "01_QC", "drop_boxplot.png"), width = 9)


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
# [1] 1735
length(no_expr) / nrow(sce) * 100
# [1] 4.740308
sce <- sce[-no_expr, ]
dim(sce)
# [1]   34866 42764

## Remove nuclei without counts
if (any(colSums(counts(sce)) == 0)) {
    message("removing cells without counts for sce")
    sce <- sce[, -which(colSums(counts(sce)) == 0)]
    dim(sce)
}

#[1]    34866 42764

#explore mitochondrial read percents
sce <- scuttle::addPerCellQC(
    sce,
    subsets = list(Mito = which(seqnames(sce) == "chrM")))

fivenum(colData(sce)$subsets_Mito_percent)
# [1]  0.00000000  0.04421320  0.09337795  0.20525799 16.99633700

#### Check for low quality nuclei ####

## High mito
sce$high_mito <- isOutlier(sce$subsets_Mito_percent, nmads = 3, type = "higher", batch = sce$Sample)
table(sce$high_mito)
# FALSE  TRUE
# 37119  5645

table(sce$Sample, sce$high_mito)
#FALSE    TRUE
#10c_dACC_SVB  4750  819
#1c_dACC_MRV   3822  613
#2c_dACC_MRV   3945  686
#3c_dACC_MRV   3832  496
#4c_dACC_MRV   3729  651
#5c_dACC_SVB   2391  374
#6c_dACC_SVB   2572  424
#7c_dACC_SVB   4037  391
#8c_dACC_SVB   3652  528
#9c_dACC_SVB   4389  663

sum(sce$subsets_Mito_percent > 0.38)
# [1] 5648

## low library size
sce$low_sum <- isOutlier(sce$sum, log = TRUE, nmads = 3, type = "lower", batch = sce$Sample)
table(sce$low_sum)
# FALSE  TRUE
# 42480   284

sum(sce$sum < 540)
# [1] 285

## low detected features
sce$low_detected <- isOutlier(sce$detected, log = TRUE, nmads = 3, type = "lower", batch = sce$Sample)
table(sce$low_detected)
# FALSE  TRUE
# 42408   356

sum(sce$detected < 500)
# [1] 368

## All low sum are also low detected
table(sce$low_sum, sce$low_detected)
#         FALSE   TRUE
# FALSE 42408    72
# TRUE      0   284

sce$discard_auto <- sce$high_mito | sce$low_sum | sce$low_detected

table(sce$discard_auto)
# FALSE   TRUE
# 37032  5732

#### QC plots ####

sce$`Mito Discard` <- sce$high_mito
sce$`Detected Discard` <- sce$low_detected
sce$`Combined Discard` <- sce$discard_auto

p1 <- plotColData(sce, x = "Sample", y = "subsets_Mito_percent", colour_by = "Mito Discard") +
    ggtitle("Mitochondrial Percent") +
    ylab("Mitochondrial Percent")

p2 <- plotColData(sce, x = "Sample", y = "detected", colour_by = "Detected Discard") +
    scale_y_log10() +
    ggtitle("Detected Genes") +
    ylab("Number of Detected Genes")

p3 <- plotColData(sce, x = "Sample", y = "detected", colour_by = "Combined Discard") +
    scale_y_log10() +
    ggtitle("Total Discarded Genes") +
    ylab("Number of Detected Genes")


png(here("plots", "snRNA-seq", "01_QC", "QC_violin_plots.png"), width = 15, height = 15, units = "in", res=300)
wrap_plots(p1,p2,p3,nrow=3) + plot_annotation(tag_levels = 'A')
dev.off()

# remove low quality nuclei
sce <- sce[, !colData(sce)$discard_auto]
dim(sce)
# [1] 34866 42408

#save drops removed and qc removed sce
save(sce,file=here("processed-data", "snRNA-seq", "01_QC", "sce_qc.rda"))

