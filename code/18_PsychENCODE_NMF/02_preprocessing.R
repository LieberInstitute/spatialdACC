# now we have the individual pseudobulked matrices
# we want to combine them into a single matrix


library(SingleCellExperiment)
library(spatialLIBD)
library(here)

pseudo_path <- here("processed-data", "18_PsychENCODE_NMF", "pseudobulk")
# find all files in this path
pseudo_files <- list.files(pseudo_path, full.names = TRUE)
sce_pseudo1 <- readRDS(pseudo_files[1])
sce_pseudo2 <- readRDS(pseudo_files[2])
sce_pseudo3 <- readRDS(pseudo_files[3])
sce_pseudo4 <- readRDS(pseudo_files[4])
sce_pseudo5 <- readRDS(pseudo_files[5])
sce_pseudo6 <- readRDS(pseudo_files[6])
sce_pseudo7 <- readRDS(pseudo_files[7])

# combine the logcounts matrices
logcounts1 <- logcounts(sce_pseudo1)
logcounts2 <- logcounts(sce_pseudo2)
logcounts3 <- logcounts(sce_pseudo3)
logcounts4 <- logcounts(sce_pseudo4)
logcounts5 <- logcounts(sce_pseudo5)
logcounts6 <- logcounts(sce_pseudo6)
logcounts7 <- logcounts(sce_pseudo7)

# problem that rownames are not the same
# just get overlaps between all datasets
rownames1 <- rownames(logcounts1)
rownames2 <- rownames(logcounts2)
rownames3 <- rownames(logcounts3)
rownames4 <- rownames(logcounts4)
rownames5 <- rownames(logcounts5)
rownames6 <- rownames(logcounts6)
rownames7 <- rownames(logcounts7)

# get overlaps
overlaps <- Reduce(intersect, list(rownames1, rownames2, rownames3, rownames4, rownames5, rownames6, rownames7))

# filter out non-overlapping genes
logcounts1 <- logcounts1[overlaps, ]
logcounts2 <- logcounts2[overlaps, ]
logcounts3 <- logcounts3[overlaps, ]
logcounts4 <- logcounts4[overlaps, ]
logcounts5 <- logcounts5[overlaps, ]
logcounts6 <- logcounts6[overlaps, ]
logcounts7 <- logcounts7[overlaps, ]

# check rownames in the same order
all(match(rownames(logcounts1), rownames(logcounts2)) == 1:nrow(logcounts1))
all(match(rownames(logcounts1), rownames(logcounts3)) == 1:nrow(logcounts1))
all(match(rownames(logcounts1), rownames(logcounts4)) == 1:nrow(logcounts1))
all(match(rownames(logcounts1), rownames(logcounts5)) == 1:nrow(logcounts1))
all(match(rownames(logcounts1), rownames(logcounts6)) == 1:nrow(logcounts1))
all(match(rownames(logcounts1), rownames(logcounts7)) == 1:nrow(logcounts1))

# combine the logcounts matrices
logcounts_combined <- cbind(logcounts1, logcounts2, logcounts3, logcounts4, logcounts5, logcounts6, logcounts7)

# create a SingleCellExperiment object
sce <- SingleCellExperiment(assays = list(logcounts = logcounts_combined))

# save
saveRDS(sce, file = here("processed-data", "18_PsychENCODE_NMF", "pseudobulk", "pseudobulk_combined.rds"))
