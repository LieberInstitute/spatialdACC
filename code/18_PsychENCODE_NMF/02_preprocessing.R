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

# get the counts matrices
counts1 <- counts(sce_pseudo1)
counts2 <- counts(sce_pseudo2)
counts3 <- counts(sce_pseudo3)
counts4 <- counts(sce_pseudo4)
counts5 <- counts(sce_pseudo5)
counts6 <- counts(sce_pseudo6)
counts7 <- counts(sce_pseudo7)

# combine the logcounts matrices
logcounts1 <- logcounts(sce_pseudo1)
logcounts2 <- logcounts(sce_pseudo2)
logcounts3 <- logcounts(sce_pseudo3)
logcounts4 <- logcounts(sce_pseudo4)
logcounts5 <- logcounts(sce_pseudo5)
logcounts6 <- logcounts(sce_pseudo6)
logcounts7 <- logcounts(sce_pseudo7)

# List of logcounts variables
logcounts_list <- list(logcounts1, logcounts2, logcounts3, logcounts4, logcounts5, logcounts6, logcounts7)

# Iterate over the logcounts list
for (i in seq_along(logcounts_list)) {
    # Get the current logcounts matrix
    logcounts_matrix <- logcounts_list[[i]]

    # Count the number of negative values
    num_negative_values <- sum(logcounts_matrix < 0)

    # Calculate the proportion of negative values
    total_values <- length(logcounts_matrix)
    proportion_negative <- num_negative_values / total_values

    # Print the results
    cat("Matrix", i, "\n")
    cat("Number of negative values:", num_negative_values, "\n")
    cat("Proportion of negative values:", proportion_negative, "\n\n")
}


# do the same for counts matrices
# List of counts variables
counts_list <- list(counts1, counts2, counts3, counts4, counts5, counts6, counts7)

# Iterate over the counts list
for (i in seq_along(counts_list)) {
    # Get the current counts matrix
    counts_matrix <- counts_list[[i]]

    # Count the number of negative values
    num_negative_values <- sum(counts_matrix == 0)

    # Calculate the proportion of negative values
    total_values <- length(counts_matrix)
    proportion_negative <- num_negative_values / total_values

    # Print the results
    cat("Matrix", i, "\n")
    cat("Number of zero values:", num_negative_values, "\n")
    cat("Proportion of zero values:", proportion_negative, "\n\n")
}

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

# filter out non-overlapping genes in counts
counts1 <- counts1[overlaps, ]
counts2 <- counts2[overlaps, ]
counts3 <- counts3[overlaps, ]
counts4 <- counts4[overlaps, ]
counts5 <- counts5[overlaps, ]
counts6 <- counts6[overlaps, ]
counts7 <- counts7[overlaps, ]

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

#combine the counts matrices
counts_combined <- cbind(counts1, counts2, counts3, counts4, counts5, counts6, counts7)

# combine the logcounts matrices
logcounts_combined <- cbind(logcounts1, logcounts2, logcounts3, logcounts4, logcounts5, logcounts6, logcounts7)

# create a SingleCellExperiment object
sce <- SingleCellExperiment(assays = list(counts = counts_combined))

# compute the logcounts
sce <- logNormCounts(sce)

# save
saveRDS(sce, file = here("processed-data", "18_PsychENCODE_NMF", "pseudobulk_combined.rds"))
