library(SingleCellExperiment)
library(zellkonverter)
library(spatialLIBD)
library(here)

# need R/4.4.x to use zell konverter

# read in meta data
# Read the tab-delimited file into R
url <- "https://brainscope.gersteinlab.org/data/sample_metadata/PEC2_sample_metadata.txt"
data <- read.delim(url, header = TRUE)

# paths to h5ad files
file_path <- "/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/psychENCODE/version5/"

data_list <- c("CMC_annotated.h5ad", "DevBrain-snRNAseq_annotated.h5ad",
               "IsoHuB_annotated.h5ad",
               "MultiomeBrain-DLPFC_annotated.h5ad", "PTSDBrainomics_annotated.h5ad",
               "SZBDMulti-Seq_annotated.h5ad", "UCLA-ASD_annotated_mismatches_removed.h5ad")

# CMC
dataset <- "CMC"
file_name = "CMC_annotated.h5ad"
file <- paste0(file_path, file_name)
sce <- readH5AD(file)

colData(sce)$cellType <- as.factor(make.names(colData(sce)$subclass))
table(sce$cellType)

rowData(sce)$gene_name <- rownames(sce) # save gene name as column of rowData
rownames(sce) <- rowData(sce)$featureid # have to make row names of object the ensembl id instead of gene names

## Logcounts
# default “X” contain the log-normalized counts
message(Sys.time(), " revert to counts")

## check for all 0s (just first 100 cols for mem)
stopifnot(any(assays(sce)$X[, 1:100] != 0))

counts(sce) <- assays(sce)$X # change to normalized counts
counts(sce)@x <- 2^(counts(sce)@x) - 1 ## remove log2(counts + 1)

# filter out Cohort = CMC and Disorder = control from metadata
# only get Individual_ID for CMC and Disorder = control
data_sub <- data[data$Cohort == "CMC" & data$Disorder == "control", c("Individual_ID", "Cohort", "Disorder")]

# only keep cells that have colData(sce)$individualID in data_sub$Individual_ID
sce <- sce[, colData(sce)$individualID %in% data_sub$Individual_ID]

#### Pseudobulk ####
message(Sys.time(), " Pseudobulk")
sce_pseudo <- registration_pseudobulk(sce,
                                      var_registration = "cellType",
                                      var_sample_id = "individualID",
                                      covars = NULL
)

message("\nSCE Pseudobulk Dimensions:")
dim(sce_pseudo)

## Save results
saveRDS(sce_pseudo,
        file = here(
            "processed-data", "18_PsychENCODE_NMF", "pseudobulk",
            paste0("pseudobulk_", dataset, ".rds")
        )
)

# UCLA-ASD
dataset <- "UCLA-ASD"
file_name = "UCLA-ASD_annotated_mismatches_removed.h5ad"
file <- paste0(file_path, file_name)
sce <- readH5AD(file)

colData(sce)$cellType <- as.factor(make.names(colData(sce)$subclass))
table(sce$cellType)

rowData(sce)$gene_name <- rownames(sce) # save gene name as column of rowData
rownames(sce) <- rowData(sce)$featureid # have to make row names of object the ensembl id instead of gene names

## Logcounts
# default “X” contain the log-normalized counts
message(Sys.time(), " revert to counts")

## check for all 0s (just first 100 cols for mem)
stopifnot(any(assays(sce)$X[, 1:100] != 0))

counts(sce) <- assays(sce)$X # change to normalized counts
counts(sce)@x <- 2^(counts(sce)@x) - 1 ## remove log2(counts + 1)

data_sub <- data[data$Cohort == "UCLA-ASD" & data$Disorder == "control", c("Individual_ID", "Cohort", "Disorder")]

# only keep cells that have colData(sce)$individualID in data_sub$Individual_ID
sce <- sce[, colData(sce)$individualID %in% data_sub$Individual_ID]

#### Pseudobulk ####
message(Sys.time(), " Pseudobulk")
sce_pseudo <- registration_pseudobulk(sce,
                                      var_registration = "cellType",
                                      var_sample_id = "individualID",
                                      covars = NULL
)

message("\nSCE Pseudobulk Dimensions:")
dim(sce_pseudo)

## Save results
saveRDS(sce_pseudo,
        file = here(
            "processed-data", "18_PsychENCODE_NMF", "pseudobulk",
            paste0("pseudobulk_", dataset, ".rds")
        )
)

