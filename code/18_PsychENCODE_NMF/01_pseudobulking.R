library(SingleCellExperiment)
library(spatialLIBD)
library(here)
library(zellkonverter)

# Read in meta data
url <- "https://brainscope.gersteinlab.org/data/sample_metadata/PEC2_sample_metadata.txt"
data <- read.delim(url, header = TRUE)

# Paths to h5ad files
file_path <- "/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/psychENCODE/version5/"

# List of dataset files
data_list <- c("CMC_annotated.h5ad", "DevBrain-snRNAseq_annotated.h5ad",
               "IsoHuB_annotated.h5ad",
               "MultiomeBrain-DLPFC_annotated.h5ad", "PTSDBrainomics_annotated.h5ad",
               "SZBDMulti-Seq_annotated.h5ad",
               "UCLA-ASD_annotated_mismatches_removed.h5ad")

# Process each dataset
for (dataset_file in data_list) {

    # Extract dataset name
    dataset <- gsub("_annotated.h5ad|_annotated_mismatches_removed.h5ad", "", dataset_file)

    # Construct file path for the h5ad file
    file <- paste0(file_path, dataset_file)

    # Read the H5AD file
    sce <- readH5AD(file)

    # Update colData and rowData
    colData(sce)$cellType <- as.factor(make.names(colData(sce)$subclass))
    table(sce$cellType)

    # Load raw counts matrix
    raw_csv_path <- here("processed-data", "18_PsychENCODE_NMF", paste0("raw_", dataset, ".csv"))
    counts <- read.csv(raw_csv_path)

    counts <- t(counts)
    rownames(counts) <- rownames(assay(sce, "X"))
    colnames(counts) <- colnames(assay(sce, "X"))

    assays(sce)$counts <- counts

    saveRDS(sce,
            file = here(
                "processed-data", "18_PsychENCODE_NMF", "raw",
                paste0(dataset, ".rds")
            )
    )

    # Subset data
    # the Cohort variable does not exactly match the dataset name
    
    # [1] "CMC"                 "DevBrain"            "Girgenti-snMultiome"
    # [4] "IsoHuB"              "LIBD"                "Ma_et_al"           
    # [7] "MultiomeBrain"       "PTSDBrainomics"      "ROSMAP"             
    # [10] "SZBDMulti-Seq"       "UCLA-ASD"            "Velmeshev_et_al"
    
    if (dataset == "CMC") {
        cohort <- "CMC"
    }
    if(dataset == "DevBrain-snRNAseq") {
        cohort <- "DevBrain"
    }
    if(dataset == "IsoHuB") {
        cohort <- "IsoHuB"
    }
    if(dataset == "MultiomeBrain-DLPFC") {
        cohort <- "MultiomeBrain"
    }
    if(dataset == "PTSDBrainomics") {
        cohort <- "PTSDBrainomics"
    }
    if(dataset == "SZBDMulti-Seq") {
        cohort <- "SZBDMulti-Seq"
    }
    if(dataset == "UCLA-ASD") {
        cohort <- "UCLA-ASD"
    }
    
    data_sub <- data[data$Cohort == cohort & data$Disorder == "control", c("Individual_ID", "Cohort", "Disorder")]

    # Filter SCE based on Individual_ID
    sce <- sce[, colData(sce)$individualID %in% data_sub$Individual_ID]

    # Pseudobulk
    message(Sys.time(), " Pseudobulk for dataset: ", dataset)
    sce_pseudo <- registration_pseudobulk(sce,
                                          var_registration = "cellType",
                                          var_sample_id = "individualID",
                                          covars = NULL
    )

    message("\nSCE Pseudobulk Dimensions for dataset ", dataset, ":")
    dim(sce_pseudo)

    # Save results
    saveRDS(sce_pseudo,
            file = here(
                "processed-data", "18_PsychENCODE_NMF", "pseudobulk",
                paste0("pseudobulk_", dataset, ".rds")
            )
    )
}
