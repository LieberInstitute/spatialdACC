
setwd("/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/")
suppressPackageStartupMessages({
    library("here")
    library("SpatialExperiment")
    library("spatialLIBD")
    library("rtracklayer")
    library("lobstr")
    library("sessioninfo")
})

## Define some info for the samples
load(here::here("code", "REDCap", "REDCap_dACC.rda"))

# sample_info <- data.frame(dateImg = as.Date(REDCap_dACC$date)) 
# sample_info$experimenterImg <- as.factor(REDCap_dACC$experimenter_img)
sample_info <- data.frame(slide = as.factor(REDCap_dACC$slide))
sample_info$array <- as.factor(REDCap_dACC$array)
sample_info$brnum <- as.factor(sapply(strsplit(REDCap_dACC$sample, "-"), `[`, 1))
sample_info$replicate <- as.factor(REDCap_dACC$serial)
sample_info$sample_id <- paste(sample_info$slide, sample_info$array, sep = "_")

sample_info$sample_path[sample_info$slide == "V12N28-331" | sample_info$slide == "V12J03-002"] <- file.path(here::here("processed-data", "01_spaceranger", "2023-03-29_Transfer10x_SPage"), sample_info$sample_id[sample_info$slide == "V12N28-331" | sample_info$slide == "V12J03-002"], "outs")
sample_info$sample_path[sample_info$sample_id == "V12N28-331_B1"] <- file.path(here::here("processed-data", "01_spaceranger", "2023-02-09_Transfer10x_SPage"), "V12N28-331_B1", "outs")
sample_info$sample_path[sample_info$sample_id == "V12J03-002_D1"] <- file.path(here::here("processed-data", "01_spaceranger", "2023-02-09_Transfer10x_SPage"), "V12J03-002_D1", "outs")

sample_info$sample_path[sample_info$slide == "V12N28-332"] <- file.path(here::here("processed-data", "01_spaceranger", "2023-02-09_Transfer10x_SPage"), sample_info$sample_id[sample_info$slide == "V12N28-332"], "outs")
sample_info$sample_path[sample_info$slide == "V12N28-334" | sample_info$slide == "V12Y31-080"] <- file.path(here::here("processed-data", "01_spaceranger", "spaceranger_2023-03-14_SPag021323"), sample_info$sample_id[sample_info$slide == "V12N28-334" | sample_info$slide == "V12Y31-080"], "outs")
sample_info$sample_path[sample_info$sample_id == "V12N28-334_A1"] <- file.path(here::here("processed-data", "01_spaceranger", "spaceranger_2023-04-03_SPag022423"), "V12N28-334_A1", "outs")

##discard barnyard samples
sample_info = sample_info[sample_info$sample_id != "V12Y31-080_D1" & sample_info$sample_id != "V12J03-002_D1",]

stopifnot(all(file.exists(sample_info$sample_path)))

## Define the donor info using information from
donor_info <- read.csv(file.path(here::here("raw-data", "sample_info", "demographicInfo_Geo.csv")), header = TRUE, stringsAsFactors = FALSE)

## check if all donor info included
setdiff(sample_info$brnum, donor_info$brnum)

## Combine sample info with the donor info
sample_info <- merge(sample_info, donor_info)

##re order to bring in-house samples 1st in the list
sample_info = sample_info[c(4:6,15:17,1:3,7:14),]

## Build basic SPE
Sys.time()
spe <- read10xVisiumWrapper(
    sample_info$sample_path,
    sample_info$sample_id,
    type = "sparse",
    data = "raw",
    images = c("lowres", "hires", "detected", "aligned"),
    load = TRUE
)
Sys.time()

## Add the study design info
add_design <- function(spe) {
    new_col <- merge(colData(spe), sample_info)
    ## Fix order
    new_col <- new_col[match(spe$key, new_col$key), ]
    stopifnot(identical(new_col$key, spe$key))
    rownames(new_col) <- rownames(colData(spe))
    colData(spe) <-
        new_col[, -which(colnames(new_col) == "sample_path")]
    return(spe)
}
spe <- add_design(spe)

# dir.create(here::here("processed-data", "pilot_data_checks"), showWarnings = FALSE)
save(spe, file = here::here("processed-data", "02_build_spe", "spe_raw.Rdata"))

## Size in Gb
lobstr::obj_size(spe_raw)
# 5.141452 B
dim(spe_raw)
# 36601 159744

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
