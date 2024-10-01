setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')
suppressPackageStartupMessages(library("spatialLIBD"))
suppressPackageStartupMessages(library("ggspavis"))
suppressPackageStartupMessages(library("gridExtra"))
suppressPackageStartupMessages(library("here"))
library("nmfLabelTransfer")
library("scran")

load(here("processed-data", "VSPG", "01_QC", "spe_QC.Rdata"), verbose = TRUE)

# remove spots with low sum and low detected features
spe <- spe[,!colData(spe)$low_sum_br & !colData(spe)$low_detected_br]
dim(spe)
# [1] 28096 17278

spe_target <- spe
spe_target <- logNormCounts(spe_target)

# load spe for k=9 without WM-CC
load(
    file = here("processed-data", "08_clustering", "PRECAST", "spe_nnSVG_PRECAST_9_labels.Rdata")
)

spe_source <- spe
dim(spe)

layer_labs <- "layer"

res <- transfer_labels(targets=spe_target,
                       source=spe_source,
                       assay="logcounts",
                       annotationsName=layer_labs,
                       technicalVarName="sample_id",
                       save_nmf=T,
                       nmf_path = here("processed-data", "VSPG", "02_label_transfer", "spaTransfer_nmf_results.rds"),
                       k=NULL, tol=1e-5)

target_with_preds <- res$targets
print(target_with_preds)
# display the results
table(target_with_preds$nmf_preds)

# save the results
save(res, file = here("processed-data", "VSPG", "02_label_transfer", "spaTransfer_target_with_preds.rds"))
