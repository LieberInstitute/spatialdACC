library("spatialLIBD")
library("markdown") ## due to a shinyapps.io bug

## spatialLIBD uses golem
options("golem.app.prod" = TRUE)

## You need this to enable shinyapps to install Bioconductor packages
options(repos = BiocManager::repositories())

## Load the data (all paths are relative to this script's location)
load("spe_nnSVG_PRECAST_9_labels.Rdata", verbose = TRUE)
load("modeling-nnSVG_PRECAST_captureArea_9.Rdata", verbose = TRUE)
load("nnSVG_PRECAST_captureArea_9.Rdata", verbose = TRUE)

sig_genes <- readRDS("nnSVG_PRECAST_captureArea_9_sig_genes_all.rds")

# Added biocGenerics specifically for finding the right colnames, also due to R update
vars <- BiocGenerics::colnames(colData(spe))

## Deploy the website
spatialLIBD::run_app(
    spe,
    sce_layer = spe_pseudo,
    modeling_results = modeling_results,
    sig_genes = sig_genes,
    title = "spatialdACC, Visium",
    spe_discrete_vars = c(vars[grep("10x_", vars)], "ManualAnnotation", "layer"),
    spe_continuous_vars = c("sum_umi", "sum_gene", "expr_chrM", "expr_chrM_ratio"),
    default_cluster = "layer",
    docs_path = "www"
)
