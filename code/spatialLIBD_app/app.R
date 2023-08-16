library("spatialLIBD")
library("markdown") ## due to a shinyapps.io bug

## spatialLIBD uses golem
options("golem.app.prod" = TRUE)

## You need this to enable shinyapps to install Bioconductor packages
options(repos = BiocManager::repositories())

## Load the data (all paths are relative to this script's location)
load("spe_raw_if.Rdata", verbose = TRUE)
# spe$CellCount <- spe$segmentation_info
# spe$CellCount <- spe$NBW
# vars <- colnames(colData(spe))

# Added biocGenerics specifically for finding the right colnames, also due to R update
vars <- BiocGenerics::colnames(colData(spe))

## Deploy the website
spatialLIBD::run_app(
    spe,
    sce_layer = NULL,
    modeling_results = NULL,
    sig_genes = NULL,
    title = "spatialdACC, Visium IF",
    spe_discrete_vars = c(vars[grep("10x_", vars)], "ManualAnnotation"),
    spe_continuous_vars = c("sum_umi", "sum_gene", "expr_chrM", "expr_chrM_ratio"),
    default_cluster = "10x_graphclust"
)
