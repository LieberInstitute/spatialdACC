library("SingleCellExperiment")
library("iSEE")
library("shiny")
library("paletteer")

load("sce.Rdata", verbose = TRUE)
#load("sce_for_iSEE_LS.rda", verbose = TRUE)

#stopifnot(all(unique(sce$cellType.final) %in% names(cell_cols.clean)))

## Don't run this on app.R since we don't want to run this every single time
# lobstr::obj_size(sce)˚
# 876.34 MB

source("initial.R", print.eval = TRUE)

## From https://github.com/LieberInstitute/10xPilot_snRNAseq-human/blob/810b47364af4c8afe426bd2a6b559bd6a9f1cc98/shiny_apps/tran2021_AMY/app.R#L10-L14
## Related to https://github.com/iSEE/iSEE/issues/568
colData(sce) <- cbind(
    colData(sce)[, !colnames(colData(sce)) %in% c("Sample", "round")],
    colData(sce)[, c("round", "Sample")]
)

sce$Sample <- as.factor(sce$Sample)

#sce <- registerAppOptions(sce, color.maxlevels = length(cell_cols.clean))
iSEE(
    sce,
    appTitle = "snRNAseq_dACC",
    initial = initial,
    colormap = ExperimentColorMap(colData = list(
        Sample = function(n) {
            cols <- paletteer::paletteer_d(
                palette = "RColorBrewer::Dark2",
                n = length(unique(sce$Sample))
            )
            cols <- as.vector(cols)
            names(cols) <- levels(sce$Sample)
            return(cols)
        }
        # ,
        # cellType.final = function(n) {
        #     return(cell_cols.clean)
        # }
    ))
)
