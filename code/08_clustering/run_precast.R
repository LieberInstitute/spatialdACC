setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')
suppressPackageStartupMessages({
    library("dplyr")
    library("purrr")
    library("Seurat")
    library("here")
    library("sessioninfo")
    library("SpatialExperiment")
    library("PRECAST")
    library("tictoc")
})

#start from spe without batch correction
load(here("processed-data", "06_preprocessing", "spe_dimred.Rdata"))

spe_to_seuratList <- function(spe){
    uniq_sample_id <- colData(spe)$sample_id |> unique()

    # Create a seurat object for each unique sample_id
    map(uniq_sample_id,
        .f = function(smpl_id, spe){
            ret_spe <- spe[, colData(spe)$sample_id == smpl_id]
            ret_seurat <- spe_to_seurat(ret_spe)

            return(ret_seurat)
        },
        spe = spe)
}

spe_to_seurat <- function(spe){

    ret <- CreateSeuratObject(
        counts=assays(spe)$counts,
        meta.data=data.frame(
            row=spatialCoords(spe)[,1],
            col=spatialCoords(spe)[,2])
    )

    return(ret)
}

seuList <- spe_to_seuratList(spe)

preobj = CreatePRECASTObject(seuList = seuList, gene.number=2000, selectGenesMethod='HVGs')
preobj@seulist

PRECASTObj <- AddAdjList(preobj, platform = "Visium")

PRECASTObj <- AddParSetting(PRECASTObj, Sigma_equal = FALSE,  maxIter = 30, verbose = TRUE)

K <- as.numeric(Sys.getenv("SGE_TASK_ID"))

tic()
set.seed(1)
PRECASTObj <- PRECAST(PRECASTObj, K = K)
toc()

save(PRECASTObj, file = here("processed-data", "08_clustering", "PRECAST", paste0("PRECASTObj_",K,".Rdata")))
