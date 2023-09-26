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

colnames(spe) <- spe$key

seuList <- unique(spe$sample_id) |>
    set_names(unique(spe$sample_id)) |>
    map(.f = function(id) {
        tmp_spe <- spe[, spe$sample_id == id]

        tmp_spe$row <- tmp_spe$array_row
        tmp_spe$col <- tmp_spe$array_col

        # browser()
        CreateSeuratObject(
            counts=as.matrix(counts(tmp_spe)),
            meta.data=data.frame(colData(tmp_spe)),
            project="dACC")
    })

set.seed(1)
preobj <- CreatePRECASTObject(seuList = seuList, gene.number=2000, selectGenesMethod='HVGs',
                              premin.spots = 1, premin.features=1, postmin.spots=1, postmin.features=1)
preobj@seulist

PRECASTObj <- AddAdjList(preobj, platform = "Visium")

PRECASTObj <- AddParSetting(PRECASTObj, Sigma_equal = FALSE,  maxIter = 30, verbose = TRUE)

tic()
PRECASTObj <- PRECAST(PRECASTObj, K = 5:20)
toc()

resList <- PRECASTObj@resList
PRECASTObj <- SelectModel(PRECASTObj)

str(PRECASTObj@resList)

true_cluster <- lapply(PRECASTObj@seulist, function(x) x$true_cluster)
str(true_cluster)
mclust::adjustedRandIndex(unlist(PRECASTObj@resList$cluster), unlist(true_cluster))

save(PRECASTObj, file = here("processed-data", "08_clustering", "PRECAST", "PRECASTObj_5_20.Rdata"))
