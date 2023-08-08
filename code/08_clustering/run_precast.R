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

spe$spot_id = spe$key

spe$row <- spe$array_row
spe$col <- spe$array_col

seu <- CreateSeuratObject(
    counts=as.matrix(counts(spe)),
    meta.data=data.frame(colData(spe)),
    project="dACC")

seuList = list()

brains = unique(spe$brnum)

for (i in seq_along(brains)){
    seuList[[i]] = subset(x=seu, subset = brnum == brains[i])
}

preobj = CreatePRECASTObject(seuList = seuList, gene.number=2000, selectGenesMethod='HVGs')
preobj@seulist

PRECASTObj <- AddAdjList(preobj, platform = "Visium")

PRECASTObj <- AddParSetting(PRECASTObj, Sigma_equal = FALSE, coreNum = 8, maxIter = 30, verbose = TRUE)

K <- as.numeric(Sys.getenv("SGE_TASK_ID"))

tic()
PRECASTObj <- PRECAST(PRECASTObj, K = K)
toc()

save(PRECASTObj, file = here("processed-data", "08_clustering", "PRECAST", paste0("PRECASTObj_",K,".Rdata")))
