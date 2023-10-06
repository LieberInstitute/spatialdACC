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

#load nnSVG results
load(file=here::here('processed-data', '08_clustering', 'nnSVG','nnSVG_1000.rda'))
#use gene ids instead of gene names
genes <- rownames(df_summaryReplicated)

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
preobj <- CreatePRECASTObject(seuList = seuList, selectGenesMethod=NULL,
                              customGenelist = genes,
                              premin.spots = 1, premin.features=1, postmin.spots=1, postmin.features=1)
preobj@seulist

PRECASTObj <- AddAdjList(preobj, platform = "Visium")

PRECASTObj <- AddParSetting(PRECASTObj, Sigma_equal = FALSE,  maxIter = 30, verbose = TRUE)

K <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

tic()
PRECASTObj <- PRECAST(PRECASTObj, K = K)
toc()

save(PRECASTObj, file = here("processed-data", "08_clustering", "PRECAST", paste0("nnSVG_PRECASTObj_",K,".Rdata")))
