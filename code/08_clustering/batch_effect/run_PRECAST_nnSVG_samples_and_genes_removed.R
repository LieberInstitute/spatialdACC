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

#remove samples with batch effect suspected
# "V12N28−332_A1" and "V12N28−332_B1"
suspected_batch_samples <- c("V12N28-332_A1", "V12N28-332_B1")
spe <- spe[, !colData(spe)$sample_id %in% suspected_batch_samples]

#load nnSVG results
load(file=here::here('processed-data', '08_clustering', 'batch_effect','nnSVG_1000_samples_removed.rda'))
#use gene ids instead of gene names
genes <- rownames(df_summaryReplicated)

biased_genes <- c("ENSG00000256618", "ENSG00000255823",
                  "ENSG00000198840",
                  "ENSG00000198886", "ENSG00000198899", "ENSG00000198804", "ENSG00000198938")

genes <- genes[!genes %in% biased_genes]

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

save(PRECASTObj, file = here("processed-data", "08_clustering", "batch_effect", paste0("nnSVG_PRECASTObj_samples_and_genes_removed",K,".Rdata")))
