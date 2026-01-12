setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')
suppressPackageStartupMessages({
    library("dplyr")
    library("purrr")
    library("Seurat")
    library("here")
    library("sessioninfo")
    library("SpatialExperiment")
    library("dplyr")
    library("purrr")
    library("tidyverse")
    library("spatialLIBD")
    library("PRECAST")
    library("tictoc")
    library("aricode")
})

#start from spe without WM-CC
load(here("processed-data", "08_clustering", "PRECAST", "spe_nnSVG_PRECAST_9_labels.Rdata"))

#load nnSVG results
load(file=here::here('processed-data', '08_clustering', 'nnSVG','nnSVG_1000.rda'))
#use gene ids instead of gene names
genes <- rownames(df_summaryReplicated)

colnames(spe) <- spe$key

leave_out <- unique(spe$sample_id)[as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))]

spe <- spe[, spe$sample_id != leave_out]

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

tic()
PRECASTObj <- PRECAST(PRECASTObj, K = 8)
toc()

save(PRECASTObj, file = here("processed-data", "08_clustering", "PRECAST", paste0("nnSVG_PRECASTObj_8_no_",leave_out,".Rdata")))

PRECASTObj <- SelectModel(PRECASTObj)
seuInt <- IntegrateSpaData(PRECASTObj, species = "Human")

# Merge with spe object
cluster_df <- seuInt@meta.data |>
    mutate(cluster = factor(cluster)) |>
    rename_with(~ paste0("PRECAST_", .x)) |>
    rownames_to_column(var = "key")

col_data_df <- colData(spe) |>
    data.frame() |>
    left_join(cluster_df, by="key")

rownames(col_data_df) <- colnames(spe)
colData(spe)$PRECAST_cluster <- col_data_df$PRECAST_cluster

# load original clusters, remove same sample, calculate ARI/NMI

spe_curr <- spe
load(here("processed-data", "08_clustering", "PRECAST", "spe_nnSVG_PRECAST_9_labels.Rdata"))
spe_orig <- spe

spe_orig <- spe_orig[, spe_orig$sample_id != leave_out]

print(dim(spe_curr))
print(dim(spe_orig))

na_idx <- is.na(spe_curr$PRECAST_cluster)

nmi_val <- NMI(spe_curr$PRECAST_cluster[!na_idx], spe_orig$nnSVG_PRECAST_captureArea_9[!na_idx])
ari_val <- ARI(spe_curr$PRECAST_cluster[!na_idx], spe_orig$nnSVG_PRECAST_captureArea_9[!na_idx])

nmi_val
ari_val

NewData <- data.frame(X1 = leave_out, X2 = nmi_val, X3 = ari_val)
write.table(NewData, file = here("processed-data", "08_clustering", "PRECAST", "nnSVG_PRECAST_8_crossval.csv"),
            append = TRUE, sep = ",",
            col.names = FALSE, row.names = FALSE)

