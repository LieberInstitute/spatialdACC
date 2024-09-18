library("SingleCellExperiment")
library("Seurat")
library("Azimuth")
library("SeuratData")
library("patchwork")
library("here")
library("sessioninfo")

plot_dir <- here("plots", "snRNA-seq", "05_azimuth")
if (!dir.exists(plot_dir)) dir.create(plot_dir)

load(file = here::here("processed-data", "snRNA-seq", "03_batch_correction", paste0("sce_harmony.Rdata")))
query <- CreateSeuratObject(
    counts = as.matrix(counts(sce)),
    meta.data = data.frame(colData(sce)),
    project = "dACC"
)

## Use Azimuth Refrence ##
set.seed(20)
query <- RunAzimuth(query, reference = "humancortexref") ## Cell annotation with Azimuth

sce$cellType_azimuth <- query$predicted.subclass
table(query$predicted.subclass)

#      Astro       Endo    L2/3 IT      L5 ET      L5 IT    L5/6 NP      L6 CT
#      3009         89       3833        164       1657        387       1031
#     L6 IT L6 IT Car3        L6b      Lamp5  Micro-PVM      Oligo        OPC
#       960        138        577        949       2616      13256       2085
#     Pvalb       Sncg        Sst  Sst Chodl        Vip       VLMC
#      1610        427       1129         34       1051        159



# replace "-" with nothing
colData(sce)$cellType_azimuth <- replace(colData(sce)$cellType_azimuth, colData(sce)$cellType_azimuth == "Micro-PVM", "MicroPVM")

# replace spaces with underscores
colData(sce)$cellType_azimuth <- replace(colData(sce)$cellType_azimuth, colData(sce)$cellType_azimuth == "L2/3 IT", "L2_3_IT")
colData(sce)$cellType_azimuth <- replace(colData(sce)$cellType_azimuth, colData(sce)$cellType_azimuth == "L5 ET", "L5_ET")
colData(sce)$cellType_azimuth <- replace(colData(sce)$cellType_azimuth, colData(sce)$cellType_azimuth == "L5 IT", "L5_IT")
colData(sce)$cellType_azimuth <- replace(colData(sce)$cellType_azimuth, colData(sce)$cellType_azimuth == "L5/6 NP", "L5_6_NP")
colData(sce)$cellType_azimuth <- replace(colData(sce)$cellType_azimuth, colData(sce)$cellType_azimuth == "L6 CT", "L6_CT")
colData(sce)$cellType_azimuth <- replace(colData(sce)$cellType_azimuth, colData(sce)$cellType_azimuth == "L6 IT", "L6_IT")
colData(sce)$cellType_azimuth <- replace(colData(sce)$cellType_azimuth, colData(sce)$cellType_azimuth == "L6 IT Car3", "L6_IT_Car3")

save(sce, file = here("processed-data", "snRNA-seq", "05_azimuth", "sce_azimuth.Rdata"))

