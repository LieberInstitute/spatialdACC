# qrsh -l bluejay,mem_free=80G,h_vmem=80G
# cd /dcs04/lieber/marmaypag/BLA_molprofile_LIBD1070/BLA_molprofile/
setwd("/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/")

library("SingleCellExperiment")
library("DropletUtils")
library("here")
library("rtracklayer")
#library("lobstr")
library("sessioninfo")
library("dplyr")

# Read in the CSV file as a data frame
tmp <- read.delim(here("raw-data","sample_info",
                "running-Chromium_summary.csv"),
                header = T,sep=',')

# View the subsetted data
#print(tmp)


##set up sample data table
sample_info <- data.frame(
  sample_id = tmp$Sample,
  sort = tmp$Tissue,
  brain = tmp$Brain,
  round = tmp$Round
)
sample_info<-sample_info[-c(1),]
stopifnot(all(!duplicated(sample_info$Sample)))

# add path to cellranger output
#sample_info$sample_path<-rep(NA,5)
sample_info$sample_path<- file.path(
  here::here("processed-data","03_cellranger"),
  sample_info$sample_id,
  "outs",
  "raw_feature_bc_matrix"
)

## Get the rest of the donor info at some point
#
# XXXXXXXXXX
# XXXXXXXXXX
# XXXXXXXXXX
#

## Build basic SCE
message("Read 10x data and create sce - ", Sys.time())
# Read 10x data and create sce - 2023-12-12 12:22:37.087479

sce <- read10xCounts(
  sample_info$sample_path,
  sample_info$sample_id,
  type = "sparse",
  col.names = TRUE
)
message("RDone - ", Sys.time())
# RDone - 2023-12-12 12:29:33.440454

# Note that the rownames are the Ensembl gene names - let's switch to the more familiar gene symbols:
# ...but also notice that length(unique(rowData(sce)$Symbol)) != nrow(sce)
#   - That's because some gene symbols are used multiple times.

## Use code from https://github.com/LieberInstitute/Visium_IF_AD/commit/08df3f7e4a3178563d6b4b1861b664b21466b395#diff-10cb35de98e2a3e5f4235cd88f6dabce5469eead2b2db1fd7121126849fcf585L100
## Read in the gene information from the annotation GTF file
gtf <-
    rtracklayer::import(
        "/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-GRCh38-2020-A/genes/genes.gtf"
    )
gtf <- gtf[gtf$type == "gene"]
names(gtf) <- gtf$gene_id

## Match the genes
match_genes <- match(rownames(sce), gtf$gene_id)
stopifnot(all(!is.na(match_genes)))

## Keep only some columns from the gtf
mcols(gtf) <- mcols(gtf)[, c("source", "type", "gene_id", "gene_version", "gene_name", "gene_type")]

## Add the gene info to our SPE object
rowRanges(sce) <- gtf[match_genes]

# -> To avoid getting into sticky situations, let's 'uniquify' these names:
rowData(sce)$Symbol.uniq <- scuttle::uniquifyFeatureNames(rowData(sce)$gene_id, rowData(sce)$gene_name)
rownames(sce) <- rowData(sce)$Symbol.uniq


# Add metadata
# sce$key <- paste0(sce$Barcode, "_", sce$Sample)
# new_col <- merge(colData(sce), sample_info[, -which(colnames(sample_info) == "sample_path")])
# new_col1 <- new_col[match(sce$key, new_col$key), ]
sce$key <- colnames(sce)
sample_info = sample_info |>
  rename(Sample = sample_id) |>
  select(-sample_path)

colData(sce) = colData(sce) |>
  as_tibble() |>
  left_join(sample_info, by = 'Sample') |>
  DataFrame()

colnames(sce) = sce$key
rownames(colData(sce)) = colnames(sce)


#stopifnot(identical(sce$key, new_col$key))
#rownames(new_col) <- colnames(sce)
#colData(sce) <- new_col

## Inspect object
sce
# class: SingleCellExperiment 
# dim: 36601 19473661 
# metadata(1): Samples
# assays(1): counts
# rownames(36601): MIR1302-2HG FAM138A ... AC007325.4 AC007325.2
# rowData names(7): source type ... gene_type Symbol.uniq
# colnames(19473661): 1_AAACCCAAGAAACCCA-1 1_AAACCCAAGAAACTCA-1 ...
#   10_TTTGTTGTCTTTGCTG-1 10_TTTGTTGTCTTTGGAG-1
# colData names(6): Sample Barcode ... brain round
# reducedDimNames(0):
# mainExpName: NULL
# altExpNames(0):


if (!dir.exists(here("processed-data", "04_build_sce"))) dir.create(here("processed-data", "04_build_sce"))
save(sce, file = here("processed-data", "04_build_sce", "1c-10c_sce_raw.rda"))


## Size in Gb
#lobstr::obj_size(sce) / 1024^3
# 7.43

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

# ─ Session info ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  setting  value
#  version  R version 4.3.2 Patched (2023-11-13 r85524)
#  os       Rocky Linux 9.2 (Blue Onyx)
#  system   x86_64, linux-gnu
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       US/Eastern
#  date     2023-12-12
#  pandoc   3.1.3 @ /jhpce/shared/community/core/conda_R/4.3.x/bin/pandoc

# ─ Packages ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package              * version   date (UTC) lib source
#  abind                  1.4-5     2016-07-21 [2] CRAN (R 4.3.2)
#  beachmat               2.18.0    2023-10-24 [2] Bioconductor
#  Biobase              * 2.62.0    2023-10-24 [2] Bioconductor
#  BiocGenerics         * 0.48.1    2023-11-01 [2] Bioconductor
#  BiocIO                 1.12.0    2023-10-24 [2] Bioconductor
#  BiocParallel           1.36.0    2023-10-24 [2] Bioconductor
#  Biostrings             2.70.1    2023-10-25 [2] Bioconductor
#  bitops                 1.0-7     2021-04-24 [2] CRAN (R 4.3.2)
#  cli                    3.6.1     2023-03-23 [2] CRAN (R 4.3.2)
#  codetools              0.2-19    2023-02-01 [3] CRAN (R 4.3.2)
#  crayon                 1.5.2     2022-09-29 [2] CRAN (R 4.3.2)
#  DelayedArray           0.28.0    2023-10-24 [2] Bioconductor
#  DelayedMatrixStats     1.24.0    2023-10-24 [2] Bioconductor
#  dplyr                * 1.1.3     2023-09-03 [2] CRAN (R 4.3.2)
#  dqrng                  0.3.1     2023-08-30 [2] CRAN (R 4.3.2)
#  DropletUtils         * 1.22.0    2023-10-24 [2] Bioconductor
#  edgeR                  4.0.1     2023-10-29 [2] Bioconductor
#  fansi                  1.0.5     2023-10-08 [2] CRAN (R 4.3.2)
#  generics               0.1.3     2022-07-05 [2] CRAN (R 4.3.2)
#  GenomeInfoDb         * 1.38.1    2023-11-08 [2] Bioconductor
#  GenomeInfoDbData       1.2.11    2023-11-15 [2] Bioconductor
#  GenomicAlignments      1.38.0    2023-10-24 [2] Bioconductor
#  GenomicRanges        * 1.54.1    2023-10-29 [2] Bioconductor
#  glue                   1.6.2     2022-02-24 [2] CRAN (R 4.3.2)
#  HDF5Array              1.30.0    2023-10-24 [2] Bioconductor
#  here                 * 1.0.1     2020-12-13 [2] CRAN (R 4.3.2)
#  IRanges              * 2.36.0    2023-10-24 [2] Bioconductor
#  lattice                0.22-5    2023-10-24 [3] CRAN (R 4.3.2)
#  lifecycle              1.0.4     2023-11-07 [2] CRAN (R 4.3.2)
#  limma                  3.58.1    2023-10-31 [2] Bioconductor
#  lobstr                 1.1.2     2022-06-22 [2] CRAN (R 4.3.2)
#  locfit                 1.5-9.8   2023-06-11 [2] CRAN (R 4.3.2)
#  magrittr               2.0.3     2022-03-30 [2] CRAN (R 4.3.2)
#  Matrix                 1.6-3     2023-11-14 [3] CRAN (R 4.3.2)
#  MatrixGenerics       * 1.14.0    2023-10-24 [2] Bioconductor
#  matrixStats          * 1.1.0     2023-11-07 [2] CRAN (R 4.3.2)
#  pillar                 1.9.0     2023-03-22 [2] CRAN (R 4.3.2)
#  pkgconfig              2.0.3     2019-09-22 [2] CRAN (R 4.3.2)
#  prettyunits            1.2.0     2023-09-24 [2] CRAN (R 4.3.2)
#  R.methodsS3            1.8.2     2022-06-13 [2] CRAN (R 4.3.2)
#  R.oo                   1.25.0    2022-06-12 [2] CRAN (R 4.3.2)
#  R.utils                2.12.2    2022-11-11 [2] CRAN (R 4.3.2)
#  R6                     2.5.1     2021-08-19 [2] CRAN (R 4.3.2)
#  Rcpp                   1.0.11    2023-07-06 [2] CRAN (R 4.3.2)
#  RCurl                  1.98-1.13 2023-11-02 [2] CRAN (R 4.3.2)
#  restfulr               0.0.15    2022-06-16 [2] CRAN (R 4.3.2)
#  rhdf5                  2.46.0    2023-10-24 [2] Bioconductor
#  rhdf5filters           1.14.1    2023-11-06 [2] Bioconductor
#  Rhdf5lib               1.24.0    2023-10-24 [2] Bioconductor
#  rjson                  0.2.21    2022-01-09 [2] CRAN (R 4.3.2)
#  rlang                  1.1.2     2023-11-04 [2] CRAN (R 4.3.2)
#  rprojroot              2.0.4     2023-11-05 [2] CRAN (R 4.3.2)
#  Rsamtools              2.18.0    2023-10-24 [2] Bioconductor
#  rtracklayer          * 1.62.0    2023-10-24 [2] Bioconductor
#  S4Arrays               1.2.0     2023-10-24 [2] Bioconductor
#  S4Vectors            * 0.40.1    2023-10-26 [2] Bioconductor
#  scuttle                1.12.0    2023-10-24 [2] Bioconductor
#  sessioninfo          * 1.2.2     2021-12-06 [2] CRAN (R 4.3.2)
#  SingleCellExperiment * 1.24.0    2023-10-24 [2] Bioconductor
#  SparseArray            1.2.2     2023-11-07 [2] Bioconductor
#  sparseMatrixStats      1.14.0    2023-10-24 [2] Bioconductor
#  statmod                1.5.0     2023-01-06 [2] CRAN (R 4.3.2)
#  SummarizedExperiment * 1.32.0    2023-10-24 [2] Bioconductor
#  tibble                 3.2.1     2023-03-20 [2] CRAN (R 4.3.2)
#  tidyselect             1.2.0     2022-10-10 [2] CRAN (R 4.3.2)
#  utf8                   1.2.4     2023-10-22 [2] CRAN (R 4.3.2)
#  vctrs                  0.6.4     2023-10-12 [2] CRAN (R 4.3.2)
#  withr                  2.5.2     2023-10-30 [2] CRAN (R 4.3.2)
#  XML                    3.99-0.15 2023-11-02 [2] CRAN (R 4.3.2)
#  XVector                0.42.0    2023-10-24 [2] Bioconductor
#  yaml                   2.3.7     2023-01-23 [2] CRAN (R 4.3.2)
#  zlibbioc               1.48.0    2023-10-24 [2] Bioconductor

#  [1] /users/rmiller/R/4.3.x
#  [2] /jhpce/shared/community/core/conda_R/4.3.x/R/lib64/R/site-library
#  [3] /jhpce/shared/community/core/conda_R/4.3.x/R/lib64/R/library

# ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
