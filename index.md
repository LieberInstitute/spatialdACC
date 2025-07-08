# Spatio-molecular gene expression reflects dorsal anterior cingulate cortex structure and function in the human brain
## Overview

Welcome to the `spatial_dACC` project! In this study, we generated
spatially-resolved transcriptomics (SRT) and single-nucleus
RNA-sequencing (snRNA-seq) data from adjacent tissue sections of the
dorsal anterior cingulate cortex across ten adult neurotypical donors. SRT
data was generated using [10x Genomics
**Visium**](https://www.10xgenomics.com/products/spatial-gene-expression)
(n=17 capture areas). snRNA-seq data was generated using [10x Genomics
**Chromium**](https://www.10xgenomics.com/products/single-cell-gene-expression)
(n=10 total snRNA-seq libraries).

Thank you for your interest in our work!

## Study design

![Experimental Overview](./img/dACC_overview.png)

Experimental design to generate paired single-nucleus RNA-sequencing
(snRNA-seq) and spatially-resolved transcriptomics (SRT) data in the
human dorsal anterior cingulate cortex. (A) dACC (blue) and dlPFC (pink) outlined at midsagittal (dACC) and lateral (dlPFC) levels as well as in a coronal hemislabs (top). dACC tissue (blue) used here was sourced from the same n=10 neurotypical control donors previously similarly profiled for dlPFC (pink) to facilitate within-donor comparisons of dACC agranular cortex vs. dlPFC granular cortex (bottom).

(B) Fresh-frozen coronal brain slab taken at the level of the anterior striatum overlaid with the outline of major landmarks from the Atlas of the Human Brain (Mai et al. 2015). Three Brodmann areas (BA) 33, 24, and 32 corresponding to the dACC are highlighted in yellow (top left). H&E staining of tissue cryosections confirmed inclusion of dACC on the tissue block (top right). Following anatomical validation, cryosections were collected for 10X Visium and Chromium assays from the same brain block for each donor (bottom). (This figure was created with
[Biorender](https://biorender.com))

## Interactive Websites

All of these interactive websites are powered by open source software,
namely:

- 🔍 [`samui`](http://dx.doi.org/10.1017/S2633903X2300017X)
- 👀 [`iSEE`](https://doi.org/10.12688%2Ff1000research.14966.1)

We provide the following interactive websites, organized by dataset with
software labeled by emojis:

- Visium (n = 17)
  - 👀 <https://libd.shinyapps.io/spatialdACC_Visium/>
    - Provides tools for visualization of pseudobulked Visium data.
  - 🔍 [dACC Samui
    browser](https://samuibrowser.com/from?url=data.libd.org/samuibrowser/&s=Br2720_V12J03-002_A1&s=Br2720_V12N28-331_A1&s=Br2720_V12N28-332_A1&s=Br2743_V12N28-334_C1&s=Br3942_V12N28-334_A1&s=Br6423_V12N28-334_D1&s=Br6432_V12J03-002_B1&s=Br6432_V12N28-331_B1&s=Br6432_V12N28-332_B1&s=Br6471_V12J03-002_C1&s=Br6471_V12N28-331_C1&s=Br6471_V12N28-332_C1&s=Br6522_V12N28-331_D1&s=Br6522_V12N28-332_D1&s=Br8325_V12Y31-080_C1&s=Br8492_V12N28-334_B1&s=Br8667_V12Y31-080_B1)
    - Provides interactive spot-level visualization of Visium data.
- snRNA-seq (n = 10)
  - 👀 <https://libd.shinyapps.io/snRNAseq_dACC/>
    - Provides tools for visualization of snRNA-seq data.
   
## Data Access

All data, including raw FASTQ files and `SpaceRanger`/`CellRanger`
processed data outputs, can be accessed via Gene Expression Omnibus
(GEO) under accessions
[GSE296731](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE296731)
(SRT) and
[GSE296789](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE296789)
(snRNA-seq).

- R objects used for analysis are hosted at a public GLOBUS endpoint [https://research.libd.org/globus/](https://research.libd.org/globus/)
- Zenodo Archive for this project can be found at [https://doi.org/10.5281/zenodo.15830481](https://doi.org/10.5281/zenodo.15830481)

## Contact

We value public questions, as they allow other users to learn from the
answers. If you have any questions, please ask them at
[LieberInstitute/spatialdACC/issues](https://github.com/LieberInstitute/spatialdACC/issues)
and refrain from emailing us. Thank you again for your interest in our
work!



