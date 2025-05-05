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

<img src="https://github.com/LieberInstitute/spatialdACC/blob/main/plots/dACC_Fig1_draft2.png" width="1000px" align="left" />

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
  - 👀 <https://libd.shinyapps.io/pseudobulk_HPC/>
    - Provides tools for visualization of pseudobulked Visium data.
  - 🔍 [dACC Samui
    browser](https://samuibrowser.com/from?url=data.libd.org/samuibrowser/&s=Br3942&s=Br8325&s=Br2720&s=Br2743&s=Br3942-VSPG&s=Br6423&s=Br6432&s=Br6471&s=Br6522&s=Br8325-VSPG&s=Br8492&s=Br8667)
    - Provides interactive spot-level visualization of Visium data.
- snRNA-seq (n = 10)
  - 👀 <https://libd.shinyapps.io/HPC_snRNAseq_data/>
    - Provides tools for visualization of snRNA-seq data.
   
## Data Access

All data, including raw FASTQ files and `SpaceRanger`/`CellRanger`
processed data outputs, can be accessed via Gene Expression Omnibus
(GEO) under accessions
[GSE264692](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE264692)
(SRT) and
[GSE264624](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE264624)
(snRNA-seq).

## Contact

We value public questions, as they allow other users to learn from the
answers. If you have any questions, please ask them at
[LieberInstitute/spatialdACC/issues](https://github.com/LieberInstitute/spatialdACC/issues)
and refrain from emailing us. Thank you again for your interest in our
work!

## Internal

JHPCE location:`/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/`
