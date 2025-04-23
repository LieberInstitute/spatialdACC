library(ggplot2)
library(tidyverse)
library(patchwork)
library(preprocessCore)

setwd("../../../Downloads/Data")

# set up dataframes
filenames <- list.files(".", pattern="*.csv", full.names=TRUE)

Br6432_dACC <- read.csv(filenames[1])
Br6432_dlPFC <- read.csv(filenames[2])
Br8325_dlPFC <- read.csv(filenames[3])
Br8325_left_dACC <- read.csv(filenames[4])
Br8325_right_dACC <- read.csv(filenames[5])

cols <- c("X25xSil.Opal.520.Copies", "X20x.Opal.570.Copies",
          "X20x.Opal.620.Copies", "X20x.Opal.690.Copies")

Br6432_dACC_nucleus_area <- Br6432_dACC$Nucleus.Area..µm..
Br6432_dACC <- Br6432_dACC[,cols]
colnames(Br6432_dACC) <- c("PCP4", "POU3F1", "SULF2", "GABRQ")

Br6432_dlPFC_nucleus_area <- Br6432_dlPFC$Nucleus.Area..µm..
Br6432_dlPFC <- Br6432_dlPFC[,cols]
colnames(Br6432_dlPFC) <- c("PCP4", "POU3F1", "SULF2", "GABRQ")

Br8325_dlPFC_nucleus_area <- Br8325_dlPFC$Nucleus.Area..µm..
Br8325_dlPFC <- Br8325_dlPFC[,cols]
colnames(Br8325_dlPFC) <- c("PCP4", "POU3F1", "SULF2", "GABRQ")

Br8325_left_dACC_nucleus_area <- Br8325_left_dACC$Nucleus.Area..µm..
Br8325_left_dACC <- Br8325_left_dACC[,cols]
colnames(Br8325_left_dACC) <- c("PCP4", "POU3F1", "SULF2", "GABRQ")

Br8325_right_dACC_nucleus_area <- Br8325_right_dACC$Nucleus.Area..µm..
Br8325_right_dACC <- Br8325_right_dACC[,cols]
colnames(Br8325_right_dACC) <- c("PCP4", "POU3F1", "SULF2", "GABRQ")

# count cells that have zero for all counts
sum(rowSums(Br6432_dACC)==0)/dim(Br6432_dACC)[1]
sum(rowSums(Br6432_dlPFC)==0)/dim(Br6432_dlPFC)[1]
sum(rowSums(Br8325_dlPFC)==0)/dim(Br8325_dlPFC)[1]
sum(rowSums(Br8325_left_dACC)==0)/dim(Br8325_left_dACC)[1]
sum(rowSums(Br8325_right_dACC)==0)/dim(Br8325_right_dACC)[1]

# remove cells that have zero for all counts - they definitely are not VENs/we have no info
idx <- rowSums(Br6432_dACC)!=0
Br6432_dACC <- Br6432_dACC[idx,]
Br6432_dACC_nucleus_area <- Br6432_dACC_nucleus_area[idx]

idx <- rowSums(Br6432_dlPFC)!=0
Br6432_dlPFC <- Br6432_dlPFC[idx,]
Br6432_dlPFC_nucleus_area <- Br6432_dlPFC_nucleus_area[idx]

idx <- rowSums(Br8325_dlPFC)!=0
Br8325_dlPFC <- Br8325_dlPFC[idx,]
Br8325_dlPFC_nucleus_area <- Br8325_dlPFC_nucleus_area[idx]

idx <- rowSums(Br8325_left_dACC)!=0
Br8325_left_dACC <- Br8325_left_dACC[idx,]
Br8325_left_dACC_nucleus_area <- Br8325_left_dACC_nucleus_area[idx]

idx <- rowSums(Br8325_right_dACC)!=0
Br8325_right_dACC <- Br8325_right_dACC[idx,]
Br8325_right_dACC_nucleus_area <- Br8325_right_dACC_nucleus_area[idx]

# boxplots of raw counts distributions
datasets <- list(
    Br6432_dACC = Br6432_dACC,
    Br6432_dlPFC = Br6432_dlPFC,
    Br8325_dlPFC = Br8325_dlPFC,
    Br8325_left_dACC = Br8325_left_dACC,
    Br8325_right_dACC = Br8325_right_dACC
)

long_df <- bind_rows(
    lapply(names(datasets), function(name) {
        datasets[[name]] %>%
            mutate(dataset = name) %>%
            pivot_longer(
                cols = -dataset,
                names_to = "gene",
                values_to = "expression"
            )
    })
)

long_df$expression <- long_df$expression +1

ggplot(long_df, aes(x = gene, y = expression, color=dataset)) +
    geom_boxplot(position = position_dodge(0.8)) +
    theme_minimal() +
    labs(x = "Gene", y = "Expression (log10(x+1) scale)", color = "Dataset") +
    scale_fill_brewer(palette = "Set2") +
    scale_y_continuous(trans='log10')

long_df$expression <- long_df$expression -1

p1 <- ggplot(long_df[which(long_df$gene=="GABRQ"),], aes(x = gene, y = expression, color = dataset)) +
    geom_boxplot(position = position_dodge(0.8)) +
    theme_minimal() +
    labs(x = "Gene", y = "Expression", color = "Dataset") +
    scale_fill_brewer(palette = "Set2")
p2 <- ggplot(long_df[which(long_df$gene=="PCP4"),], aes(x = gene, y = expression, color = dataset)) +
    geom_boxplot(position = position_dodge(0.8)) +
    theme_minimal() +
    labs(x = "Gene", y = "Expression", color = "Dataset") +
    scale_fill_brewer(palette = "Set2")
p3 <- ggplot(long_df[which(long_df$gene=="POU3F1"),], aes(x = gene, y = expression, color = dataset)) +
    geom_boxplot(position = position_dodge(0.8)) +
    theme_minimal() +
    labs(x = "Gene", y = "Expression", color = "Dataset") +
    scale_fill_brewer(palette = "Set2")
p4 <- ggplot(long_df[which(long_df$gene=="SULF2"),], aes(x = gene, y = expression, color = dataset)) +
    geom_boxplot(position = position_dodge(0.8)) +
    theme_minimal() +
    labs(x = "Gene", y = "Expression", color = "Dataset") +
    scale_fill_brewer(palette = "Set2")

wrap_plots(p1,p2,p3,p4, nrow=2, guides="collect")

# boxplots of nucleus area
df <- data.frame(
    Area = c(Br6432_dACC_nucleus_area,
             Br6432_dlPFC_nucleus_area,
             Br8325_dlPFC_nucleus_area,
             Br8325_left_dACC_nucleus_area,
             Br8325_right_dACC_nucleus_area),
    Dataset = c(rep("Br6432_dACC", length(Br6432_dACC_nucleus_area)),
                rep("Br6432_dlPFC", length(Br6432_dlPFC_nucleus_area)),
                rep("Br8325_dlPFC", length(Br8325_dlPFC_nucleus_area)),
                rep("Br8325_left_dACC", length(Br8325_left_dACC_nucleus_area)),
                rep("Br8325_right_dACC", length(Br8325_right_dACC_nucleus_area)))
)

ggplot(df, aes(x = Dataset, y = Area, color=Dataset)) +
    geom_boxplot(position = position_dodge(0.8)) +
    theme_minimal() +
    labs(x = "Dataset", y = "Nucleus Area", color = "Dataset") +
    scale_fill_brewer(palette = "Set2") +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())

# similar to block 9
# https://github.com/LylaAtta123/normalization-analyses/blob/main/R/cosmx.ipynb
# for nucleus area normalization, scale by current nucleus area / median nucleus area
scale_Br6432_dACC <- Br6432_dACC_nucleus_area/median(Br6432_dACC_nucleus_area)
scale_Br6432_dlPFC <- Br6432_dlPFC_nucleus_area/median(Br6432_dlPFC_nucleus_area)
scale_Br8325_dlPFC <- Br8325_dlPFC_nucleus_area/median(Br8325_dlPFC_nucleus_area)
scale_Br8325_left_dACC <- Br8325_left_dACC_nucleus_area/median(Br8325_left_dACC_nucleus_area)
scale_Br8325_right_dACC <- Br8325_right_dACC_nucleus_area/median(Br8325_right_dACC_nucleus_area)

Br6432_dACC <- Br6432_dACC/scale_Br6432_dACC
Br6432_dlPFC <- Br6432_dlPFC/scale_Br6432_dlPFC
Br8325_dlPFC <- Br8325_dlPFC/scale_Br8325_dlPFC
Br8325_left_dACC <- Br8325_left_dACC/scale_Br8325_left_dACC
Br8325_right_dACC <- Br8325_right_dACC/scale_Br8325_right_dACC

# boxplots
long_df <- bind_rows(
    lapply(names(datasets), function(name) {
        datasets[[name]] %>%
            mutate(dataset = name) %>%
            pivot_longer(
                cols = -dataset,
                names_to = "gene",
                values_to = "expression"
            )
    })
)

long_df$expression <- long_df$expression +1

ggplot(long_df, aes(x = gene, y = expression, color=dataset)) +
    geom_boxplot(position = position_dodge(0.8)) +
    theme_minimal() +
    labs(x = "Gene", y = "Expression (log10(x+1) scale)", color = "Dataset") +
    scale_fill_brewer(palette = "Set2") +
    scale_y_continuous(trans='log10')

long_df$expression <- long_df$expression -1

p1 <- ggplot(long_df[which(long_df$gene=="GABRQ"),], aes(x = gene, y = expression, color = dataset)) +
    geom_boxplot(position = position_dodge(0.8)) +
    theme_minimal() +
    labs(x = "Gene", y = "Expression", color = "Dataset") +
    scale_fill_brewer(palette = "Set2")
p2 <- ggplot(long_df[which(long_df$gene=="PCP4"),], aes(x = gene, y = expression, color = dataset)) +
    geom_boxplot(position = position_dodge(0.8)) +
    theme_minimal() +
    labs(x = "Gene", y = "Expression", color = "Dataset") +
    scale_fill_brewer(palette = "Set2")
p3 <- ggplot(long_df[which(long_df$gene=="POU3F1"),], aes(x = gene, y = expression, color = dataset)) +
    geom_boxplot(position = position_dodge(0.8)) +
    theme_minimal() +
    labs(x = "Gene", y = "Expression", color = "Dataset") +
    scale_fill_brewer(palette = "Set2")
p4 <- ggplot(long_df[which(long_df$gene=="SULF2"),], aes(x = gene, y = expression, color = dataset)) +
    geom_boxplot(position = position_dodge(0.8)) +
    theme_minimal() +
    labs(x = "Gene", y = "Expression", color = "Dataset") +
    scale_fill_brewer(palette = "Set2")

wrap_plots(p1,p2,p3,p4, nrow=2, guides="collect")

# follow nucleus-area normalization with quantile normalization
# for each gene, quantiles by each sample

# create a df with genes expression of 1 gene for all 5 samples
genes <- colnames(Br6432_dACC)

# boxplots of raw counts distributions
datasets <- list(
    Br6432_dACC = Br6432_dACC,
    Br6432_dlPFC = Br6432_dlPFC,
    Br8325_dlPFC = Br8325_dlPFC,
    Br8325_left_dACC = Br8325_left_dACC,
    Br8325_right_dACC = Br8325_right_dACC
)

get_padded_matrix <- function(gene) {
    gene_list <- lapply(datasets, function(df) df[[gene]])
    max_len <- max(sapply(gene_list, length))
    gene_list_padded <- lapply(gene_list, function(x) {
        length(x) <- max_len  # pads with NA
        x
    })
    mat <- do.call(cbind, gene_list_padded)
    colnames(mat) <- names(datasets)
    return(mat)
}

gene_matrices <- setNames(lapply(genes, get_padded_matrix), genes)

# normalize
PCP4_norm <- normalize.quantiles(gene_matrices$PCP4)
POU3F1_norm <- normalize.quantiles(gene_matrices$POU3F1)
SULF2_norm <- normalize.quantiles(gene_matrices$SULF2)
GABRQ_norm <- normalize.quantiles(gene_matrices$GABRQ)

# put back into dfs
Br6432_dACC_norm <- data.frame(PCP4 = PCP4_norm[,1],
                               POU3F1 = POU3F1_norm[,1],
                               SULF2 = SULF2_norm[,1],
                               GABRQ = GABRQ_norm[,1])
Br6432_dACC_norm <- na.omit(Br6432_dACC_norm)

Br6432_dlPFC_norm <- data.frame(PCP4 = PCP4_norm[,2],
                               POU3F1 = POU3F1_norm[,2],
                               SULF2 = SULF2_norm[,2],
                               GABRQ = GABRQ_norm[,2])
Br6432_dlPFC_norm <- na.omit(Br6432_dlPFC_norm)

Br8325_dlPFC_norm <- data.frame(PCP4 = PCP4_norm[,3],
                               POU3F1 = POU3F1_norm[,3],
                               SULF2 = SULF2_norm[,3],
                               GABRQ = GABRQ_norm[,3])
Br8325_dlPFC_norm <- na.omit(Br8325_dlPFC_norm)

Br8325_left_dACC_norm <- data.frame(PCP4 = PCP4_norm[,4],
                                POU3F1 = POU3F1_norm[,4],
                                SULF2 = SULF2_norm[,4],
                                GABRQ = GABRQ_norm[,4])
Br8325_left_dACC_norm <- na.omit(Br8325_left_dACC_norm)

Br8325_right_dACC_norm <- data.frame(PCP4 = PCP4_norm[,5],
                                    POU3F1 = POU3F1_norm[,5],
                                    SULF2 = SULF2_norm[,5],
                                    GABRQ = GABRQ_norm[,5])
Br8325_right_dACC_norm <- na.omit(Br8325_right_dACC_norm)

# vis boxplots of expression after quantile normalization
datasets <- list(
    Br6432_dACC = Br6432_dACC_norm,
    Br6432_dlPFC = Br6432_dlPFC_norm,
    Br8325_dlPFC = Br8325_dlPFC_norm,
    Br8325_left_dACC = Br8325_left_dACC_norm,
    Br8325_right_dACC = Br8325_right_dACC_norm
)

long_df <- bind_rows(
    lapply(names(datasets), function(name) {
        datasets[[name]] %>%
            mutate(dataset = name) %>%
            pivot_longer(
                cols = -dataset,
                names_to = "gene",
                values_to = "expression"
            )
    })
)

long_df$expression <- long_df$expression +1

ggplot(long_df, aes(x = gene, y = expression, color=dataset)) +
    geom_boxplot(position = position_dodge(0.8)) +
    theme_minimal() +
    labs(x = "Gene", y = "Expression (log10(x+1) scale)", color = "Dataset") +
    scale_fill_brewer(palette = "Set2") +
    scale_y_continuous(trans='log10')

long_df$expression <- long_df$expression +1

p1 <- ggplot(long_df[which(long_df$gene=="GABRQ"),], aes(x = gene, y = expression, color = dataset)) +
    geom_boxplot(position = position_dodge(0.8)) +
    theme_minimal() +
    labs(x = "Gene", y = "Expression", color = "Dataset") +
    scale_fill_brewer(palette = "Set2") +
    scale_y_continuous(trans="log2")
p2 <- ggplot(long_df[which(long_df$gene=="PCP4"),], aes(x = gene, y = expression, color = dataset)) +
    geom_boxplot(position = position_dodge(0.8)) +
    theme_minimal() +
    labs(x = "Gene", y = "Expression", color = "Dataset") +
    scale_fill_brewer(palette = "Set2")+
    scale_y_continuous(trans="log2")
p3 <- ggplot(long_df[which(long_df$gene=="POU3F1"),], aes(x = gene, y = expression, color = dataset)) +
    geom_boxplot(position = position_dodge(0.8)) +
    theme_minimal() +
    labs(x = "Gene", y = "Expression", color = "Dataset") +
    scale_fill_brewer(palette = "Set2")+
    scale_y_continuous(trans="log2")
p4 <- ggplot(long_df[which(long_df$gene=="SULF2"),], aes(x = gene, y = expression, color = dataset)) +
    geom_boxplot(position = position_dodge(0.8)) +
    theme_minimal() +
    labs(x = "Gene", y = "Expression", color = "Dataset") +
    scale_fill_brewer(palette = "Set2")+
    scale_y_continuous(trans="log2")

wrap_plots(p1,p2,p3,p4, nrow=2, guides="collect")

save(Br6432_dACC_norm, Br6432_dlPFC_norm,
     Br8325_dlPFC_norm, Br8325_left_dACC_norm, Br8325_right_dACC_norm,
     file="normalized_dfs.Rdata")
