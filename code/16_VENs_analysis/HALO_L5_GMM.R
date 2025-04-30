library(ggplot2)
library(tidyverse)
library(patchwork)
library(preprocessCore)
library(SingleCellExperiment)
library(factoextra)
library(cluster)
library(mclust)
library(scater)
library(ComplexHeatmap)
library(circlize)

setwd("../../../Downloads/Data")

df <- readxl::read_xlsx("Vens-Selected PCP4 POU3F1 SULF2 GABRQ selected.xlsx")

# create separate dfs

Br6432_dACC <- df[which(df$Algorithm.Name=="Vens dACC 6432.2"),]
Br6432_dlPFC <- df[which(df$Algorithm.Name=="Vens dlPFC 6432.2"),]
Br8325_dlPFC <- df[which(df$Algorithm.Name=="Vens dlPFC 8325.2"),]
Br8325_left_dACC <- df[which(df$Algorithm.Name=="Vens dACC left 8325.2"),]
Br8325_right_dACC <- df[which(df$Algorithm.Name=="Vens dACC right 8325.2"),]

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

# use GMM for each gene in each sample
set.seed(190)
for (gene in c("POU3F1","SULF2","GABRQ")) {
    print(gene)
    Br6432_dACC_mod <- Mclust(Br6432_dACC[,gene],3)
    Br6432_dACC[,paste0(gene,"_class")] <- Br6432_dACC_mod$classification
    print(summary(Br6432_dACC_mod, parameters=T))
    print(plot.Mclust(Br6432_dACC_mod, what="classification"), xlab=gene)

    Br6432_dlPFC_mod <- Mclust(Br6432_dlPFC[,gene],3)
    Br6432_dlPFC[,paste0(gene,"_class")] <- Br6432_dlPFC_mod$classification
    print(summary(Br6432_dlPFC_mod, parameters=T))
    print(plot.Mclust(Br6432_dlPFC_mod, what="classification"), xlab=gene)

    Br8325_dlPFC_mod <- Mclust(Br8325_dlPFC[,gene],3)
    Br8325_dlPFC[,paste0(gene,"_class")] <- Br8325_dlPFC_mod$classification
    print(summary(Br8325_dlPFC_mod, parameters=T))
    print(plot.Mclust(Br8325_dlPFC_mod, what="classification"), xlab=gene)

    Br8325_left_dACC_mod <- Mclust(Br8325_left_dACC[,gene],3)
    Br8325_left_dACC[,paste0(gene,"_class")] <- Br8325_left_dACC_mod$classification
    print(summary(Br8325_left_dACC_mod, parameters=T))
    print(plot.Mclust(Br8325_left_dACC_mod, what="classification"), xlab=gene)

    Br8325_right_dACC_mod <- Mclust(Br8325_right_dACC[,gene],)
    Br8325_right_dACC[,paste0(gene,"_class")] <- Br8325_right_dACC_mod$classification
    print(summary(Br8325_right_dACC_mod, parameters=T))
    print(plot.Mclust(Br8325_right_dACC_mod, what="classification"), xlab=gene)
}

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
            select(POU3F1_class, SULF2_class, GABRQ_class) %>%
            mutate(dataset = name)
    }),
    .id = "source"
) %>%
    pivot_longer(
        cols = c(POU3F1_class, SULF2_class, GABRQ_class),
        names_to = "gene",
        values_to = "class"
    )

prop_df <- long_df %>%
    group_by(dataset, gene, class) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(dataset, gene) %>%
    mutate(proportion = count / sum(count)) %>%
    ungroup()

prop_df$class <- as.factor(prop_df$class)

ggplot(prop_df[which(prop_df$class==2),], aes(x = gene, y = proportion, color=dataset)) +
    geom_jitter() +
    theme_bw() +
    labs(x = "Gene", y = "Proportion", title = "Classification Proportions Across Datasets and Genes") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

# plot distribution of expression for classification = 2

expr_list <- lapply(names(datasets), function(name) {
    dat <- datasets[[name]]

    # Create filtered data for each gene
    bind_rows(
        dat %>%
            filter(POU3F1_class == 3) %>%
            select(expression = POU3F1) %>%
            mutate(gene = "POU3F1", dataset = name),

        dat %>%
            filter(SULF2_class == 3) %>%
            select(expression = SULF2) %>%
            mutate(gene = "SULF2", dataset = name),

        dat %>%
            filter(GABRQ_class == 3) %>%
            select(expression = GABRQ) %>%
            mutate(gene = "GABRQ", dataset = name)
    )
})

expr_df <- bind_rows(expr_list)

ggplot(expr_df, aes(x = gene, y = expression, fill = dataset)) +
    geom_boxplot() +
    theme_bw() +
    labs(x = "Gene", y = "Expression", title = "Gene Expression (``Expressed`` Cells Only)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))


# compute coexpression summary
coexpression_summary <- data.frame(
    dataset = character(),
    POU3F1_SULF2 = integer(),
    POU3F1_GABRQ = integer(),
    GABRQ_SULF2 = integer(),
    all_3 = integer(),
    POU3F1_SULF2_prop = numeric(),
    POU3F1_GABRQ_prop = numeric(),
    GABRQ_SULF2_prop = numeric(),
    all_3_prop = numeric(),
    stringsAsFactors = FALSE
)
for (name in names(datasets)) {
    dat <- datasets[[name]]

    total_cells <- nrow(dat)

    # coexpression counts (both == 3)
    pou3f1_sulf2 <- sum(dat$POU3F1_class == 3 & dat$SULF2_class == 3, na.rm = TRUE)
    pou3f1_gabrq <- sum(dat$POU3F1_class == 3 & dat$GABRQ_class == 3, na.rm = TRUE)
    gabrq_sulf2 <- sum(dat$GABRQ_class == 3 & dat$SULF2_class == 3, na.rm = TRUE)
    all_three <- sum(dat$POU3F1_class == 3 & dat$SULF2_class == 3 & dat$GABRQ_class == 3, na.rm = TRUE)

    coexpression_summary <- rbind(coexpression_summary, data.frame(
        dataset = name,
        POU3F1_SULF2 = pou3f1_sulf2,
        POU3F1_GABRQ = pou3f1_gabrq,
        GABRQ_SULF2 = gabrq_sulf2,
        all_3 = all_three,
        POU3F1_SULF2_prop = round(pou3f1_sulf2 / total_cells, 5),
        POU3F1_GABRQ_prop = round(pou3f1_gabrq / total_cells,5),
        GABRQ_SULF2_prop = round(gabrq_sulf2 / total_cells,5),
        all_3_prop = round(all_three / total_cells,5)
    ))
}

print(coexpression_summary)

coexpression_summary_long <- pivot_longer(coexpression_summary,
                                          cols=POU3F1_SULF2_prop:all_3_prop)

ggplot(coexpression_summary_long, aes(x = name, y = value, color=dataset)) +
    geom_jitter() +
    theme_bw() +
    labs(x = "Group", y = "Proportion Coexpression", title = "Coexpression Proportions Across Datasets and Genes") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

# visualize coexpression as heatmap separated for dACC and dlPFC

SULF2_POU3F1 <- mean(c(dim(Br6432_dACC %>% filter(SULF2_class == 3) %>% filter(POU3F1_class == 3))[1] / dim(Br6432_dACC %>% filter(SULF2_class == 3))[1],
                       dim(Br8325_left_dACC %>% filter(SULF2_class == 3) %>% filter(POU3F1_class == 3))[1] / dim(Br8325_left_dACC %>% filter(SULF2_class == 3))[1],
                       dim(Br8325_right_dACC %>% filter(SULF2_class == 3) %>% filter(POU3F1_class == 3))[1] / dim(Br8325_right_dACC %>% filter(SULF2_class == 3))[1]))
GABRQ_POU3F1 <- mean(c(dim(Br6432_dACC %>% filter(GABRQ_class == 3) %>% filter(POU3F1_class == 3))[1] / dim(Br6432_dACC %>% filter(GABRQ_class == 3))[1],
                       dim(Br8325_left_dACC %>% filter(GABRQ_class == 3) %>% filter(POU3F1_class == 3))[1] / dim(Br8325_left_dACC %>% filter(GABRQ_class == 3))[1],
                       dim(Br8325_right_dACC %>% filter(GABRQ_class == 3) %>% filter(POU3F1_class == 3))[1] / dim(Br8325_right_dACC %>% filter(GABRQ_class == 3))[1]))
GABRQ_SULF2 <- mean(c(dim(Br6432_dACC %>% filter(GABRQ_class == 3) %>% filter(SULF2_class == 3))[1] / dim(Br6432_dACC %>% filter(GABRQ_class == 3))[1],
                      dim(Br8325_left_dACC %>% filter(GABRQ_class == 3) %>% filter(SULF2_class == 3))[1] / dim(Br8325_left_dACC %>% filter(GABRQ_class == 3))[1],
                      dim(Br8325_right_dACC %>% filter(GABRQ_class == 3) %>% filter(SULF2_class == 3))[1] / dim(Br8325_right_dACC %>% filter(GABRQ_class == 3))[1]))

POU3F1_SULF2 <- mean(c(dim(Br6432_dACC %>% filter(POU3F1_class == 3) %>% filter(SULF2_class == 3))[1] / dim(Br6432_dACC %>% filter(POU3F1_class == 3))[1],
                       dim(Br8325_left_dACC %>% filter(POU3F1_class == 3) %>% filter(SULF2_class == 3))[1] / dim(Br8325_left_dACC %>% filter(POU3F1_class == 3))[1],
                       dim(Br8325_right_dACC %>% filter(POU3F1_class == 3) %>% filter(SULF2_class == 3))[1] / dim(Br8325_right_dACC %>% filter(POU3F1_class == 3))[1]))
POU3F1_GABRQ <- mean(c(dim(Br6432_dACC %>% filter(POU3F1_class == 3) %>% filter(GABRQ_class == 3))[1] / dim(Br6432_dACC %>% filter(POU3F1_class == 3))[1],
                       dim(Br8325_left_dACC %>% filter(POU3F1_class == 3) %>% filter(GABRQ_class == 3))[1] / dim(Br8325_left_dACC %>% filter(POU3F1_class == 3))[1],
                       dim(Br8325_right_dACC %>% filter(POU3F1_class == 3) %>% filter(GABRQ_class == 3))[1] / dim(Br8325_right_dACC %>% filter(POU3F1_class == 3))[1]))
SULF2_GABRQ <- mean(c(dim(Br6432_dACC %>% filter(SULF2_class == 3) %>% filter(GABRQ_class == 3))[1] / dim(Br6432_dACC %>% filter(SULF2_class == 3))[1],
                      dim(Br8325_left_dACC %>% filter(SULF2_class == 3) %>% filter(GABRQ_class == 3))[1] / dim(Br8325_left_dACC %>% filter(SULF2_class == 3))[1],
                      dim(Br8325_right_dACC %>% filter(SULF2_class == 3) %>% filter(GABRQ_class == 3))[1] / dim(Br8325_right_dACC %>% filter(SULF2_class == 3))[1]))

dACC_mat <- matrix(c(NA, SULF2_POU3F1, GABRQ_POU3F1,
                     POU3F1_SULF2, NA, GABRQ_SULF2,
                     POU3F1_GABRQ, SULF2_GABRQ, NA), ncol=3)
rownames(dACC_mat) <- c("POU3F1","SULF2","GABRQ")
colnames(dACC_mat) <- c("POU3F1","SULF2","GABRQ")

col_fun <- colorRamp2(c(0, 0.5, 1), c("red", "green", "purple"))

Heatmap(dACC_mat, na_col = "black",
        cluster_rows = FALSE, cluster_columns = FALSE,
        col = col_fun, row_names_side = "left", column_names_side = "top",
        column_names_rot = 0, column_names_centered = T,
        heatmap_legend_param = list(title="dACC\nproportion\nof overlaps"),
        cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(sprintf("%.2f", dACC_mat[i, j]), x, y, gp = gpar(fontsize = 10))
        })

SULF2_POU3F1 <- mean(c(dim(Br6432_dlPFC %>% filter(SULF2_class == 3) %>% filter(POU3F1_class == 3))[1] / dim(Br6432_dlPFC %>% filter(SULF2_class == 3))[1],
                       dim(Br8325_dlPFC %>% filter(SULF2_class == 3) %>% filter(POU3F1_class == 3))[1] / dim(Br8325_dlPFC %>% filter(SULF2_class == 3))[1]))
GABRQ_POU3F1 <- mean(c(dim(Br6432_dlPFC %>% filter(GABRQ_class == 3) %>% filter(POU3F1_class == 3))[1] / dim(Br6432_dlPFC %>% filter(GABRQ_class == 3))[1],
                       dim(Br8325_dlPFC %>% filter(GABRQ_class == 3) %>% filter(POU3F1_class == 3))[1] / dim(Br8325_dlPFC %>% filter(GABRQ_class == 3))[1]))
GABRQ_SULF2 <- mean(c(dim(Br6432_dlPFC %>% filter(GABRQ_class == 3) %>% filter(SULF2_class == 3))[1] / dim(Br6432_dlPFC %>% filter(GABRQ_class == 3))[1],
                      dim(Br8325_dlPFC %>% filter(GABRQ_class == 3) %>% filter(SULF2_class == 3))[1] / dim(Br8325_dlPFC %>% filter(GABRQ_class == 3))[1]))
POU3F1_SULF2 <- mean(c(dim(Br6432_dlPFC %>% filter(POU3F1_class == 3) %>% filter(SULF2_class == 3))[1] / dim(Br6432_dlPFC %>% filter(POU3F1_class == 3))[1],
                       dim(Br8325_dlPFC %>% filter(POU3F1_class == 3) %>% filter(SULF2_class == 3))[1] / dim(Br8325_dlPFC %>% filter(POU3F1_class == 3))[1]))
POU3F1_GABRQ <- mean(c(dim(Br6432_dlPFC %>% filter(POU3F1_class == 3) %>% filter(GABRQ_class == 3))[1] / dim(Br6432_dlPFC %>% filter(POU3F1_class == 3))[1],
                       dim(Br8325_dlPFC %>% filter(POU3F1_class == 3) %>% filter(GABRQ_class == 3))[1] / dim(Br8325_dlPFC %>% filter(POU3F1_class == 3))[1]))
SULF2_GABRQ <- mean(c(dim(Br6432_dlPFC %>% filter(SULF2_class == 3) %>% filter(GABRQ_class == 3))[1] / dim(Br6432_dlPFC %>% filter(SULF2_class == 3))[1],
                      dim(Br8325_dlPFC %>% filter(SULF2_class == 3) %>% filter(GABRQ_class == 3))[1] / dim(Br8325_dlPFC %>% filter(SULF2_class == 3))[1]))


dlPFC_mat <- matrix(c(NA, SULF2_POU3F1, GABRQ_POU3F1,
                      POU3F1_SULF2, NA, GABRQ_SULF2,
                      POU3F1_GABRQ, SULF2_GABRQ, NA), ncol=3)
rownames(dlPFC_mat) <- c("POU3F1","SULF2","GABRQ")
colnames(dlPFC_mat) <- c("POU3F1","SULF2","GABRQ")

col_fun <- colorRamp2(c(0, 0.5, 1), c("red", "green", "purple"))

Heatmap(dlPFC_mat, na_col = "black",
        cluster_rows = FALSE, cluster_columns = FALSE,
        col = col_fun, row_names_side = "left", column_names_side = "top",
        column_names_rot = 0, column_names_centered = T,
        heatmap_legend_param = list(title="dlPFC\nproportion\nof overlaps"),
        cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(sprintf("%.2f", dlPFC_mat[i, j]), x, y, gp = gpar(fontsize = 10))
        })



# compute correlations
# make df of correlations
corr_list <- c(cor(Br6432_dACC$POU3F1, Br6432_dACC$SULF2, method="spearman"),
               cor(Br6432_dACC$POU3F1, Br6432_dACC$GABRQ, method="spearman"),
               cor(Br6432_dACC$SULF2, Br6432_dACC$GABRQ, method="spearman"),
               cor(Br6432_dlPFC$POU3F1, Br6432_dlPFC$SULF2, method="spearman"),
               cor(Br6432_dlPFC$POU3F1, Br6432_dlPFC$GABRQ, method="spearman"),
               cor(Br6432_dlPFC$SULF2, Br6432_dlPFC$GABRQ, method="spearman"),
               cor(Br8325_dlPFC$POU3F1, Br8325_dlPFC$SULF2, method="spearman"),
               cor(Br8325_dlPFC$POU3F1, Br8325_dlPFC$GABRQ, method="spearman"),
               cor(Br8325_dlPFC$SULF2, Br8325_dlPFC$GABRQ, method="spearman"),
               cor(Br8325_left_dACC$POU3F1, Br8325_left_dACC$SULF2, method="spearman"),
               cor(Br8325_left_dACC$POU3F1, Br8325_left_dACC$GABRQ, method="spearman"),
               cor(Br8325_left_dACC$SULF2, Br8325_left_dACC$GABRQ, method="spearman"),
               cor(Br8325_right_dACC$POU3F1, Br8325_right_dACC$SULF2, method="spearman"),
               cor(Br8325_right_dACC$POU3F1, Br8325_right_dACC$GABRQ, method="spearman"),
               cor(Br8325_right_dACC$SULF2, Br8325_right_dACC$GABRQ, method="spearman"))

region_list <- c(rep("dACC",3), rep("dlPFC",3), rep("dlPFC",3),
                 rep("dACC",3), rep("dACC",3))

genes_list <- c(rep(c("POU3F1 & SULF2","POU3F1 & GABRQ","GABRQ & SULF2"),5))

df <- data.frame(
    Correlation = corr_list,
    Region = region_list,
    Genes = genes_list
)

ggplot(df, aes(x=Genes, y=Correlation, color=Region)) +
    geom_point() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust=1))

