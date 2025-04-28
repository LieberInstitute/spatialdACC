library(ggplot2)
library(tidyverse)
library(patchwork)
library(preprocessCore)
library(SingleCellExperiment)
library(factoextra)
library(cluster)
library(mclust)
library(scater)

setwd("../../../Downloads/Data")

load(file="nucleus_normalized_dfs.Rdata")
set.seed(12)
# use GMM for each gene in each sample
for (gene in c("POU3F1","SULF2","GABRQ")) {
    print(gene)
    Br6432_dACC_mod <- Mclust(Br6432_dACC[,gene],2)
    Br6432_dACC[,paste0(gene,"_class")] <- Br6432_dACC_mod$classification
    print(summary(Br6432_dACC_mod, parameters=T))
    print(plot.Mclust(Br6432_dACC_mod, what="classification"), xlab=gene)

    Br6432_dlPFC_mod <- Mclust(Br6432_dlPFC[,gene],2)
    Br6432_dlPFC[,paste0(gene,"_class")] <- Br6432_dlPFC_mod$classification
    print(summary(Br6432_dlPFC_mod, parameters=T))
    print(plot.Mclust(Br6432_dlPFC_mod, what="classification"), xlab=gene)

    Br8325_dlPFC_mod <- Mclust(Br8325_dlPFC[,gene],2)
    Br8325_dlPFC[,paste0(gene,"_class")] <- Br8325_dlPFC_mod$classification
    print(summary(Br8325_dlPFC_mod, parameters=T))
    print(plot.Mclust(Br8325_dlPFC_mod, what="classification"), xlab=gene)

    Br8325_left_dACC_mod <- Mclust(Br8325_left_dACC[,gene],2)
    Br8325_left_dACC[,paste0(gene,"_class")] <- Br8325_left_dACC_mod$classification
    print(summary(Br8325_left_dACC_mod, parameters=T))
    print(plot.Mclust(Br8325_left_dACC_mod, what="classification"), xlab=gene)

    Br8325_right_dACC_mod <- Mclust(Br8325_right_dACC[,gene],2)
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
            filter(POU3F1_class == 2) %>%
            select(expression = POU3F1) %>%
            mutate(gene = "POU3F1", dataset = name),

        dat %>%
            filter(SULF2_class == 2) %>%
            select(expression = SULF2) %>%
            mutate(gene = "SULF2", dataset = name),

        dat %>%
            filter(GABRQ_class == 2) %>%
            select(expression = GABRQ) %>%
            mutate(gene = "GABRQ", dataset = name)
    )
})

# Combine into one long data frame
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

    # coexpression counts (both == 2)
    pou3f1_sulf2 <- sum(dat$POU3F1_class == 2 & dat$SULF2_class == 2, na.rm = TRUE)
    pou3f1_gabrq <- sum(dat$POU3F1_class == 2 & dat$GABRQ_class == 2, na.rm = TRUE)
    gabrq_sulf2 <- sum(dat$GABRQ_class == 2 & dat$SULF2_class == 2, na.rm = TRUE)
    all_three <- sum(dat$POU3F1_class == 2 & dat$SULF2_class == 2 & dat$GABRQ_class == 2, na.rm = TRUE)

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
    geom_point() +
    theme_bw() +
    labs(x = "Group", y = "Proportion Coexpression", title = "Coexpression Proportions Across Datasets and Genes") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))


