setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')
library(biomaRt)
library(zellkonverter)
library(SingleCellExperiment)
library(RcppML)
library(here)
library(tidyr)
library(forcats)
library(ggplot2)
library(dplyr)

##load the data
mch<-readH5AD(file=here::here('processed-data','snRNA-seq',
                              '06_NMF','rs2_mch_matrix.h5ad'))

##get gene names
mart <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
symb <- getBM(attributes = c("ensembl_gene_id","mgi_symbol"),
              filters = "ensembl_gene_id", values = rownames(mch),
              mart = mart)
symbs <- symb$mgi_symbol[match(rownames(mch), symb$ensembl_gene_id, nomatch = NA)]
sum(is.na(symbs))
#[1] 524

rowData(mch)$gene_name<-symbs
rowData(mch)$start<-NULL
rowData(mch)$end<-NULL

##drop genes with no names (need these for ortholog matching)
mch<-mch[!is.na(rowData(mch)$gene_name),]
rownames(mch)<-rowData(mch)$gene_name

# Translate from one mchcies to the other using the orthology
orthology<-read.csv(file=here::here('raw-data','Retro-seq',
                                    'human_mouse_orthologs.csv'))
names <- orthology[orthology$Column3 %in% rownames(mch),]

names <- names[match(rownames(mch), names$Column3),]

setdiff(names$Column3, rownames(mch))

sum(is.na(names$Column1))
# [1] 12406

rownames(mch) <- names$Column1
dim(mch)
# [1] 31684  3265

# remove rownames that are NA
mch <- mch[!is.na(rownames(mch)),]
dim(mch)
# [1] 19278  3265

# there are some repeated rownames in mch
duplicated_rows <- duplicated(rownames(mch))
mch <- mch[!duplicated_rows, ]
dim(mch)
# [1] 17473  3265

##load nmf patterns
x <- readRDS(file = here("processed-data", "snRNA-seq", "06_NMF", "nmf_results.RDS"))
loadings<-x@w

no_expr <- which(rowSums(loadings) == 0)
loadings <- loadings[-no_expr, ]
dim(loadings)
# 18565    75

set.seed(101)
i<-intersect(rownames(mch),rownames(loadings))
length(i)
# [1] 13343

loadings<-loadings[rownames(loadings) %in% i,]
mch<-mch[rownames(mch) %in% i,]

mch<-mch[,!is.na(mch$Subclass)]

## projection
loadings<-loadings[match(rownames(mch),rownames(loadings)),]
proj<-project(loadings,assay(mch,'X'),L1=0)
proj<-t(proj)

proj<-apply(proj,2,function(x){x/sum(x)})

colData(mch)<-cbind(colData(mch),proj)

saveRDS(mch, file = here("processed-data", "snRNA-seq", "06_NMF", "mch_projection.RDS"))

# visualize results

create_custom_dot_plot <- function(data, category_col, features_cols,
                                   plot_title, x_axis_title, y_axis_title,
                                   legend_size_title, legend_color_title) {

    # Ensuring that features_cols is a character vector
    features_cols <- as.character(features_cols)

    # Melting the data into long format
    long_data <- data %>%
        dplyr::select(!!sym(category_col), all_of(features_cols)) %>%
        pivot_longer(cols = -!!sym(category_col), names_to = "Feature", values_to = "Value")

    # Calculating sum and proportion of non-zero cells
    stats <- long_data %>%
        group_by(!!sym(category_col), Feature) %>%
        summarize(
            Sum = sum(Value),
            NonZeroProportion = sum(Value != 0) / n()  # Explicit proportion calculation
        ) %>%
        ungroup()

    # Adjusting Feature as a factor with the specified order and reversing the category_col order
    stats$Feature <- factor(stats$Feature, levels = features_cols)
    stats[[category_col]] <- fct_relevel(stats[[category_col]], rev)
    print(stats$NonZeroProportion)

    # Creating the plot
    ggplot(stats, aes(x = Feature, y = !!sym(category_col), size = NonZeroProportion, color = Sum)) +
        geom_point() +
        scale_size_continuous(range = c(0,10)) + # Set minimum size to zero
        scale_color_gradient(low = "white", high = "black") + # Greyscale color scale
        theme_minimal() +
        labs(
            title = plot_title,
            x = x_axis_title,
            y = y_axis_title,
            size = legend_size_title,
            color = legend_color_title
        ) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Adjust x-axis text angle for readability
}

mch <- readRDS(file = here("processed-data", "snRNA-seq", "06_NMF", "mch_projection.RDS"))
data<-as.data.frame(colData(mch))
# keep cols starting with "nmf"
nmf_cols <- grep("^nmf", colnames(data), value = TRUE)
data <- data[, c("Subclass", "Target",nmf_cols)]

# subset to only some nmf columns in the data df
# keep subclass and target columns
# we also want to rename the patterns to be more informative, such as "Astro-NMF32"
# Oligo: 26, 23, 27, 13, 43, 40, 36, 28, 9, 33, 39
# L5_6_NP: 46
# L6_b: 35
# L5_ET: 61
# L6_CT: 15
# L6_IT_Car3: 68
# L2_3_IT: 3, 11
# L5_IT: 38
# L6_IT: 32
# Pvalb: 10, 63
# SST: 52, 56
# SST Chodl: 51
# LAMP5: 37, 60
# Sncg: 55, 58
# Vip: 44, 47
# Endo: 75, 49
# Astro: 14, 21, 53, 65
# OPC: 17, 24
# VLMC: 59
# microPVM: 19, 54, 57
# misc: 64

# first subset the data to only the columns we want
nmf_cols <- c("nmf26", "nmf23", "nmf27", "nmf13", "nmf43", "nmf40", "nmf36", "nmf28", "nmf9", "nmf33", "nmf39",
              "nmf46", "nmf35", "nmf61", "nmf15", "nmf68", "nmf3", "nmf11", "nmf38", "nmf32", "nmf10", "nmf63",
              "nmf52", "nmf56", "nmf51", "nmf37", "nmf60", "nmf55", "nmf58", "nmf44", "nmf47", "nmf75", "nmf49",
              "nmf14", "nmf21", "nmf53", "nmf65", "nmf17", "nmf24", "nmf59", "nmf19", "nmf54", "nmf57",
              "nmf64")

data <- data[, c("Subclass", "Target", nmf_cols)]

# rename the columns
colnames(data) <- c("Subclass", "Target",
                    "Oligo-NMF26", "Oligo-NMF23", "Oligo-NMF27", "Oligo-NMF13", "Oligo-NMF43", "Oligo-NMF40", "Oligo-NMF36", "Oligo-NMF28", "Oligo-NMF9", "Oligo-NMF33", "Oligo-NMF39",
                    "L5_6_NP-NMF46", "L6_b-NMF35", "L5_ET-NMF61", "L6_CT-NMF15", "L6_IT_Car3-NMF68", "L2_3_IT-NMF3", "L2_3_IT-NMF11", "L5_IT-NMF38", "L6_IT-NMF32", "Pvalb-NMF10", "Pvalb-NMF63",
                    "SST-NMF52", "SST-NMF56", "SST Chodl-NMF51", "LAMP5-NMF37", "LAMP5-NMF60", "Sncg-NMF55", "Sncg-NMF58", "Vip-NMF44", "Vip-NMF47", "Endo-NMF75", "Endo-NMF49",
                    "Astro-NMF14", "Astro-NMF21", "Astro-NMF53", "Astro-NMF65", "OPC-NMF17", "OPC-NMF24", "VLMC-NMF59", "microPVM-NMF19", "microPVM-NMF54", "microPVM-NMF57",
                    "misc-NMF64")

# remove columns that are not excitatory
data[,c("Oligo-NMF26", "Oligo-NMF23", "Oligo-NMF27", "Oligo-NMF13", "Oligo-NMF43", "Oligo-NMF40", "Oligo-NMF36", "Oligo-NMF28", "Oligo-NMF9", "Oligo-NMF33", "Oligo-NMF39",
        "Pvalb-NMF10", "Pvalb-NMF63",
        "SST-NMF52", "SST-NMF56", "SST Chodl-NMF51", "LAMP5-NMF37", "LAMP5-NMF60", "Sncg-NMF55", "Sncg-NMF58", "Vip-NMF44", "Vip-NMF47", "Endo-NMF75", "Endo-NMF49",
        "Astro-NMF14", "Astro-NMF21", "Astro-NMF53", "Astro-NMF65", "OPC-NMF17", "OPC-NMF24", "VLMC-NMF59", "microPVM-NMF19", "microPVM-NMF54", "microPVM-NMF57",
        "misc-NMF64")] <- NULL

# remove columns that are not excitatory
data_subclass <- data[-which(data$Subclass %in% c("Lamp5 Gaba","Pvalb Gaba","Sst Gaba","Vip Gaba")),]

feat_cols <- colnames(data)[-c(1,2)]
feat_cols <- feat_cols[-c(7)]
feat_cols <- feat_cols[c(7,3,6,1,4,8,5,2)]

pdf(file=here::here('plots','snRNA-seq','06_NMF','mch_subclass_dotplot.pdf'),h=11,w=15)
create_custom_dot_plot(data_subclass, "Subclass", feat_cols, "", "NMF pattern",
                       "Allen subclass", "proportion nuclei\nwith nonzero\nweight",
                       "aggregate\nnuclei-level\nweights")+
    theme(axis.text=element_text(size=32,color='black'),text=element_text(size=32,color='black'))
dev.off()

pdf(file=here::here('plots','snRNA-seq','06_NMF','mch_target_dotplot.pdf'),h=11,w=15)
create_custom_dot_plot(data, "Target", feat_cols, "", "NMF pattern",
                       "Target", "proportion nuclei\nwith nonzero\nweight",
                       "aggregate\nnuclei-level\nweights")+
    theme(axis.text=element_text(size=32,color='black'),text=element_text(size=32,color='black'))
dev.off()


