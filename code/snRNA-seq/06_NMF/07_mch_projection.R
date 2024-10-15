setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')
library(biomaRt)
library(zellkonverter)
library(SingleCellExperiment)
library(RcppML)
library(here)
library(tidyr)
library(forcats)
library(ggplot2)

##load the data
mch<-readH5AD(file=here::here('processed-data','snRNA-seq',
                              '06_NMF','rs2_mch_matrix.h5ad'))

##get gene names
mart <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
symb <- getBM(attributes = c("ensembl_gene_id","mgi_symbol"),
              filters = "ensembl_gene_id", values = rownames(mch),
              mart = mart)
symbs <- symb$mgi_symbol[match(rownames(mch), symb$ensembl_gene_id, nomatch = NA)]
count(is.na(symbs))
#[1] 165

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

count(is.na(names$Column1))
# [1] 12677

rownames(mch) <- names$Column1
dim(mch)
# [1] 32043  3265

# remove rownames that are NA
mch <- mch[!is.na(rownames(mch)),]
dim(mch)
# [1] 19366  3265

# there are some repeated rownames in mch
duplicated_rows <- duplicated(rownames(mch))
mch <- mch[!duplicated_rows, ]
dim(mch)
# [1] 17490  3265

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
# [1] 13354

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

pdf(file=here::here('plots','snRNA-seq','06_NMF','mch_subclass_dotplot.pdf'),h=10,w=45)
create_custom_dot_plot(data, "Subclass", nmf_cols, "", "NMF pattern",
                       "Allen subclass", "proportion nuclei\nwith nonzero\nweight",
                       "aggregate\nnuclei-level\nweights")+
    theme(axis.text=element_text(size=32,color='black'),text=element_text(size=32,color='black'))
dev.off()

pdf(file=here::here('plots','snRNA-seq','06_NMF','mch_target_dotplot.pdf'),h=10,w=40)
create_custom_dot_plot(data, "Target", nmf_cols, "", "NMF pattern",
                       "Target", "proportion nuclei\nwith nonzero\nweight",
                       "aggregate\nnuclei-level\nweights")+
    theme(axis.text=element_text(size=32,color='black'),text=element_text(size=32,color='black'))
dev.off()


