setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

# Load libraries
library(SpatialExperiment)
library(RcppML)
library(here)
library(UpSetR)

load(file = here("processed-data", "snRNA-seq", "03_batch_correction", "sce_harmony.Rdata"))
x <- readRDS(file = here("processed-data", "snRNA-seq", "06_NMF", "nmf_results.RDS"))

# function for getting top n genes for each pattern
top_genes <- function(W, n=10){
    top_genes <- apply(W, 2, function(x) names(sort(x, decreasing=TRUE)[1:n]))
    return(top_genes)
}

# get top 10 genes
top10 <- top_genes(x$w, 10)
write.csv(top10, file = here("processed-data", "snRNA-seq", "06_NMF", "top10_genes.csv"))

#check if there are overlapping genes in columns of "nmf", 26, 23, 27, 13, 43, 40, 36, 28, 9, 33, 39
cols <- c("nmf26", "nmf23", "nmf27", "nmf13", "nmf43", "nmf40", "nmf36", "nmf28", "nmf9", "nmf33", "nmf39")

# create UpSet plot
subset_top10 <- top10[cols]
upset_data <- fromList(subset_top10)
pdf(file = here("plots", "snRNA-seq", "06_NMF", "upset_top10_genes.pdf"), width = 10, height = 10)
upset(upset_data, sets = names(subset_top10), order.by = "freq", main.bar.color = "blue", matrix.color = "darkred")
grid.text("UpSet of dACC single nucleus oligo patterns",x = 0.65, y=0.95, gp=gpar(fontsize=20))
dev.off()

# get the overlapping genes for each pair of patterns
overlapping_genes <- list()
for(i in 1:length(venn_cols)){
    for(j in 1:length(venn_cols)){
        if(i != j){
            overlapping_genes[[paste(venn_cols[i], venn_cols[j], sep="_")]] <- intersect(top10[[venn_cols[i]]], top10[[venn_cols[j]]])
        }
    }
}


# heatmap of top gene for each excitatory NMF pattern
loads<-x@w
no_expr <- which(rowSums(loads) == 0)
loads <- loads[-no_expr, ]

# create heatmap of selected genes
select.nmfs = c("nmf11","nmf61","nmf38","nmf46","nmf15","nmf32","nmf68","nmf35")
nmf.genes = c("ASIC2", "RBFOX1", "ZNF385D",
              "ROBO2", "DPP10", "KIAA1217",
              "KCNIP4", "ROBO2", "MALAT1",
              "ROBO2", "DPP10", "LRP1B",
              "KCNIP4", "LSAMP", "CSMD1",
              "KCNIP4", "EPHA6", "LRRTM4",
              "ROBO2", "KCNIP4", "AC117453.1",
              "AC109466.1", "KCNIP4", "IL1RAPL2",
              "ROBO2", "CDH13", "LRP1B",
              "VAT1L")

nmf.genes <- unique(nmf.genes)

m1 = loads[nmf.genes,select.nmfs]

colnames(m1) <- c("L2_3_IT-NMF11", "L5_ET-NMF61", "L5_IT-NMF38", "L5_6_NP-NMF46", "L6_CT-NMF15",
                  "L6_IT-NMF32", "L6_IT_Car3-NMF68","L6b-NMF35")

m1 <- t(scale(t(m1)))

pdf(here("plots", "snRNA-seq", "06_NMF", "NMF_top_genes_heatmap.pdf"), height = 4, width = 4)
ComplexHeatmap::Heatmap(m1, show_column_dend = FALSE, show_row_dend = FALSE,
                        column_names_rot = 45, row_names_side = "left",
                        heatmap_legend_param = list(
                            title = "loadings\nscaled by\ngene", at = c(-4, 0, 4),
                            labels = c("-4", "0", "4")
                        ))
dev.off()



