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
select.nmfs = c("nmf46","nmf35","nmf61","nmf15","nmf68","nmf3","nmf11","nmf38","nmf32")
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


nmf.genes = c("ASIC2",
              "KCNIP4",
              "VAT1L",
              "ROBO2",
              "AC109466.1")

nmf.genes <- unique(nmf.genes)

m1 = loads[nmf.genes,select.nmfs]
pdf(here("plots", "snRNA-seq", "06_NMF", "NMF_top_genes_heatmap.pdf"))
pheatmap::pheatmap(m1, scale="row", cluster_rows = F, cluster_cols = F,
                   angle_col=0)
dev.off()



