setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

# Load libraries
library(SpatialExperiment)
library(RcppML)
library(here)
library(UpSetR)
library(ComplexHeatmap)

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
select.nmfs = c("nmf3","nmf61","nmf38","nmf46","nmf15","nmf32","nmf68","nmf35")
nmf.genes = c("ASIC2", "RBFOX1", "ZNF385D",
              "ROBO2", "DPP10", "KIAA1217",
              "KCNIP4", "ROBO2",
              "ROBO2", "DPP10", "LRP1B",
              "KCNIP4", "LSAMP", "CSMD1",
              "KCNIP4", "EPHA6", "LRRTM4",
              "ROBO2", "KCNIP4", "AC117453.1",
              "AC109466.1", "KCNIP4", "IL1RAPL2",
              "ROBO2", "CDH13", "LRP1B",
              "VAT1L")

nmf.genes <- unique(nmf.genes)

m1 = loads[nmf.genes,select.nmfs]

colnames(m1) <- c("L2_3_IT-NMF3", "L5_ET-NMF61", "L5_IT-NMF38", "L5_6_NP-NMF46", "L6_CT-NMF15",
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


# heatmap of top gene for each excitatory NMF pattern
loads<-x@w
no_expr <- which(rowSums(loads) == 0)
loads <- loads[-no_expr, ]

# create heatmap of selected genes
select.nmfs = c("nmf3","nmf61","nmf38","nmf46","nmf15","nmf32","nmf68","nmf35")

# i want to get the genes that are unique (not duplicated)
top <- top_genes(x$w, 100)
subset_genes <- top[,select.nmfs]
list_subset_genes <- as.vector(subset_genes)
nmf.genes <- names(table(list_subset_genes))[table(list_subset_genes) == 1]
nmf.genes <- append(c("VAT1L", "ADRA1A", "GABRQ", "POU3F1"), nmf.genes)

group_vec <- c()

for (gene in nmf.genes) {
    for (factor_var in select.nmfs) {
        print(factor_var)
        if(gene %in% subset_genes[,factor_var]){
            group_vec <- append(group_vec, factor_var)
        }
    }
}

group_vec <- append(rep("nmf61",4), group_vec)

m1 = loads[nmf.genes,select.nmfs]

colnames(m1) <- c("L2_3_IT-NMF3", "L5_ET-NMF61", "L5_IT-NMF38", "L5_6_NP-NMF46", "L6_CT-NMF15",
                  "L6_IT-NMF32", "L6_IT_Car3-NMF68","L6b-NMF35")

m1 <- t(scale(t(m1)))

dend = cluster_within_group(t(m1),group_vec)

pdf(here("plots", "snRNA-seq", "06_NMF", "NMF_top_genes_heatmap_100.pdf"), height = 3, width = 20)
ComplexHeatmap::Heatmap(t(m1), show_column_dend = FALSE, show_row_dend = FALSE,
                        cluster_columns = dend, cluster_rows = F,
                        column_names_rot = 90, row_names_side = "left",
                        top_annotation = HeatmapAnnotation(group=group_vec,
                                                           col = list(group=c("nmf15"="#FB6496",
                                                                      "nmf3"="#1F78C8",
                                                                      "nmf32"="#FDBF6F",
                                                                      "nmf35"="#33a02c",
                                                                      "nmf38"="#6A33C2",
                                                                      "nmf46"="#0000FF",
                                                                      "nmf61"="#C8308C",
                                                                      "nmf68"="#C814FA"))),
                        heatmap_legend_param = list(
                            title = "loadings\nscaled by\ngene", at = c(-4, 0, 4),
                            labels = c("-4", "0", "4")
                        ))
dev.off()

# manually subsetting

loads<-x@w
no_expr <- which(rowSums(loads) == 0)
loads <- loads[-no_expr, ]

# create heatmap of selected genes
select.nmfs = c("nmf3","nmf61","nmf38","nmf46","nmf15","nmf32","nmf68","nmf35")

nmf.genes <- c(
                "GNAL","ROBO1","SGCD","GRIA2","CACNB2",
                "VAT1L", "ADRA1A", "GABRQ", "POU3F1","COL5A2",
                "LINC02055","IL1RAPL2","RORB","CASC15","AC008415.1",
                "ITGA8","HTR2C","NPSR1-AS1","TSHZ2","DCC",
                "ADAMTSL1","AC011246.1","SEMA5A","TRPM3","TMEFF2",
                "AC007368.1","NCAM2","SGCZ","SLIT3","PTPRK",
                "SEMA6D","RGS12","CUX1","SYNPR","ZNF804B",
                "AL136456.1","PCSK5","FUT9","PCDH11X","ZFHX3")

group_vec <- c(
               rep("NMF3",5), rep("NMF61",5),
               rep("NMF38",5), rep("NMF46",5),
               rep("NMF15",5), rep("NMF32",5),
               rep("NMF68",5), rep("NMF35",5))

m1 <- loads[nmf.genes,select.nmfs]

colnames(m1) <- c("L2_3_IT-NMF3", "L5_ET-NMF61", "L5_IT-NMF38", "L5_6_NP-NMF46", "L6_CT-NMF15",
                  "L6_IT-NMF32", "L6_IT_Car3-NMF68","L6b-NMF35")

m1 <- t(scale(t(m1)))

dend <- cluster_within_group(t(m1),group_vec)

pdf(here("plots", "snRNA-seq", "06_NMF", "NMF_top_genes_heatmap_5.pdf"), height = 3, width = 10)
ComplexHeatmap::Heatmap(t(m1), show_column_dend = FALSE, show_row_dend = FALSE,
                        cluster_columns = F, cluster_rows = F,
                        column_names_rot = 90, row_names_side = "left",
                        top_annotation = HeatmapAnnotation(factor=group_vec,
                                                           col = list(factor=c("NMF15"="#FB6496",
                                                                              "NMF3"="#1F78C8",
                                                                              "NMF32"="#FDBF6F",
                                                                              "NMF35"="#33a02c",
                                                                              "NMF38"="#6A33C2",
                                                                              "NMF46"="#0000FF",
                                                                              "NMF61"="#C8308C",
                                                                              "NMF68"="#C814FA"))),
                        heatmap_legend_param = list(
                            title = "loadings\nscaled by\ngene", at = c(-4, 0, 4),
                            labels = c("-4", "0", "4")
                        ))
dev.off()
