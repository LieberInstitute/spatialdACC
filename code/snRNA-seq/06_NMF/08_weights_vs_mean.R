setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

# Load libraries
library(SpatialExperiment)
library(RcppML)
library(here)
library(dplyr)
library(ggplot2)
library(scran)

load(file = here("processed-data", "snRNA-seq", "03_batch_correction", "sce_harmony.Rdata"))
x <- readRDS(file = here("processed-data", "snRNA-seq", "06_NMF", "nmf_results.RDS"))

sce <- logNormCounts(sce)

avg.expr = rowMeans(logcounts(sce))
loads<-x@w
no_expr <- which(rowSums(loads) == 0)
loads <- loads[-no_expr, ]
avg.expr = avg.expr[rownames(loads)]

df1 = cbind.data.frame(loads, avg.expr)

plist1 <- lapply(colnames(loads), function(x) {
    tmp = df1[,c(x,"avg.expr")]
    colnames(tmp) = c("nmf","avg.expr")
    tmp$rank = (1+nrow(loads))-rank(tmp$nmf)
    tmp$top186 = tmp$rank<=186
    ggplot(tmp, aes(x=avg.expr, y=nmf, color=top186))+
        geom_point(size=.1)+theme_bw()+
        geom_point(data=filter(tmp, top186==T), size=.3)+
        scale_color_manual(values=c("black","red"))+
        labs(title=x, y="gene weights", x="avg. expr")+
        theme(text=element_text(size=8), #axis.title.y=element_blank()),
              axis.text=element_blank(), legend.position="none")
})

ggsave(filename=here("plots", "snRNA-seq", "06_NMF","nmf-weight_gene-expr_scatter_top186-red-enlarge.png"),
       plot=do.call(gridExtra::grid.arrange, c(plist1, ncol=10)),
       width=14, height=14, units="in", bg="white")


plist1 <- lapply(colnames(loads), function(x) {
    tmp = df1[,c(x,"avg.expr")]
    colnames(tmp) = c("nmf","avg.expr")
    tmp$rank = (1+nrow(loads))-rank(tmp$nmf)
    tmp$top928 = tmp$rank<=928
    ggplot(tmp, aes(x=avg.expr, y=nmf, color=top928))+
        geom_point(size=.1)+theme_bw()+
        geom_point(data=filter(tmp, top928==T), size=.3)+
        scale_color_manual(values=c("black","red"))+
        labs(title=x, y="gene weights", x="avg. expr")+
        theme(text=element_text(size=8), #axis.title.y=element_blank()),
              axis.text=element_blank(), legend.position="none")
})

ggsave(filename=here("plots", "snRNA-seq", "06_NMF","nmf-weight_gene-expr_scatter_top928-red-enlarge.png"),
       plot=do.call(gridExtra::grid.arrange, c(plist1, ncol=10)),
       width=14, height=14, units="in", bg="white")

setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

# Load libraries
library(RcppML)
library(here)
library(tidyverse)

loads <- readRDS(file = here("processed-data", "snRNA-seq", "06_NMF", "nmf_patterns_subset.RDS"))
no_expr <- which(rowSums(loads) == 0)
loads <- loads[-no_expr, ]

ldsc.score <- as.matrix(read.csv(here::here("processed-data", "17_LDSC", "NMF_score.csv"),
                                 row.names = 1))

# make the last underscore - in colnames(ldsc.score) to .
colnames(loads) <- gsub("-", ".", colnames(loads))

colnames(ldsc.score)[25] <- "SST_Chodl.NMF51"

#### new nmf_score_top-rank.csv
identical(rownames(loads), rownames(ldsc.score))
identical(colnames(loads), colnames(ldsc.score))
topN.mat = matrix(0, nrow=nrow(ldsc.score), ncol=ncol(ldsc.score))


rownames(topN.mat) <- rownames(ldsc.score)
colnames(topN.mat) <- colnames(ldsc.score)

for(i in colnames(ldsc.score)) {
    tmp = loads[,i]
    rank1 = (1+length(tmp))-rank(tmp)
    tmp.bin = ifelse(rank1<=928, 1, 0)
    stopifnot(identical(names(tmp.bin), rownames(topN.mat)))
    topN.mat[,i] = tmp.bin
}

fivenum(rowSums(topN.mat))

pdf(here("plots", "snRNA-seq", "06_NMF","nmf_score_top928.pdf"))
hist(rowSums(topN.mat), main="Number of repeats for each gene", sub = "11,911 zeroes")
hist(rowSums(topN.mat)[which(rowSums(topN.mat)>0)], breaks= 10, main="Number of repeats for each gene - excluding 0",
     sub = "1,738 ones, 3,626 < 5 repeats")
dev.off()

# in addition to exploring whether or not there are genes that overlap across patterns
# could just be a yes/no all 1 category or multiple categories
# stacked barplots showing the percentages of y/n as the number of overlaps increases
cell_types <- sub("\\.NMF\\d+$", "", colnames(topN.mat))

result <- topN.mat %>%
    as.data.frame() %>%
    mutate(Gene = rownames(.)) %>%
    pivot_longer(-Gene, names_to = "Pattern", values_to = "Count") %>%
    mutate(CellType = sub("\\.NMF\\d+$", "", Pattern)) %>%
    group_by(Gene, CellType) %>%
    summarize(Sum = sum(Count), .groups = "drop") %>%
    pivot_wider(names_from = CellType, values_from = Sum, values_fill = 0)

result <- result %>%
    mutate(TotalOverlaps = rowSums(across(-Gene)),
           CategoryCount = rowSums(select(., -Gene) > 0))

result <- result %>%
    filter(TotalOverlaps > 1)

result <- result %>%
    mutate(Classification = ifelse(CategoryCount == 1, "Single Category", "Multiple Categories"))

plot_data <- result %>%
    group_by(TotalOverlaps, Classification) %>%
    summarize(NumGenes = n(), .groups = "drop")

p <- ggplot(plot_data, aes(x = factor(TotalOverlaps), y = NumGenes, fill = Classification)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = c("Single Category" = "skyblue", "Multiple Categories" = "orange")) +
    labs(
        x = "Number of Overlaps",
        y = "Number of Genes",
        fill = "Overlap Classification",
        title = "Distribution of Gene Overlaps by Category"
    ) +
    theme_minimal()

pdf(here("plots", "snRNA-seq", "06_NMF","nmf_categories_top928.pdf"), height=15,width=15)
print(p)
dev.off()

