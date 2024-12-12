setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

# Load libraries
library(here)
library(dplyr)
library(ggplot2)
library(fgsea)
library(reactome.db)
library(org.Hs.eg.db)

x <- readRDS(file = here("processed-data", "snRNA-seq", "06_NMF", "nmf_results.RDS"))

loads<-x@w
no_expr <- which(rowSums(loads) == 0)
loads <- loads[-no_expr, ]

# prep reactome
xx <- as.list(reactomePATHID2NAME)
reactome.h <- xx[grep("^Homo",xx)]
x <- as.list(reactomePATHID2EXTID)
reactome.h = reactome.h[intersect(names(reactome.h), names(x))]
x.h <- x[names(reactome.h)]
identical(names(x.h), names(reactome.h))
reactome.h = gsub("Homo sapiens: ","",reactome.h)
names(x.h) = reactome.h

# Oligo: 26, 23, 27, 13, 43, 40, 36, 28, 9, 33, 39

nmf_pattern <- "nmf26"

non0.nmf26 = rownames(loads)[loads[,"nmf26"]>0]

# can also get the top 928 genes to do GSEA analysis on
#tmp = loads[,nmf_pattern]
#rank1 = (1+length(tmp))-rank(tmp)
#non0.nmf26 <- names(rank1)[rank1 <= 928]

non0.26.id = mapIds(org.Hs.eg.db, keys=non0.nmf26, keytype="SYMBOL", column="ENTREZID", multiVals = "first")
names(non0.26.id) = non0.nmf26
non0.26.id = non0.26.id[!is.na(non0.26.id)]
pathways.26 <- reactomePathways(non0.26.id)
pathways.26 <- x.h[names(pathways.26)]

nmf26.stats = loads[names(non0.26.id),"nmf26"]
names(nmf26.stats) = non0.26.id
#in order to get reproducible results, must call set.seed every time before you run the function
### the function has an internal set seed that makes this necessary
set.seed(123)
olig.results.26 = fgseaMultilevel(pathways.26, stats=nmf26.stats, scoreType="pos", minSize=15, maxSize=500)
olig.results.26$leadingEdge2 = sapply(olig.results.26$leadingEdge, paste, collapse="/")
small_dat <- olig.results.26[,c(1,3)]
print(small_dat[order(small_dat$padj),][c(1:10),])
print(small_dat[order(small_dat$padj),][c(1:3),]$pathway)
write.csv(olig.results.26[,c(1:7,9)], here("processed-data", "snRNA-seq", "06_NMF", paste0(nmf_pattern,"_reactome_results_sig.csv")))


# plot the top 3 pathways for nmf26
pdf(here("plots", "snRNA-seq", "06_NMF", "GSEA", paste0(nmf_pattern,".pdf")))
nmf26.terms = small_dat[order(small_dat$padj),][c(1:3),]$pathway
plotGseaTable(x.h[nmf26.terms],
              nmf26.stats, olig.results.26,
              gseaParam = 0.5)
dev.off()

for (i in c(26, 23, 27, 13, 43, 40, 36, 28, 9, 33, 39)){
  nmf_pattern <- paste0("nmf",i)
  non0.nmf = rownames(loads)[loads[,nmf_pattern]>0]
  print(paste0(nmf_pattern, " - ", length(non0.nmf)))
  non0.id = mapIds(org.Hs.eg.db, keys=non0.nmf, keytype="SYMBOL", column="ENTREZID", multiVals = "first")
  names(non0.id) = non0.nmf
  non0.id = non0.id[!is.na(non0.id)]
  pathways <- reactomePathways(non0.id)
  pathways <- x.h[names(pathways)]
  nmf.stats = loads[names(non0.id),nmf_pattern]
  names(nmf.stats) = non0.id
  set.seed(123)
  olig.results = fgseaMultilevel(pathways, stats=nmf.stats, scoreType="pos", minSize=15, maxSize=500)
  olig.results$leadingEdge2 = sapply(olig.results$leadingEdge, paste, collapse="/")
  small_dat <- olig.results[,c(1,3)]
  #print(small_dat[order(small_dat$padj),][c(1:10),])
  write.csv(olig.results[,c(1:7,9)], here("processed-data", "snRNA-seq", "06_NMF", paste0(nmf_pattern,"_reactome_results_sig.csv")))

  # plot the top 3 pathways
  nmf.terms = small_dat[order(small_dat$padj),][c(1:3),]$pathway
  print(nmf.terms)
  p <- plotGseaTable(x.h[nmf.terms],
                     nmf.stats, olig.results,
                     gseaParam = 0.5)
  pdf(here("plots", "snRNA-seq", "06_NMF", "GSEA", paste0(nmf_pattern,".pdf")))
  print(p)
  dev.off()
}

# L1CAM
olig.results.9 <- read.csv(here("processed-data", "snRNA-seq", "06_NMF", paste0("nmf9","_reactome_results_sig.csv")))
l1cam.genes = strsplit(filter(olig.results.9, pathway=="L1CAM interactions")$leadingEdge, split="/")[[1]]
l1cam.genes2 = mapIds(org.Hs.eg.db, keys=l1cam.genes, keytype="ENTREZID", column="SYMBOL", multiVals = "first")
rank.olig = (1+nrow(loads))-rank(loads[,"nmf9"])
top.olig = names(rank.olig)[rank.olig<=186]
intersect(top.olig, l1cam.genes2) #"ALCAM" "SHTN1" "DNM3"  "DLG1"  "ANK2"  "NCAM1"

olig.results.13 <- read.csv(here("processed-data", "snRNA-seq", "06_NMF", paste0("nmf13","_reactome_results_sig.csv")))
l1cam.genes = strsplit(filter(olig.results.13, pathway=="L1CAM interactions")$leadingEdge, split="/")[[1]]
l1cam.genes2 = mapIds(org.Hs.eg.db, keys=l1cam.genes, keytype="ENTREZID", column="SYMBOL", multiVals = "first")
rank.olig = (1+nrow(loads))-rank(loads[,"nmf13"])
top.olig = names(rank.olig)[rank.olig<=50]
intersect(top.olig, l1cam.genes2) #"SHTN1" "DNM3"  "ANK3"  "NCAM1"

olig.results.27 <- read.csv(here("processed-data", "snRNA-seq", "06_NMF", paste0("nmf27","_reactome_results_sig.csv")))
l1cam.genes = strsplit(filter(olig.results.27, pathway=="L1CAM interactions")$leadingEdge, split="/")[[1]]
l1cam.genes2 = mapIds(org.Hs.eg.db, keys=l1cam.genes, keytype="ENTREZID", column="SYMBOL", multiVals = "first")
rank.olig = (1+nrow(loads))-rank(loads[,"nmf27"])
top.olig = names(rank.olig)[rank.olig<=50]
intersect(top.olig, l1cam.genes2) #"SHTN1" "DNM3"  "ANK3"  "DLG1"

# Dev Bio
olig.results.9 <- read.csv(here("processed-data", "snRNA-seq", "06_NMF", paste0("nmf9","_reactome_results_sig.csv")))
l1cam.genes = strsplit(filter(olig.results.9, pathway=="Developmental Biology")$leadingEdge, split="/")[[1]]
l1cam.genes2 = mapIds(org.Hs.eg.db, keys=l1cam.genes, keytype="ENTREZID", column="SYMBOL", multiVals = "first")
rank.olig = (1+nrow(loads))-rank(loads[,"nmf9"])
top.olig = names(rank.olig)[rank.olig<=50]
intersect(top.olig, l1cam.genes2) #"MBP"   "UNC5C" "GAB1"  "DNM3"  "PCSK6"

olig.results.23 <- read.csv(here("processed-data", "snRNA-seq", "06_NMF", paste0("nmf23","_reactome_results_sig.csv")))
l1cam.genes = strsplit(filter(olig.results.23, pathway=="Developmental Biology")$leadingEdge, split="/")[[1]]
l1cam.genes2 = mapIds(org.Hs.eg.db, keys=l1cam.genes, keytype="ENTREZID", column="SYMBOL", multiVals = "first")
rank.olig = (1+nrow(loads))-rank(loads[,"nmf23"])
top.olig = names(rank.olig)[rank.olig<=50]
intersect(top.olig, l1cam.genes2) #"KAZN"    "PCSK6"   "DSCAML1" "POLR2F"  "SEMA4D"  "PKP4"    "DPYSL5" "NEO1"    "MYO9B"   "DNM2"

olig.results.26 <- read.csv(here("processed-data", "snRNA-seq", "06_NMF", paste0("nmf26","_reactome_results_sig.csv")))
l1cam.genes = strsplit(filter(olig.results.26, pathway=="Developmental Biology")$leadingEdge, split="/")[[1]]
l1cam.genes2 = mapIds(org.Hs.eg.db, keys=l1cam.genes, keytype="ENTREZID", column="SYMBOL", multiVals = "first")
rank.olig = (1+nrow(loads))-rank(loads[,"nmf26"])
top.olig = names(rank.olig)[rank.olig<=50]
intersect(top.olig, l1cam.genes2) #"UNC5C"  "ALCAM"  "COL4A5" "ITGA2"

olig.results.36 <- read.csv(here("processed-data", "snRNA-seq", "06_NMF", paste0("nmf36","_reactome_results_sig.csv")))
l1cam.genes = strsplit(filter(olig.results.36, pathway=="Developmental Biology")$leadingEdge, split="/")[[1]]
l1cam.genes2 = mapIds(org.Hs.eg.db, keys=l1cam.genes, keytype="ENTREZID", column="SYMBOL", multiVals = "first")
rank.olig = (1+nrow(loads))-rank(loads[,"nmf36"])
top.olig = names(rank.olig)[rank.olig<=50]
intersect(top.olig, l1cam.genes2) #"MBP"   "UNC5C" "TCF12" "ANK3"  "LAMA2" "PTK2"  "PKP4"

# Nervous System Development
olig.results.13 <- read.csv(here("processed-data", "snRNA-seq", "06_NMF", paste0("nmf13","_reactome_results_sig.csv")))
l1cam.genes = strsplit(filter(olig.results.13, pathway=="Nervous system development")$leadingEdge, split="/")[[1]]
l1cam.genes2 = mapIds(org.Hs.eg.db, keys=l1cam.genes, keytype="ENTREZID", column="SYMBOL", multiVals = "first")
rank.olig = (1+nrow(loads))-rank(loads[,"nmf13"])
top.olig = names(rank.olig)[rank.olig<=50]
intersect(top.olig, l1cam.genes2) #"MBP"   "UNC5C" "SHTN1" "DNM3"  "ANK3"  "NCAM1"

olig.results.40 <- read.csv(here("processed-data", "snRNA-seq", "06_NMF", paste0("nmf40","_reactome_results_sig.csv")))
l1cam.genes = strsplit(filter(olig.results.40, pathway=="Nervous system development")$leadingEdge, split="/")[[1]]
l1cam.genes2 = mapIds(org.Hs.eg.db, keys=l1cam.genes, keytype="ENTREZID", column="SYMBOL", multiVals = "first")
rank.olig = (1+nrow(loads))-rank(loads[,"nmf40"])
top.olig = names(rank.olig)[rank.olig<=50]
intersect(top.olig, l1cam.genes2) #"ANK3"  "DOCK1" "SCD5"

# Neuronal System
olig.results.28 <- read.csv(here("processed-data", "snRNA-seq", "06_NMF", paste0("nmf28","_reactome_results_sig.csv")))
l1cam.genes = strsplit(filter(olig.results.28, pathway=="Neuronal System")$leadingEdge, split="/")[[1]]
l1cam.genes2 = mapIds(org.Hs.eg.db, keys=l1cam.genes, keytype="ENTREZID", column="SYMBOL", multiVals = "first")
rank.olig = (1+nrow(loads))-rank(loads[,"nmf28"])
top.olig = names(rank.olig)[rank.olig<=50]
intersect(top.olig, l1cam.genes2) #"IL1RAPL1" "NRXN3"    "NLGN1"    "KCNH8"

olig.results.40 <- read.csv(here("processed-data", "snRNA-seq", "06_NMF", paste0("nmf40","_reactome_results_sig.csv")))
l1cam.genes = strsplit(filter(olig.results.40, pathway=="Neuronal System")$leadingEdge, split="/")[[1]]
l1cam.genes2 = mapIds(org.Hs.eg.db, keys=l1cam.genes, keytype="ENTREZID", column="SYMBOL", multiVals = "first")
rank.olig = (1+nrow(loads))-rank(loads[,"nmf40"])
top.olig = names(rank.olig)[rank.olig<=50]
intersect(top.olig, l1cam.genes2) #"ERBB4"  "NRXN3"  "NLGN1"  "PPFIA2" "GRIA2"

#RHO GTPase Cycle
olig.results.28 <- read.csv(here("processed-data", "snRNA-seq", "06_NMF", paste0("nmf28","_reactome_results_sig.csv")))
l1cam.genes = strsplit(filter(olig.results.28, pathway=="RHO GTPase cycle")$leadingEdge, split="/")[[1]]
l1cam.genes2 = mapIds(org.Hs.eg.db, keys=l1cam.genes, keytype="ENTREZID", column="SYMBOL", multiVals = "first")
rank.olig = (1+nrow(loads))-rank(loads[,"nmf28"])
top.olig = names(rank.olig)[rank.olig<=50]
intersect(top.olig, l1cam.genes2) #"DOCK10" "DLC1"

olig.results.36 <- read.csv(here("processed-data", "snRNA-seq", "06_NMF", paste0("nmf36","_reactome_results_sig.csv")))
l1cam.genes = strsplit(filter(olig.results.36, pathway=="RHO GTPase cycle")$leadingEdge, split="/")[[1]]
l1cam.genes2 = mapIds(org.Hs.eg.db, keys=l1cam.genes, keytype="ENTREZID", column="SYMBOL", multiVals = "first")
rank.olig = (1+nrow(loads))-rank(loads[,"nmf36"])
top.olig = names(rank.olig)[rank.olig<=50]
intersect(top.olig, l1cam.genes2) #"PKP4"

olig.results.39 <- read.csv(here("processed-data", "snRNA-seq", "06_NMF", paste0("nmf39","_reactome_results_sig.csv")))
l1cam.genes = strsplit(filter(olig.results.39, pathway=="RHO GTPase cycle")$leadingEdge, split="/")[[1]]
l1cam.genes2 = mapIds(org.Hs.eg.db, keys=l1cam.genes, keytype="ENTREZID", column="SYMBOL", multiVals = "first")
rank.olig = (1+nrow(loads))-rank(loads[,"nmf39"])
top.olig = names(rank.olig)[rank.olig<=50]
intersect(top.olig, l1cam.genes2) #"PLD1"     "PREX1"    "STARD13"  "PICALM"   "ADD3"     "ARHGAP32"

olig.results.33 <- read.csv(here("processed-data", "snRNA-seq", "06_NMF", paste0("nmf33","_reactome_results_sig.csv")))
l1cam.genes = strsplit(filter(olig.results.33, pathway=="Signaling by Rho GTPases")$leadingEdge, split="/")[[1]]
l1cam.genes2 = mapIds(org.Hs.eg.db, keys=l1cam.genes, keytype="ENTREZID", column="SYMBOL", multiVals = "first")
rank.olig = (1+nrow(loads))-rank(loads[,"nmf33"])
top.olig = names(rank.olig)[rank.olig<=50]
intersect(top.olig, l1cam.genes2) #"PREX1"    "CALM1"    "HSP90AA1" "TUBA1A"   "PPP1R14A" "ACTB"     "TUBB4A"

#CDC42 GTPase Cycle
olig.results.36 <- read.csv(here("processed-data", "snRNA-seq", "06_NMF", paste0("nmf36","_reactome_results_sig.csv")))
l1cam.genes = strsplit(filter(olig.results.36, pathway=="CDC42 GTPase cycle")$leadingEdge, split="/")[[1]]
l1cam.genes2 = mapIds(org.Hs.eg.db, keys=l1cam.genes, keytype="ENTREZID", column="SYMBOL", multiVals = "first")
rank.olig = (1+nrow(loads))-rank(loads[,"nmf36"])
top.olig = names(rank.olig)[rank.olig<=100]
intersect(top.olig, l1cam.genes2) #"DOCK10"  "FMNL2"   "FNBP1"   "STARD13"

olig.results.39 <- read.csv(here("processed-data", "snRNA-seq", "06_NMF", paste0("nmf39","_reactome_results_sig.csv")))
l1cam.genes = strsplit(filter(olig.results.39, pathway=="CDC42 GTPase cycle")$leadingEdge, split="/")[[1]]
l1cam.genes2 = mapIds(org.Hs.eg.db, keys=l1cam.genes, keytype="ENTREZID", column="SYMBOL", multiVals = "first")
rank.olig = (1+nrow(loads))-rank(loads[,"nmf39"])
top.olig = names(rank.olig)[rank.olig<=50]
intersect(top.olig, l1cam.genes2) #"PLD1"     "PREX1"    "STARD13"  "ARHGAP32"

# create heatmap of selected genes
select.nmfs = c("nmf33","nmf28","nmf40","nmf39","nmf43","nmf26","nmf36","nmf23","nmf27","nmf13","nmf9")
nmf.genes = c("SHTN1", "DNM3","ANK3", #L1CAM
              "UNC5C", "MBP", "PCSK6", "PKP4", #Dev Bio
              "NRXN3", "NLGN1", #Neuronal System
              "PREX1", #RHO GTPase
              "STARD13") #CDC42 GTPase

m1 = loads[nmf.genes,select.nmfs]
pdf(here("plots", "snRNA-seq", "06_NMF", "GSEA", paste0("GSEA_heatmap",".pdf")))
pheatmap::pheatmap(m1, scale="column", cluster_rows = F, cluster_cols = F,
                   angle_col=0)
dev.off()
