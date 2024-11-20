setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

####code for correlating categorical technical variables with NMF patterns
library(here)
library(SpatialExperiment)
library(pheatmap)
library(reshape2)

load(file = here("processed-data", "snRNA-seq", "03_batch_correction", "sce_harmony.Rdata"))
x <- readRDS(file = here("processed-data", "snRNA-seq", "06_NMF", "nmf_results.RDS"))

####onehot encode brain
data<-as.data.frame(sce$brain)
colnames(data)<-'brain'
onehot_brain <-  dcast(data = data, rownames(data) ~ brain, length)
rownames(onehot_brain)<-onehot_brain[,1]
onehot_brain[,1]<-as.numeric(onehot_brain[,1])
onehot_brain<-onehot_brain[order(onehot_brain[,1],decreasing=F),]
onehot_brain[,1]<-NULL

###correlate with nmf patterns
pdf(here("plots", "snRNA-seq", "06_NMF", "nmf_brain_correlation_heatmap.pdf"))
pheatmap(cor(t(x@h),onehot_brain), fontsize_row = 5)
dev.off()

# add sex to sce for each donor
colData(sce)$sex[colData(sce)$brain=="Br2743"] <- "M"
colData(sce)$sex[colData(sce)$brain=="Br2720"] <- "F"
colData(sce)$sex[colData(sce)$brain=="Br6432"] <- "M"
colData(sce)$sex[colData(sce)$brain=="Br6471"] <- "M"
colData(sce)$sex[colData(sce)$brain=="Br6522"] <- "M"
colData(sce)$sex[colData(sce)$brain=="Br3942"] <- "M"
colData(sce)$sex[colData(sce)$brain=="Br6423"] <- "M"
colData(sce)$sex[colData(sce)$brain=="Br8325"] <- "F"
colData(sce)$sex[colData(sce)$brain=="Br8492"] <- "F"
colData(sce)$sex[colData(sce)$brain=="Br8667"] <- "F"

####onehot encode sex
data<-as.data.frame(sce$sex)
colnames(data)<-'sex'
onehot_sex <-  dcast(data = data, rownames(data) ~ sex, length)
rownames(onehot_sex)<-onehot_sex[,1]
onehot_sex[,1]<-as.numeric(onehot_sex[,1])
onehot_sex<-onehot_sex[order(onehot_sex[,1],decreasing=F),]
onehot_sex[,1]<-NULL

###correlate with nmf patterns
pdf(here("plots", "snRNA-seq", "06_NMF", "nmf_sex_correlation_heatmap.pdf"))
pheatmap(cor(t(x@h),onehot_sex), fontsize_row = 5)
dev.off()

###for continuous tech vars do this:

vars<-colData(sce)[,c('detected','sum','subsets_Mito_percent')]

# fix error cor(t(x@h), vars) : 'y' must be numeric
vars$detected<-as.numeric(vars$detected)
vars$sum<-as.numeric(vars$sum)
vars$mito_percent<-as.numeric(vars$subsets_Mito_percent)
vars$subsets_Mito_percent<-NULL
vars <- as.matrix(vars)

pdf(here("plots", "snRNA-seq", "06_NMF", "nmf_qc_correlation_heatmap.pdf"))
pheatmap(cor(t(x@h),vars), fontsize_row = 5)
dev.off()

head(colData(spe)$age[colData(spe)$brnum=="Br8667"])

# add age to sce for each donor
colData(sce)$age[colData(sce)$brain=="Br2743"] <- "61.54"
colData(sce)$age[colData(sce)$brain=="Br2720"] <- "48.23"
colData(sce)$age[colData(sce)$brain=="Br6432"] <- "48.88"
colData(sce)$age[colData(sce)$brain=="Br6471"] <- "55.46"
colData(sce)$age[colData(sce)$brain=="Br6522"] <- "33.39"
colData(sce)$age[colData(sce)$brain=="Br3942"] <- "47.53"
colData(sce)$age[colData(sce)$brain=="Br6423"] <- "51.73"
colData(sce)$age[colData(sce)$brain=="Br8325"] <- "57.62"
colData(sce)$age[colData(sce)$brain=="Br8492"] <- "53.4"
colData(sce)$age[colData(sce)$brain=="Br8667"] <- "37.33"

vars<-colData(sce)[,c('age','scDblFinder.score')]
vars$age<-as.numeric(vars$age)
vars$scDblFinder.score<-as.numeric(vars$scDblFinder.score)
vars <- as.matrix(vars)

pdf(here("plots", "snRNA-seq", "06_NMF", "nmf_age_correlation_heatmap.pdf"))
pheatmap(cor(t(x@h),vars), fontsize_row = 5)
dev.off()
