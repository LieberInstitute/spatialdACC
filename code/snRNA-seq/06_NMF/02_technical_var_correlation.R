setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

####code for correlating categorical technical variables with NMF patterns
library(here)
library(SpatialExperiment)
library(pheatmap)
library(reshape2)

load(file = here("processed-data", "snRNA-seq", "03_batch_correction", "sce_harmony.Rdata"))
x <- readRDS(file = here("processed-data", "snRNA-seq", "06_NMF", "nmf_results.RDS"))

####onehot encode both brain as an example
data<-as.data.frame(sce$brain)
colnames(data)<-'brain'
onehot_brain <-  dcast(data = data, rownames(data) ~ brain, length)
rownames(onehot_brain)<-onehot_brain[,1]
onehot_brain[,1]<-as.numeric(onehot_brain[,1])
onehot_brain<-onehot_brain[order(onehot_brain[,1],decreasing=F),]
onehot_brain[,1]<-NULL

###correlate with nmf patterns
pheatmap(cor(t(x@h),onehot_brain))

###repeat with any categorical vars

###for continuous tech vars do this:

vars<-colData(sce)[,c('detected','sum','subsets_Mito_percent')]
vars$detected<-as.numeric(vars$detected)
pheatmap(cor(t(x@h),vars))
