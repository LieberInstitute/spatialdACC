# https://github.com/jac-thom/findBiasedFeatures/blob/main/spatialLIBD.qmd

setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

library(SpatialExperiment)
library(dplyr)
library(ggplot2)
library(ggspavis)
library(gridExtra)
library(ggrepel)
library(here)
library(scran)
library(scry)

#start from spe without batch correction
load(here("processed-data", "06_preprocessing", "spe_dimred.Rdata"))

#remove samples with batch effect suspected
# "V12N28−332_A1" and "V12N28−332_B1"
suspected_batch_samples <- c("V12N28-332_A1", "V12N28-332_B1")
spe <- spe[, !colData(spe)$sample_id %in% suspected_batch_samples]

table(colData(spe)[,c("sample_id","brnum")])

mv <- modelGeneVarByPoisson(logcounts(spe))
mv$ensembl = rownames(mv)
mv$rank = (nrow(mv)+1)-rank(mv$bio)
top_hvgs <- getTopHVGs(mv, n = 3000)

bd <- devianceFeatureSelection(spe, fam="binomial")
bd.df = cbind.data.frame("gene"=rownames(bd),"gene_name"=rowData(bd)$gene_name,
                         "dev"= rowData(bd)$binomial_deviance,
                         "rank"=(nrow(bd)+1)-rank(rowData(bd)$binomial_deviance))
rownames(bd.df) = bd.df$gene

load(file=here::here('processed-data', '08_clustering', 'batch_effect','nnSVG_1000_samples_removed.rda'))
genes <- rownames(df_summaryReplicated)

mv$is_svg = factor(mv$ensembl %in% genes,
                   levels=c(TRUE,FALSE), labels=c("SVGs","not SVGs"))
bd.df$is_svg = factor(bd.df$gene %in% genes,
                      levels=c(TRUE,FALSE), labels=c("SVGs","not SVGs"))

var1 <- ggplot(mv, aes(x=rank, y=bio, color=is_svg))+
    geom_point(size=.5)+scale_color_manual(values=c("tomato","black"))+
    labs(title="mean-variance",
         y="variance", color="")+
    theme_bw()+theme(legend.position="inside", legend.position.inside=c(.8,.8))
dev1 <- ggplot(bd.df, aes(x=rank, y=dev, color=is_svg))+
    geom_point(size=.5)+scale_color_manual(values=c("tomato","black"))+
    labs(title="binomial deviance",
         y="deviance", color="")+
    theme_bw()+theme(legend.position="inside", legend.position.inside=c(.8,.8),
                     axis.text.y=element_text(size=6))
pdf(here("plots", '08_clustering', 'batch_effect', "mean-variance_binomial_deviance_models.pdf"), width=8, height=4)
grid.arrange(var1, dev1, ncol=2)
dev.off()

# examine for brnum first
bd.batch <- devianceFeatureSelection(spe, fam="binomial",
                                     batch=as.factor(spe$brnum))
bd.batch.df = cbind.data.frame("gene"=rownames(bd.batch),"gene_name"=rowData(bd.batch)$gene_name,
                               "dev"= rowData(bd.batch)$binomial_deviance,
                               "rank"=(nrow(bd.batch)+1)-rank(rowData(bd.batch)$binomial_deviance))
rownames(bd.batch.df) = bd.batch.df$gene
bd.batch.df$is_svg = factor(bd.batch.df$gene %in% genes,
                            levels=c(TRUE,FALSE), labels=c("SVGs","not SVGs"))
dev3 <- ggplot(bd.batch.df, aes(x=rank, y=dev, color=is_svg))+
    geom_point(size=.5)+scale_color_manual(values=c("red","black"))+
    labs(title="batch = subject",
         y="deviance", color="")+
    theme_bw()+theme(legend.position="inside", legend.position.inside=c(.8,.8),
                     axis.text.y=element_text(size=6))

subject.df <- left_join(bd.df, bd.batch.df,
                        by=c("gene", "gene_name"),
                        suffix=c("_default","_subject"))
subject.df$d.diff = (subject.df$dev_default-subject.df$dev_subject)/subject.df$dev_subject
subject.df$r.diff = subject.df$rank_subject-subject.df$rank_default
top.df <- filter(subject.df, rank_default<=3000, rank_subject<=3000)

delta1 <- ggplot(top.df, aes(x=dev_default, y=dev_subject, color=d.diff))+
    geom_point()+
    geom_text_repel(data=filter(top.df, d.diff>1), aes(label=gene_name))+
    scale_x_log10()+scale_y_log10()+
    scale_color_viridis_c(option="F", direction=-1)+
    geom_abline(aes(slope=1, intercept=0), lty=2)+
    labs(x="deviance (no batch)", y="deviance (batch)", title="\u0394 deviance")+
    theme_bw()+theme(legend.position="inside", legend.position.inside=c(.1,.8))+
    theme(legend.title = element_text(size = 10),
          legend.text  = element_text(size = 10),
          legend.key.size = unit(.6, "lines"))
delta2 <- ggplot(top.df, aes(x=rank_default, y=rank_subject, color=r.diff))+
    geom_point()+scale_y_reverse()+
    geom_text_repel(data=filter(top.df, r.diff>1000), aes(label=gene_name))+
    scale_color_viridis_c(option="F", direction=-1)+
    geom_abline(aes(slope=-1, intercept=0), lty=2)+
    labs(x="rank (no batch)",y="rank (batch)", title="\u0394 rank")+
    theme_bw()+theme(legend.position="inside", legend.position.inside=c(.8,.8))+
    theme(legend.title = element_text(size = 10),
          legend.text  = element_text(size = 10),
          legend.key.size = unit(.6, "lines"))

top.svg.df = filter(subject.df, gene %in% genes)
delta3 <- ggplot(top.svg.df, aes(x=dev_default, y=dev_subject, color=d.diff))+
    geom_point()+
    geom_text_repel(data=filter(top.svg.df, d.diff>1),# & is_svg=="SVGs"),
                    aes(label=gene_name), color="black")+
    scale_x_log10()+scale_y_log10()+
    scale_color_viridis_c(option="F", direction=-1)+
    geom_abline(aes(slope=1, intercept=0), lty=2)+
    labs(x="deviance (no batch)", y="deviance (batch)", title="\u0394 deviance")+
    theme_bw()+
    theme(legend.position="none")
delta4 <- ggplot(top.svg.df, aes(x=rank_default, y=rank_subject, color=r.diff))+
    geom_point()+scale_y_reverse()+xlim(0,6000)+
    geom_text_repel(data=filter(top.svg.df, r.diff>1000),# & is_svg=="SVGs"),
                    aes(label=gene_name), color="black")+
    scale_color_viridis_c(option="F", direction=-1)+
    geom_abline(aes(slope=-1, intercept=0), lty=2)+
    labs(x="rank (no batch)",y="rank (batch)", title="\u0394 rank")+
    theme_bw()+
    theme(legend.position="none")

top.df$is_mito = top.df$gene_name %in% grep("^MT-", top.df$gene_name, value=T)
top.df$is_biased = top.df$d.diff>1 | top.df$r.diff>1000
top.df$other = top.df$is_mito==F & top.df$is_biased==F & top.df$d.diff>.1
top.df$group = factor(paste(top.df$is_mito, top.df$is_biased),
                      levels=c("FALSE FALSE","FALSE TRUE","TRUE FALSE"))

mt <- ggplot(top.df, aes(x=r.diff, y=d.diff, color=group))+
    geom_point()+
    geom_text_repel(data= filter(top.df, d.diff>.21 | r.diff>500),
                    aes(label=gene_name), size=3)+
    ggbreak::scale_y_break(c(.5, 3), scales=.25, ticklabels=c(seq(0,.5,.1),3.0,3.1,3.2))+
    scale_color_manual(values=c("grey40","black","red"))+
    labs(x="\u0394 rank", y="\u0394 deviance", title="chrM in red; biased genes black")+
    theme_bw()+theme(legend.position="none", title=element_text(size=10))

ml1 <- plotSpots(spe, annotate="ENSG00000256618", assay="logcounts",
                 sample_id="sample_id", point_size=.1)+ggtitle("MTRNR2L1")+
    scale_color_gradient(low='grey90', high='black')
ml8 <- plotSpots(spe, annotate="ENSG00000255823", assay="logcounts",
                 sample_id="sample_id", point_size=.1)+ggtitle("MTRNR2L8")+
    scale_color_gradient(low='grey90', high='black')


pdf(here("plots", "08_clustering", "batch_effect", "subject_binomial_deviance_models.pdf"), width=18, height=14)
grid.arrange(dev1+ggtitle("batch = NULL"), dev3, ncol=2)
grid.arrange(delta1, delta2, ncol=2)
grid.arrange(delta3, delta4, ncol=2)
grid.arrange(mt, ncol=1)
grid.arrange(ml1, ml8, ncol=2)
dev.off()

# examine for slide next
bd.batch <- devianceFeatureSelection(spe, fam="binomial",
                                     batch=as.factor(spe$slide))
bd.batch.df = cbind.data.frame("gene"=rownames(bd.batch),"gene_name"=rowData(bd.batch)$gene_name,
                               "dev"= rowData(bd.batch)$binomial_deviance,
                               "rank"=(nrow(bd.batch)+1)-rank(rowData(bd.batch)$binomial_deviance))
rownames(bd.batch.df) = bd.batch.df$gene
bd.batch.df$is_svg = factor(bd.batch.df$gene %in% genes,
                            levels=c(TRUE,FALSE), labels=c("SVGs","not SVGs"))
dev3 <- ggplot(bd.batch.df, aes(x=rank, y=dev, color=is_svg))+
    geom_point(size=.5)+scale_color_manual(values=c("red","black"))+
    labs(title="batch = slide",
         y="deviance", color="")+
    theme_bw()+theme(legend.position="inside", legend.position.inside=c(.8,.8),
                     axis.text.y=element_text(size=6))

slide.df <- left_join(bd.df, bd.batch.df,
                        by=c("gene", "gene_name"),
                        suffix=c("_default","_slide"))
slide.df$d.diff = (slide.df$dev_default-slide.df$dev_slide)/slide.df$dev_slide
slide.df$r.diff = slide.df$rank_slide-slide.df$rank_default
top.df <- filter(slide.df, rank_default<=3000, rank_slide<=3000)

delta1 <- ggplot(top.df, aes(x=dev_default, y=dev_slide, color=d.diff))+
    geom_point()+
    geom_text_repel(data=filter(top.df, d.diff>1), aes(label=gene_name))+
    scale_x_log10()+scale_y_log10()+
    scale_color_viridis_c(option="F", direction=-1)+
    geom_abline(aes(slope=1, intercept=0), lty=2)+
    labs(x="deviance (no batch)", y="deviance (batch)", title="\u0394 deviance")+
    theme_bw()+theme(legend.position="inside", legend.position.inside=c(.1,.8))+
    theme(legend.title = element_text(size = 10),
          legend.text  = element_text(size = 10),
          legend.key.size = unit(.6, "lines"))
delta2 <- ggplot(top.df, aes(x=rank_default, y=rank_slide, color=r.diff))+
    geom_point()+scale_y_reverse()+
    geom_text_repel(data=filter(top.df, r.diff>1000), aes(label=gene_name))+
    scale_color_viridis_c(option="F", direction=-1)+
    geom_abline(aes(slope=-1, intercept=0), lty=2)+
    labs(x="rank (no batch)",y="rank (batch)", title="\u0394 rank")+
    theme_bw()+theme(legend.position="inside", legend.position.inside=c(.8,.8))+
    theme(legend.title = element_text(size = 10),
          legend.text  = element_text(size = 10),
          legend.key.size = unit(.6, "lines"))

top.svg.df = filter(slide.df, gene %in% genes)
delta3 <- ggplot(top.svg.df, aes(x=dev_default, y=dev_slide, color=d.diff))+
    geom_point()+
    geom_text_repel(data=filter(top.svg.df, d.diff>1),# & is_svg=="SVGs"),
                    aes(label=gene_name), color="black")+
    scale_x_log10()+scale_y_log10()+
    scale_color_viridis_c(option="F", direction=-1)+
    geom_abline(aes(slope=1, intercept=0), lty=2)+
    labs(x="deviance (no batch)", y="deviance (batch)", title="\u0394 deviance")+
    theme_bw()+
    theme(legend.position="none")
delta4 <- ggplot(top.svg.df, aes(x=rank_default, y=rank_slide, color=r.diff))+
    geom_point()+scale_y_reverse()+xlim(0,6000)+
    geom_text_repel(data=filter(top.svg.df, r.diff>1000),# & is_svg=="SVGs"),
                    aes(label=gene_name), color="black")+
    scale_color_viridis_c(option="F", direction=-1)+
    geom_abline(aes(slope=-1, intercept=0), lty=2)+
    labs(x="rank (no batch)",y="rank (batch)", title="\u0394 rank")+
    theme_bw()+
    theme(legend.position="none")

top.df$is_mito = top.df$gene_name %in% grep("^MT-", top.df$gene_name, value=T)
top.df$is_biased = top.df$d.diff>1 | top.df$r.diff>1000
top.df$other = top.df$is_mito==F & top.df$is_biased==F & top.df$d.diff>.1
top.df$group = factor(paste(top.df$is_mito, top.df$is_biased),
                      levels=c("FALSE FALSE","FALSE TRUE","TRUE FALSE"))

mt <- ggplot(top.df, aes(x=r.diff, y=d.diff, color=group))+
    geom_point()+
    geom_text_repel(data= filter(top.df, d.diff>.21 | r.diff>500),
                    aes(label=gene_name), size=3)+
    ggbreak::scale_y_break(c(.5, 3), scales=.25, ticklabels=c(seq(0,.5,.1),3.0,3.1,3.2))+
    scale_color_manual(values=c("grey40","black","red"))+
    labs(x="\u0394 rank", y="\u0394 deviance", title="chrM in red; biased genes black")+
    theme_bw()+theme(legend.position="none", title=element_text(size=10))

ml1 <- plotSpots(spe, annotate="ENSG00000198840", assay="logcounts",
                 sample_id="sample_id", point_size=.1)+ggtitle("MT−ND3")+
    scale_color_gradient(low='grey90', high='black')


pdf(here("plots", "08_clustering", "batch_effect", "slide_binomial_deviance_models.pdf"), width=9, height=14)
grid.arrange(dev1+ggtitle("batch = NULL"), dev3, ncol=2)
grid.arrange(delta1, delta2, ncol=2)
grid.arrange(delta3, delta4, ncol=2)
grid.arrange(mt, ncol=1)
grid.arrange(ml1, ncol=1)
dev.off()

# examine for sample_id
bd.batch <- devianceFeatureSelection(spe, fam="binomial",
                                     batch=as.factor(spe$sample_id))
bd.batch.df = cbind.data.frame("gene"=rownames(bd.batch),"gene_name"=rowData(bd.batch)$gene_name,
                               "dev"= rowData(bd.batch)$binomial_deviance,
                               "rank"=(nrow(bd.batch)+1)-rank(rowData(bd.batch)$binomial_deviance))
rownames(bd.batch.df) = bd.batch.df$gene
bd.batch.df$is_svg = factor(bd.batch.df$gene %in% genes,
                            levels=c(TRUE,FALSE), labels=c("SVGs","not SVGs"))
dev3 <- ggplot(bd.batch.df, aes(x=rank, y=dev, color=is_svg))+
    geom_point(size=.5)+scale_color_manual(values=c("red","black"))+
    labs(title="batch = sample",
         y="deviance", color="")+
    theme_bw()+theme(legend.position="inside", legend.position.inside=c(.8,.8),
                     axis.text.y=element_text(size=6))

sample.df <- left_join(bd.df, bd.batch.df,
                        by=c("gene", "gene_name"),
                        suffix=c("_default","_sample"))
sample.df$d.diff = (sample.df$dev_default-sample.df$dev_sample)/sample.df$dev_sample
sample.df$r.diff = sample.df$rank_sample-sample.df$rank_default
top.df <- filter(sample.df, rank_default<=3000, rank_sample<=3000)

delta1 <- ggplot(top.df, aes(x=dev_default, y=dev_sample, color=d.diff))+
    geom_point()+
    geom_text_repel(data=filter(top.df, d.diff>1), aes(label=gene_name))+
    scale_x_log10()+scale_y_log10()+
    scale_color_viridis_c(option="F", direction=-1)+
    geom_abline(aes(slope=1, intercept=0), lty=2)+
    labs(x="deviance (no batch)", y="deviance (batch)", title="\u0394 deviance")+
    theme_bw()+theme(legend.position="inside", legend.position.inside=c(.1,.8))+
    theme(legend.title = element_text(size = 10),
          legend.text  = element_text(size = 10),
          legend.key.size = unit(.6, "lines"))
delta2 <- ggplot(top.df, aes(x=rank_default, y=rank_sample, color=r.diff))+
    geom_point()+scale_y_reverse()+
    geom_text_repel(data=filter(top.df, r.diff>1000), aes(label=gene_name))+
    scale_color_viridis_c(option="F", direction=-1)+
    geom_abline(aes(slope=-1, intercept=0), lty=2)+
    labs(x="rank (no batch)",y="rank (batch)", title="\u0394 rank")+
    theme_bw()+theme(legend.position="inside", legend.position.inside=c(.8,.8))+
    theme(legend.title = element_text(size = 10),
          legend.text  = element_text(size = 10),
          legend.key.size = unit(.6, "lines"))

top.svg.df = filter(sample.df, gene %in% genes)
delta3 <- ggplot(top.svg.df, aes(x=dev_default, y=dev_sample, color=d.diff))+
    geom_point()+
    geom_text_repel(data=filter(top.svg.df, d.diff>1),# & is_svg=="SVGs"),
                    aes(label=gene_name), color="black")+
    scale_x_log10()+scale_y_log10()+
    scale_color_viridis_c(option="F", direction=-1)+
    geom_abline(aes(slope=1, intercept=0), lty=2)+
    labs(x="deviance (no batch)", y="deviance (batch)", title="\u0394 deviance")+
    theme_bw()+
    theme(legend.position="none")
delta4 <- ggplot(top.svg.df, aes(x=rank_default, y=rank_sample, color=r.diff))+
    geom_point()+scale_y_reverse()+xlim(0,6000)+
    geom_text_repel(data=filter(top.svg.df, r.diff>1000),# & is_svg=="SVGs"),
                    aes(label=gene_name), color="black")+
    scale_color_viridis_c(option="F", direction=-1)+
    geom_abline(aes(slope=-1, intercept=0), lty=2)+
    labs(x="rank (no batch)",y="rank (batch)", title="\u0394 rank")+
    theme_bw()+
    theme(legend.position="none")

top.df$is_mito = top.df$gene_name %in% grep("^MT-", top.df$gene_name, value=T)
top.df$is_biased = top.df$d.diff>1 | top.df$r.diff>1000
top.df$other = top.df$is_mito==F & top.df$is_biased==F & top.df$d.diff>.1
top.df$group = factor(paste(top.df$is_mito, top.df$is_biased),
                      levels=c("FALSE FALSE","FALSE TRUE","TRUE FALSE"))

mt <- ggplot(top.df, aes(x=r.diff, y=d.diff, color=group))+
    geom_point()+
    geom_text_repel(data= filter(top.df, d.diff>.21 | r.diff>500),
                    aes(label=gene_name), size=3)+
    ggbreak::scale_y_break(c(.5, 3), scales=.25, ticklabels=c(seq(0,.5,.1),3.0,3.1,3.2))+
    scale_color_manual(values=c("grey40","black","red"))+
    labs(x="\u0394 rank", y="\u0394 deviance", title="chrM in red; biased genes black")+
    theme_bw()+theme(legend.position="none", title=element_text(size=10))

# which(grepl("MT[-−–]ND4", rowData(spe)$gene_name))

ml1 <- plotSpots(spe, annotate="ENSG00000256618", assay="logcounts",
                 sample_id="sample_id", point_size=0.05)+ggtitle("MTRNR2L1")+
    scale_color_gradient(low='grey90', high='black')
ml2 <- plotSpots(spe, annotate="ENSG00000255823", assay="logcounts",
                 sample_id="sample_id", point_size=0.05)+ggtitle("MTRNR2L8")+
    scale_color_gradient(low='grey90', high='black')
ml3 <- plotSpots(spe, annotate="ENSG00000198840", assay="logcounts",
                 sample_id="sample_id", point_size=0.05)+ggtitle("MT-ND3")+
    scale_color_gradient(low='grey90', high='black')
ml4 <- plotSpots(spe, annotate="ENSG00000198886", assay="logcounts",
                 sample_id="sample_id", point_size=0.05)+ggtitle("MT-ND4")+
    scale_color_gradient(low='grey90', high='black')
ml5 <- plotSpots(spe, annotate="ENSG00000198899", assay="logcounts",
                 sample_id="sample_id", point_size=0.05)+ggtitle("MT-ATP6")+
    scale_color_gradient(low='grey90', high='black')
ml6 <- plotSpots(spe, annotate="ENSG00000198804", assay="logcounts",
                 sample_id="sample_id", point_size=0.05)+ggtitle("MT-CO1")+
    scale_color_gradient(low='grey90', high='black')
ml7 <- plotSpots(spe, annotate="ENSG00000198938", assay="logcounts",
                 sample_id="sample_id", point_size=0.05)+ggtitle("MT-CO3")+
    scale_color_gradient(low='grey90', high='black')

pdf(here("plots", "08_clustering", "batch_effect", "sample_binomial_deviance_models.pdf"), width=18, height=14)
grid.arrange(dev1+ggtitle("batch = NULL"), dev3, ncol=2)
grid.arrange(delta1, delta2, ncol=2)
grid.arrange(delta3, delta4, ncol=2)
grid.arrange(mt, ncol=1)
grid.arrange(ml1, ncol=1)
grid.arrange(ml2, ncol=1)
grid.arrange(ml3, ncol=1)
grid.arrange(ml4, ncol=1)
grid.arrange(ml5, ncol=1)
grid.arrange(ml6, ncol=1)
grid.arrange(ml7, ncol=1)
dev.off()
