setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')
suppressPackageStartupMessages({
    library("dplyr")
    library("here")
    library("sessioninfo")
    library("ggplot2")
})

df_bayesspace_harmony <- read.table(file=here("processed-data", "08_clustering", "cluster_diagnostics", "fasthplus", "fasthplus_results_bayesSpace_harmony.csv"), header=T)
#problem with writing to csv file, added header every time
row_odd <- seq_len(nrow(df_bayesspace_harmony)) %% 2
df_bayesspace_harmony <- df_bayesspace_harmony[row_odd == 1, ]
df_bayesspace_harmony$k<-as.numeric(df_bayesspace_harmony$k)
df_bayesspace_harmony$fasthplus<-as.numeric(df_bayesspace_harmony$fasthplus)

df_bayesspace_mnn <- read.table(file=here("processed-data", "08_clustering", "cluster_diagnostics", "fasthplus", "fasthplus_results_bayesSpace_mnn.csv"), header=T)
#problem with writing to csv file, added header every time
row_odd <- seq_len(nrow(df_bayesspace_mnn)) %% 2
df_bayesspace_mnn <- df_bayesspace_mnn[row_odd == 1, ]
df_bayesspace_mnn$k<-as.numeric(df_bayesspace_mnn$k)
df_bayesspace_mnn$fasthplus<-as.numeric(df_bayesspace_mnn$fasthplus)

df_precast <- read.table(file=here("processed-data", "08_clustering", "cluster_diagnostics", "fasthplus", "fasthplus_results_precast.csv"), header=T)
#problem with writing to csv file, added header every time
row_odd <- seq_len(nrow(df_precast)) %% 2
df_precast <- df_precast[row_odd == 1, ]
df_precast$k<-as.numeric(df_precast$k)
df_precast$fasthplus<-as.numeric(df_precast$fasthplus)

df_nnSVG_precast <- read.table(file=here("processed-data", "08_clustering", "cluster_diagnostics", "fasthplus", "fasthplus_results_nnSVG_precast.csv"), header=T)

pdf(here("plots", "08_clustering", "fasthplus.pdf"), width = 21)


ggplot(data = df_bayesspace_harmony, aes(x = k, y = fasthplus, group = 1)) +
    geom_line() +
    geom_point() +
    ylab(expression(H^{
        "+"
    })) +
    theme_bw(base_size = 20) +
    ggtitle("Bayes Space post Harmony")

ggplot(data = df_bayesspace_mnn, aes(x = k, y = fasthplus, group = 1)) +
    geom_line() +
    geom_point() +
    ylab(expression(H^{
        "+"
    })) +
    theme_bw(base_size = 20) +
    ggtitle("Bayes Space post MNN")

ggplot(data = df_precast, aes(x = k, y = fasthplus, group = 1)) +
    geom_line() +
    geom_point() +
    ylab(expression(H^{
        "+"
    })) +
    theme_bw(base_size = 20) +
    ggtitle("PRECAST")

ggplot(data = df_nnSVG_precast, aes(x = k, y = fasthplus, group = 1)) +
    geom_line() +
    geom_point() +
    ylab(expression(H^{
        "+"
    })) +
    theme_bw(base_size = 20) +
    ggtitle("nnSVG PRECAST")

dev.off()
