library(ggplot2)

# set up dataframes
filenames <- list.files(".", pattern="*.csv", full.names=TRUE)

Br6432_dACC <- read.csv(filenames[1])
Br6432_dlPFC <- read.csv(filenames[2])
Br8325_dlPFC <- read.csv(filenames[3])
Br8325_left_dACC <- read.csv(filenames[4])
Br8325_right_dACC <- read.csv(filenames[5])

cols <- c("X20x.Opal.570.Copies",
          "X20x.Opal.620.Copies", "X20x.Opal.690.Copies")
Br6432_dACC <- Br6432_dACC[,cols]
colnames(Br6432_dACC) <- c("POU3F1", "SULF2", "GABRQ")

Br6432_dlPFC <- Br6432_dlPFC[,cols]
colnames(Br6432_dlPFC) <- c("POU3F1", "SULF2", "GABRQ")

Br8325_dlPFC <- Br8325_dlPFC[,cols]
colnames(Br8325_dlPFC) <- c("POU3F1", "SULF2", "GABRQ")

Br8325_left_dACC <- Br8325_left_dACC[,cols]
colnames(Br8325_left_dACC) <- c("POU3F1", "SULF2", "GABRQ")

Br8325_right_dACC <- Br8325_right_dACC[,cols]
colnames(Br8325_right_dACC) <- c("POU3F1", "SULF2", "GABRQ")

# compute correlations
# make df of correlations
corr_list <- c(cor(Br6432_dACC$POU3F1, Br6432_dACC$SULF2),
               cor(Br6432_dACC$POU3F1, Br6432_dACC$GABRQ),
               cor(Br6432_dACC$SULF2, Br6432_dACC$GABRQ),
               cor(Br6432_dlPFC$POU3F1, Br6432_dlPFC$SULF2),
               cor(Br6432_dlPFC$POU3F1, Br6432_dlPFC$GABRQ),
               cor(Br6432_dlPFC$SULF2, Br6432_dlPFC$GABRQ),
               cor(Br8325_dlPFC$POU3F1, Br8325_dlPFC$SULF2),
               cor(Br8325_dlPFC$POU3F1, Br8325_dlPFC$GABRQ),
               cor(Br8325_dlPFC$SULF2, Br8325_dlPFC$GABRQ),
               cor(Br8325_left_dACC$POU3F1, Br8325_left_dACC$SULF2),
               cor(Br8325_left_dACC$POU3F1, Br8325_left_dACC$GABRQ),
               cor(Br8325_left_dACC$SULF2, Br8325_left_dACC$GABRQ),
               cor(Br8325_right_dACC$POU3F1, Br8325_right_dACC$SULF2),
               cor(Br8325_right_dACC$POU3F1, Br8325_right_dACC$GABRQ),
               cor(Br8325_right_dACC$SULF2, Br8325_right_dACC$GABRQ))

region_list <- c(rep("dACC",3), rep("dlPFC",3), rep("dlPFC",3),
                 rep("dACC",3), rep("dACC",3))

genes_list <- c(rep(c("POU3F1 & SULF2","POU3F1 & GABRQ","GABRQ & SULF2"),5))

df <- data.frame(
    Correlation = corr_list,
    Region = region_list,
    Genes = genes_list
)

ggplot(df, aes(x=Genes, y=Correlation, color=Region)) +
    geom_point() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust=1))

# compute distributions of counts of 3 genes

# put these into sce and then use library normalization to hopefully integrate

# make PCA plots and color by expr of each gene, region, donor?

# identify VENs in each region
