library("rsconnect")

# source("token.R")

options(repos = BiocManager::repositories())
rsconnect::deployApp(
    appFiles = c("app.R", "sce_for_iSEE_LS.rda", "initial.R"),
    appName = "snRNAseq_lateral_septum",
    account = "libd",
    server = "shinyapps.io"
)
