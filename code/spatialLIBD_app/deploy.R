library("rsconnect")

## Locate app_dir. Edit as needed
app_dir <- here::here("code", "spatialLIBD_app")

## Or you can go to your shinyapps.io account and copy this
## Here we do this to keep our information hidden.
#load(file.path(app_dir, ".deploy_info.Rdata"), verbose = TRUE)

source(file.path(app_dir, "token.R"))

## You need this to enable shinyapps to install Bioconductor packages
options(repos = BiocManager::repositories())

## Deploy the app, that is, upload it to shinyapps.io
## Note that appFiles has to be relative to app_dir.
## Drop the www directory if you didn't customize the documentation files and
## edit app.R accordingly.
rsconnect::deployApp(
    appDir = app_dir,
    appFiles = c(
        "app.R",
        "spe_nnSVG_PRECAST_9_labels.Rdata",
        "nnSVG_PRECAST_captureArea_9.Rdata",
        "nnSVG_PRECAST_captureArea_9_sig_genes_all.rds",
        "modeling-nnSVG_PRECAST_captureArea_9.Rdata"
    ),
    appName = "spatialdACC_Visium",
    account = "libd",
    server = "shinyapps.io"
)
