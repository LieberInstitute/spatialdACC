#extract sample names from the master excel sheet

setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')
library("here")
library("readxl")
test <- as.data.frame(read_excel(file.path(here::here("raw-data", "sample_info", "dACC_Visium_summary_20230227_Maddy.xlsx"))))
colnames(test)
# [1] "Sample #"              "Tissue"                "Brain"                
# [4] "Slide #"               "Array #"               "Ct"                   
# [7] "cDNA Amp Cycle"        "Agilent [cDNA] pg/ul"  "Dilution Factor...9"  
# [10] "Final [cDNA] pg/ul"    "Total cDNA ng yield"   "cDNA Input"           
# [13] "SI cycles"             "Ave frag length"       "Agilent [lib] pg/ul"  
# [16] "Dilution Factor...16"  "Final [lib] pg/ul"     "index_name"           
# [19] "index(i7)"             "index2_workflow_a(i5)" "index2_workflow_b(i5)"
# [22] "% Coverage Array"      "Est Read Pairs"        "Replicate?"           
# [25] "Replicate number"      "H&E?"                  "10x names"            
# [28] "spaceranger output"    "images"               

test[,28][1:12]
# [1] "/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2023-03-29_Transfer10x_SPage/Transfer_to_LIBD_20230323/1444617/outs/"
# [2] "/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2023-02-09_Transfer10x_SPage/LIBD/libd_1/outs/"                      
# [3] "/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2023-03-29_Transfer10x_SPage/Transfer_to_LIBD_20230323/1444618/outs/"
# [4] "/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2023-03-29_Transfer10x_SPage/Transfer_to_LIBD_20230323/1444619/outs/"
# [5] "/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2023-03-29_Transfer10x_SPage/Transfer_to_LIBD_20230323/1444620/outs/"
# [6] "/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2023-03-29_Transfer10x_SPage/Transfer_to_LIBD_20230323/1444621/outs/"
# [7] "/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2023-03-29_Transfer10x_SPage/Transfer_to_LIBD_20230323/1444622/outs/"
# [8] "/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2023-02-09_Transfer10x_SPage/LIBD/libd_3/outs/"                      
# [9] "/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2023-02-09_Transfer10x_SPage/LIBD/libd_4/outs/"                      
# [10] "/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2023-02-09_Transfer10x_SPage/LIBD/libd_5/outs/"                      
# [11] "/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2023-02-09_Transfer10x_SPage/LIBD/libd_6/outs/"                      
# [12] "/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2023-02-09_Transfer10x_SPage/LIBD/libd_7/outs/"                      

samples = paste0(test[,4],"_",test[,5])[1:12]
samples = as.data.frame(samples)
samples$path = test[,28][1:12]
samples$directory = sapply(strsplit(samples$path,"/"), `[`, 7)
unique(samples$directory)
# [1] "2023-03-29_Transfer10x_SPage" "2023-02-09_Transfer10x_SPage"

dir.create(here::here("processed-data", "01_spaceranger", unique(samples$directory)[1]))
subfolders = samples[samples$directory == "2023-03-29_Transfer10x_SPage",]
for (i in seq_along(subfolders$samples)) {
  dir.create(here::here("processed-data", "01_spaceranger", "2023-03-29_Transfer10x_SPage", subfolders$samples[i]))
  file.symlink(subfolders$path[i], (here::here("processed-data", "01_spaceranger", "2023-03-29_Transfer10x_SPage", subfolders$samples[i])))
}

dir.create(here::here("processed-data", "01_spaceranger", unique(samples$directory)[2]))
subfolders = samples[samples$directory == "2023-02-09_Transfer10x_SPage",]
for (i in seq_along(subfolders$samples)) {
  dir.create(here::here("processed-data", "01_spaceranger", "2023-02-09_Transfer10x_SPage", subfolders$samples[i]))
  file.symlink(subfolders$path[i], (here::here("processed-data", "01_spaceranger", "2023-02-09_Transfer10x_SPage", subfolders$samples[i])))
}

