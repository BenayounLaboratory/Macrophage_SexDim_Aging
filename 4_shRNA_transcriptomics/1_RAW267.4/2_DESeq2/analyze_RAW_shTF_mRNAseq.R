setwd('/Volumes/BB_HQ_1/Immune_sex_dimorphism_Aging/Macrophages/Clean_Macrophage_Data/Intervention_Data/shRNA_analyses/DEseq2/')
options(stringsAsFactors = F)

library(DESeq2)
library(pheatmap)
library(beeswarm)
library(pvclust)
library(bitops)

source('shTF_process_function.R')

# 2024-05-23
# Analyze RAW shRNA data
# Meis1 and Mef2c KD
# non targetting and luc shRNA as controls
# analyze the overlap between 2 hairpins

# 2024-08-22
# correct MEf2c column order

##################################################################
# 1. Read STAR counts and format data for DESeq2 processing
data.irf2     <- read.table('../STAR/2024-08-15_RAW_shIRF2_CLEAN_counts_unstranded.txt', sep = "\t", header = T)
data.tbl1xr1  <- read.table('../STAR/2024-08-15_RAW_shTBL1XR1_CLEAN_counts_unstranded.txt', sep = "\t", header = T)
data.tal1     <- read.table('../STAR/2024-08-15_RAW_shTAL1_CLEAN_counts_unstranded.txt', sep = "\t", header = T)
data.meis1    <- read.table('../STAR/22024-08-15_RAW_shMeis1_CLEAN_counts_unstranded.txt', sep = "\t", header = T)
data.mef2c    <- read.table('../STAR/22024-08-15_RAW_shMef2c_CLEAN_counts_unstranded.txt', sep = "\t", header = T)

# reorder controls
data.mef2c <- data.mef2c[,c(1:6,16:18,7:15)]

# 2. process RNAseq
process_shTF(data.irf2     , "Irf2"   , shNon = 1:3, shLuc = 4:6, shTF1 = 7:9, shTF2 = 10:12)
process_shTF(data.tbl1xr1  , "Tbl1xr1", shNon = 1:3, shLuc = 4:6, shTF1 = 7:9, shTF2 = 10:11)
process_shTF(data.tal1     , "Tal1"   , shNon = 1:3, shLuc = 4:6, shTF1 = 7:9, shTF2 = 10:12)
process_shTF(data.meis1    , "Meis1"  , shNon = 1:3, shLuc = 4:6, shTF1 = 7:9, shTF2 = 10:12)
process_shTF(data.mef2c    , "Mef2c"  , shNon = 1:3, shLuc = 4:6, shTF1 = 7:9, shTF2 = 10:12)



#######################
sink(file = paste(Sys.Date(),"shTF_RAWcells_RNAseq_session_Info.txt", sep =""))
sessionInfo()
sink()




