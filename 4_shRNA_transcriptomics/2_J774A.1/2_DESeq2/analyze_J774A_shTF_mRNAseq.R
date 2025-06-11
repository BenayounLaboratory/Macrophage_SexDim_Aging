setwd('/Volumes/BB_HQ_1/Immune_sex_dimorphism_Aging/Macrophages/Clean_Macrophage_Data/Intervention_Data/shRNA_analyses/J774_cells/DESeq2')
options(stringsAsFactors = F)

library(DESeq2)
library(pheatmap)
library(beeswarm)
library(pvclust)
library(bitops)

source('shTF_process_function_J774.R')

# 2025-01-30
# Analyze J774 sh IRF2 RNAseq

##################################################################
# 1. Read STAR counts and format data for DESeq2 processing
data.irf2     <- read.table('../STAR/2025-01-23_J774A1_shIRF2_CLEAN_counts_unstranded.txt', sep = "\t", header = T)

# 2. process RNAseq
process_shTF(data.irf2     , "Irf2"   , shNon = 1:3, shLuc = 4:6, shTF1 = 7:9, shTF2 = 10:12)



#######################
sink(file = paste(Sys.Date(),"shTF_J774cells_RNAseq_session_Info.txt", sep =""))
sessionInfo()
sink()
