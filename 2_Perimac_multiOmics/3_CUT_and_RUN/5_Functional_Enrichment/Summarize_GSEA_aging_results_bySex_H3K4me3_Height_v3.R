setwd('/Volumes/BB_HQ_1/Immune_sex_dimorphism_Aging/Macrophages/Clean_Macrophage_Data/CUT_and_RUN/H3K4me3_Aging/H3K4me3_height/Functional_Enrichment/GSEA_v2/Summarize')
options(stringsAsFactors = F)

# 2023-11-01
# Analysis to compare male and female aging (H3K4me3 height)
# Summarize common and specifc pathways

# 2024-07-11
# Update plotting parameters (less pathways shown, preset scale to harmonize across omics for figures, etc.)

library(ggplot2)
library(scales) 

source('function_sexAging_GSEA_summary_H3K4me3_Height_v2.R')

##############################################################################################################################
my.m.react <- read.csv("../2024-07-11_H3K4me3_height_PeriMac_Aging_Males_MSigDB_Reactome_FDR_5_Phenotest_GSEA_Analysis_table_187_significant_OUTPUT_ALL.txt", header = T, sep = "\t")
my.f.react <- read.csv("../2024-07-11_H3K4me3_height_PeriMac_Aging_Females_MSigDB_Reactome_FDR_5_Phenotest_GSEA_Analysis_table_423_significant_OUTPUT_ALL.txt"  , header = T, sep = "\t")
summarize_SexAge_peak( my.m.react, my.f.react, "REACTOME")

my.m.GO    <- read.csv("../2024-07-11_H3K4me3_height_PeriMac_Aging_Males_EMSEMBL_GO_FDR_5_Phenotest_GSEA_Analysis_table_941_significant_OUTPUT_ALL.txt", header = T, sep = "\t")
my.f.GO    <- read.csv("../2024-07-11_H3K4me3_height_PeriMac_Aging_Females_EMSEMBL_GO_FDR_5_Phenotest_GSEA_Analysis_table_1349_significant_OUTPUT_ALL.txt"  , header = T, sep = "\t")
summarize_SexAge_peak( my.m.GO, my.f.GO, "GO_ALL")

my.m.hall    <- read.csv("../2024-07-11_H3K4me3_height_PeriMac_Aging_Males_MSigDB_Hallmarks_FDR_5_Phenotest_GSEA_Analysis_table_14_significant_OUTPUT_ALL.txt", header = T, sep = "\t")
my.f.hall    <- read.csv("../2024-07-11_H3K4me3_height_PeriMac_Aging_Females_MSigDB_Hallmarks_FDR_5_Phenotest_GSEA_Analysis_table_21_significant_OUTPUT_ALL.txt"  , header = T, sep = "\t")
summarize_SexAge_peak( my.m.hall, my.f.hall, "MsigDB_Hallmarks")

#######################
sink(file = paste(Sys.Date(),"_GSEA_AgingSex_H3K4me3_height_session_Info_STRINGENT.txt", sep =""))
sessionInfo()
sink()
