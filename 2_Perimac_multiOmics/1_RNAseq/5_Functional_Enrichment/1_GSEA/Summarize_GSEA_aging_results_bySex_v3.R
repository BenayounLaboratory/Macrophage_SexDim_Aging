setwd('/Volumes/BB_HQ_1//Immune_sex_dimorphism_Aging/Macrophages/Clean_Macrophage_Data/RNAseq/Functional_enrichment/GSEA_v2/Summarize/')
options(stringsAsFactors = F)

source('function_sexAging_GSEA_summary_v5.R')

##############################################################################################################################
my.m.GO    <- read.csv("../2024-07-11_PeriMac_Aging_Males_EMSEMBL_GO_FDR_5_Phenotest_GSEA_Analysis_table_773_significant_OUTPUT_ALL.txt", header = T, sep = "\t")
my.f.GO    <- read.csv("../2024-07-11_PeriMac_Aging_Females_EMSEMBL_GO_FDR_5_Phenotest_GSEA_Analysis_table_868_significant_OUTPUT_ALL.txt"  , header = T, sep = "\t")
summarize_SexAge( my.m.GO, my.f.GO, "GO_ALL")

my.m.react <- read.csv("../2024-07-11_PeriMac_Aging_Males_MSigDB_Reactome_FDR_5_Phenotest_GSEA_Analysis_table_329_significant_OUTPUT_ALL.txt", header = T, sep = "\t")
my.f.react <- read.csv("../2024-07-11_PeriMac_Aging_Females_MSigDB_Reactome_FDR_5_Phenotest_GSEA_Analysis_table_174_significant_OUTPUT_ALL.txt"  , header = T, sep = "\t")
summarize_SexAge( my.m.react, my.f.react, "REACTOME")


#######################
sink(file = paste(Sys.Date(),"_GSEA_AgingSex_session_Info.txt", sep =""))
sessionInfo()
sink()


