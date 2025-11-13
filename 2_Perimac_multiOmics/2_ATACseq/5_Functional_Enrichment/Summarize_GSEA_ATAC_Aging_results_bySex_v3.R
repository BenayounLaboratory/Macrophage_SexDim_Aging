setwd('/Volumes/BB_HQ_1/Immune_sex_dimorphism_Aging/Macrophages/Clean_Macrophage_Data/ATACseq/Functional_Enrichment_Analysis/GSEA_v2//Summarize/')
options(stringsAsFactors = F)

# 2023-11-01
# Analysis to compare male and female aging (H3K4me3 height)
# Summarize common and specifc pathways

# 2024-07-11
# Update plotting parameters (less pathways shown, preset scale to harmonize across omics for figures, etc.)

# 2024-12-06
# replot

source('function_sexAging_GSEA_summary_ATAC_v4.R')

##############################################################################################################################
my.m.GO    <- read.csv("../2024-12-06_ATAC_PeriMac_Aging_Males_EMSEMBL_GO_FDR_5_Phenotest_GSEA_Analysis_table_1151_significant_OUTPUT_ALL.txt", header = T, sep = "\t")
my.f.GO    <- read.csv("../2024-12-06_ATAC_PeriMac_Aging_Females_EMSEMBL_GO_FDR_5_Phenotest_GSEA_Analysis_table_755_significant_OUTPUT_ALL.txt"  , header = T, sep = "\t")
summarize_SexAge( my.m.GO, my.f.GO, "GO_ALL")

my.m.react <- read.csv("../2024-12-06_ATAC_PeriMac_Aging_Males_MSigDB_Reactome_FDR_5_Phenotest_GSEA_Analysis_table_286_significant_OUTPUT_ALL.txt", header = T, sep = "\t")
my.f.react <- read.csv("../2024-12-06_ATAC_PeriMac_Aging_Females_MSigDB_Reactome_FDR_5_Phenotest_GSEA_Analysis_table_84_significant_OUTPUT_ALL.txt"  , header = T, sep = "\t")
summarize_SexAge( my.m.react, my.f.react, "REACTOME")

my.m.hall    <- read.csv("../2024-12-06_ATAC_PeriMac_Aging_Males_MSigDB_Hallmarks_FDR_5_Phenotest_GSEA_Analysis_table_9_significant_OUTPUT_ALL.txt", header = T, sep = "\t")
my.f.hall    <- read.csv("../2024-12-06_ATAC_PeriMac_Aging_Females_MSigDB_Hallmarks_FDR_5_Phenotest_GSEA_Analysis_table_2_significant_OUTPUT_ALL.txt"  , header = T, sep = "\t")
summarize_SexAge( my.m.hall, my.f.hall, "MSigDB_hallmarks")

#######################
sink(file = paste(Sys.Date(),"_GSEA_ATAC_AgingSex_session_Info_STRIGENT.txt", sep =""))
sessionInfo()
sink()


