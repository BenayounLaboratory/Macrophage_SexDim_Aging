setwd('/Volumes/BB_HQ_1/Immune_sex_dimorphism_Aging/Macrophages/Clean_Macrophage_Data/RNAseq/DEseq2_analysis/DESeq2')
options(stringsAsFactors = FALSE)

# load libraries for analysis
library(DESeq2)
library(pheatmap)
library('bitops')
library('limma')
library(RColorBrewer)
library(fields)


# 2024-11-01
# plot DE heatmaps

my.cts <- read.table('2022-03-30_Peritoneal_Macrophages_DESeq2_model_with_AGE_SEX_log2_counts_matrix.txt')

M.age.spe <- read.table('2022-03-30_Peritoneal_Macrophages_AGING_Male.NOT.Female_FDR5_genes_statistics.txt')
F.age.spe <- read.table('2022-03-30_Peritoneal_Macrophages_AGING_Female.NOT.Male_FDR5_genes_statistics.txt')



pdf(paste(Sys.Date(),"Male_Specific_AGING_Heatmap_FDR5_RNAseq.pdf", sep = "_"), onefile = F, height = 6, width = 6)
my.heatmap.title <- paste("Peritoneal_Macrophages M specific aging genes", nrow(M.age.spe),sep="")
pheatmap(my.cts[rownames(M.age.spe),],
         cluster_cols = F,
         cluster_rows = T,
         colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
         show_rownames = F, scale="row",
         main = my.heatmap.title, 
         cellwidth = 10,
         cellheight = 0.2)
dev.off()




pdf(paste(Sys.Date(),"Female_Specific_AGING_Heatmap_FDR5_RNAseq.pdf", sep = "_"), onefile = F, height = 6, width = 6)
my.heatmap.title <- paste("Peritoneal_Macrophages F specific aging genes", nrow(F.age.spe),sep="")
pheatmap(my.cts[rownames(F.age.spe),],
         cluster_cols = F,
         cluster_rows = T,
         colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
         show_rownames = F, scale="row",
         main = my.heatmap.title, 
         cellwidth = 10,
         cellheight = 0.2)
dev.off()


#######################
sink(file = paste(Sys.Date(),"Macrophage_aging_sex_RNAseq_HEATMAPS_session_Info.txt", sep =""))
sessionInfo()
sink()

