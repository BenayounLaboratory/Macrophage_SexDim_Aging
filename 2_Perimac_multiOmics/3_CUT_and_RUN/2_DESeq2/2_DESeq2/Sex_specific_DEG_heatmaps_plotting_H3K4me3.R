setwd('/Volumes/BB_HQ_1/Immune_sex_dimorphism_Aging/Macrophages/Clean_Macrophage_Data/CUT_and_RUN/H3K4me3_Aging/H3K4me3_height/DEseq2_H3K4me3_height')
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

my.cts <- read.table('2021-09-30_PeriMac_H3K4me3_height_DESeq2_Analysis__log2_counts_matrix_DEseq2.txt')

# reorder
my.cts <- my.cts[,c("K4me3_PeriMac_YF_1",
                    "K4me3_PeriMac_YF_2",
                    "K4me3_PeriMac_YF_3",
                    "K4me3_PeriMac_YF_4",
                    "K4me3_PeriMac_OF_1",
                    "K4me3_PeriMac_OF_2",
                    "K4me3_PeriMac_OF_3",
                    "K4me3_PeriMac_YM_1",
                    "K4me3_PeriMac_YM_2",
                    "K4me3_PeriMac_YM_3",
                    "K4me3_PeriMac_YM_4",
                    "K4me3_PeriMac_YM_5",
                    "K4me3_PeriMac_OM_1",
                    "K4me3_PeriMac_OM_2",
                    "K4me3_PeriMac_OM_3",
                    "K4me3_PeriMac_OM_4",
                    "K4me3_PeriMac_OM_5")]

M.age.spe <- read.table('2021-09-30_PeriMac_H3K4me3_height_AGING_Male.NOT.Female_FDR5_genes_statistics.txt')
F.age.spe <- read.table('2021-09-30_PeriMac_H3K4me3_height_AGING_Female.NOT.Male_FDR5_genes_statistics.txt')




pdf(paste(Sys.Date(),"Male_Specific_AGING_Heatmap_FDR5_K4me3seq.pdf", sep = "_"), onefile = F, height = 6, width = 6) 
my.heatmap.title <- paste("Peritoneal_Macrophages M specific aging peaks ", nrow(M.age.spe),sep="")
pheatmap(my.cts[rownames(M.age.spe),],
         cluster_cols = F,
         cluster_rows = T,
         colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
         show_rownames = F, scale="row",
         main = my.heatmap.title, 
         cellwidth = 10,
         cellheight = 0.2)
dev.off()




pdf(paste(Sys.Date(),"Female_Specific_AGING_Heatmap_FDR5_K4me3seq.pdf", sep = "_"), onefile = F, height = 12, width = 6)
my.heatmap.title <- paste("Peritoneal_Macrophages F specific aging peaks ", nrow(F.age.spe),sep="")
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
sink(file = paste(Sys.Date(),"Macrophage_aging_sex_K4me3seq_HEATMAPS_session_Info.txt", sep =""))
sessionInfo()
sink()

