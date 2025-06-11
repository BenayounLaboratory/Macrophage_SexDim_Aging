setwd('/Volumes/BB_HQ_1/Immune_sex_dimorphism_Aging/Macrophages/Clean_Macrophage_Data/RNAseq/Functional_enrichment/Top_TF_Validation/Irf2_ChIPseq_BMDM//')
options(stringsAsFactors = F)

# load libraries for analysis
library(bitops)
library(ggvenn)
library(ggplot2)

# 2024-08-19
# Compare Irf2 BMDM TF ChIP and RAW shRNA
# Use consistent peak seats from MSPC

# 2025-01-30
# compare with Female aging

############################################################################################################################################
# Load TF targets peaks/genes
Irf2_peaks   <- read.csv('/Volumes/BB_HQ_1/Immune_sex_dimorphism_Aging/Macrophages/Clean_Macrophage_Data/ChIP-seq/MSPC/HOMER_MSPC_Irf2_BMDM_ConsensusPeaks.xls' , header = T, sep = "\t")

# gene is called target if within 2kb of TSS
Irf2_targets_500bp   <- unique(Irf2_peaks$Gene.Name  [abs(Irf2_peaks  $Distance.to.TSS) < 500 ]) # 82
Irf2_targets_1kb     <- unique(Irf2_peaks$Gene.Name  [abs(Irf2_peaks  $Distance.to.TSS) < 1000]) # 104
Irf2_targets_2kb     <- unique(Irf2_peaks$Gene.Name  [abs(Irf2_peaks  $Distance.to.TSS) < 2000]) # 116
Irf2_targets_5kb     <- unique(Irf2_peaks$Gene.Name  [abs(Irf2_peaks  $Distance.to.TSS) < 5000]) # 146

# Load female all aging
fem.age.all <- read.table('/Volumes/BB_HQ_1/Immune_sex_dimorphism_Aging/Macrophages/Clean_Macrophage_Data/RNAseq/DEseq2_analysis/DESeq2/2022-03-30_Peritoneal_Macrophages_AGING_Female_ALL_genes_statistics.txt')
fem.all.up       <- rownames(fem.age.all)[bitAnd(fem.age.all$log2FoldChange_F > 0, fem.age.all$padj_F < 0.05) >0]
fem.all.dwn      <- rownames(fem.age.all)[bitAnd(fem.age.all$log2FoldChange_F < 0, fem.age.all$padj_F < 0.05) >0]

fem.all.up.fdr10      <- rownames(fem.age.all)[bitAnd(fem.age.all$log2FoldChange_F > 0, fem.age.all$padj_F < 0.1) >0]
fem.all.dwn.fdr10     <- rownames(fem.age.all)[bitAnd(fem.age.all$log2FoldChange_F < 0, fem.age.all$padj_F < 0.1) >0]


#### load expression data
tissue.cts <- read.table('/Volumes/BB_HQ_1/Immune_sex_dimorphism_Aging/Macrophages/Clean_Macrophage_Data/RNAseq/DEseq2_analysis/DEseq2/2022-03-30_Peritoneal_Macrophages_DESeq2_model_with_AGE_SEX_log2_counts_matrix.txt', sep = "\t", header = T)
###############################################################################################################################################


####################################################################################################
##### Overlap analysis with female aging

intersect(fem.all.up,Irf2_targets_5kb)
# [1] "Xaf1"    "Dhx58"   "Mx2"     "H2-T24"  "Ifit3"   "B2m"     "Oasl2"   "Oas3"    "Usp18"   "Trim34a" "Ddx60"  

setdiff(intersect(fem.all.up.fdr10,Irf2_targets_5kb),intersect(fem.all.up,Irf2_targets_5kb))
# "Stat2"  "Cmpk2"  "Rtp4"   "H2-T10" "Adar"  

intersect(fem.all.dwn,Irf2_targets_5kb)
# "Tm6sf1"

setdiff(intersect(fem.all.dwn.fdr10,Irf2_targets_5kb),intersect(fem.all.dwn,Irf2_targets_5kb))
# 0


########################
my.ups <- list("Irf2 targets (5 kb)"      = Irf2_targets_5kb  ,
               "Aging F up (FDR5)"        = fem.all.up )

overlap.up.test <- fisher.test(matrix(c(11,135,230,nrow(fem.age.all)-11-135-230),2,2), alternative = "greater")
overlap.up.test$p.value ### 2.721893e-05

pdf(paste0(Sys.Date(),"_Irf2_direct_targets_5kb_and_aging_F_up.pdf"))
ggvenn(my.ups, 
       fill_color = c("peru", "deeppink4"),
       stroke_size = 0.5, set_name_size = 4,
       auto_scale = T, show_stats = "c") + ggtitle(paste0("p = ",signif(overlap.up.test$p.value,3)))
dev.off()

########################
my.downs <- list("Irf2 targets (5 kb)"      = Irf2_targets_5kb  ,
                 "Aging F down (FDR5)"       = fem.all.dwn)


overlap.dwn.test <- fisher.test(matrix(c(1,145,113, nrow(fem.age.all)-1-145-113),2,2), alternative = "greater")
overlap.dwn.test$p.value ### 0.6788519

pdf(paste0(Sys.Date(),"_Irf2_direct_targets_5kb_and_aging_F_down.pdf"))
ggvenn(my.downs, 
       fill_color = c("peru", "deeppink"),
       stroke_size = 0.5, set_name_size = 4,
       auto_scale = T, show_stats = "c") + ggtitle(paste0("p = ",signif(overlap.dwn.test$p.value,3)))
dev.off()



#### plot heatmap
fem.irf2 <- intersect(c(fem.all.up,fem.all.dwn),Irf2_targets_5kb)

pdf(paste0(Sys.Date(),"_Female_Aging_FDR5_and_Irf2_ChIPseq_targets.pdf"))
pheatmap::pheatmap(tissue.cts[fem.irf2,], scale = "row", cluster_cols = F, 
                   colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
                   main = "Female age-regulation, Irf2 ChIP targets", cellheight = 12, cellwidth = 12)
dev.off()


###############################################################################################################################################

#######################
sink(file = paste(Sys.Date(),"_F_aging_Irf2_ChIP_session_Info.txt", sep =""))
sessionInfo()
sink()


