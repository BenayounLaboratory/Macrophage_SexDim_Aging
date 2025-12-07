setwd('/Users/berenice/Desktop/Irf2_CUTRUN/Comparison_RNAseq')
options(stringsAsFactors = F)

# load libraries for analysis
library(bitops)
library(ggvenn)
library(ggplot2)

# 2025-12-01
# Compare Irf2 CUT and RUN TF and Female aging

############################################################################################################################################
# Load TF targets peaks/genes
Irf2_peaks   <- read.csv('/Users/berenice/Desktop/Irf2_CUTRUN/MSPC/HOMER_MSPC_Irf2_PeriMac_4samples_ConsensusPeaks.xls' , header = T, sep = "\t")

# gene is called target if within 2kb of TSS
Irf2_targets_5kb     <- unique(Irf2_peaks$Gene.Name  [abs(Irf2_peaks  $Distance.to.TSS) < 5000])
Irf2_targets_5kb     <- Irf2_targets_5kb[!is.na(Irf2_targets_5kb)]

# Load female all aging
fem.age.all      <- read.table('/Volumes/BB_Travel_2/Mph_Aging/RNAseq/DEseq2_analysis/DESeq2/2022-03-30_Peritoneal_Macrophages_AGING_Female_ALL_genes_statistics.txt')
fem.all.up       <- rownames(fem.age.all)[bitAnd(fem.age.all$log2FoldChange_F > 0, fem.age.all$padj_F < 0.05) >0]
fem.all.dwn      <- rownames(fem.age.all)[bitAnd(fem.age.all$log2FoldChange_F < 0, fem.age.all$padj_F < 0.05) >0]

#### load expression data
tissue.cts <- read.table('/Volumes/BB_Travel_2/Mph_Aging/RNAseq/DEseq2_analysis/DESeq2//2022-03-30_Peritoneal_Macrophages_DESeq2_model_with_AGE_SEX_log2_counts_matrix.txt', sep = "\t", header = T)
###############################################################################################################################################


####################################################################################################
##### Overlap analysis with female aging

intersect(fem.all.up,Irf2_targets_5kb)
# [1] "Xaf1"    "Slfn5"   "Dhx58"   "Rsad2"   "Phf11d"  "Mx2"     "H2-D1"   "H2-T24"  "Enpp5"  
# [10] "Ifit3b"  "Ifih1"   "B2m"     "Oasl2"   "Usp18"   "Dnaaf3"  "Trim34a" "Trim30a" "Trim30d"
# [19] "Ddx60" 

intersect(fem.all.dwn,Irf2_targets_5kb)
#  "Jdp2"   "Mef2c"  "Vcan"   "Tm6sf1" "Irf2"  


########################
my.ups <- list("Irf2 targets (5 kb)"      = Irf2_targets_5kb  ,
               "Aging F up (FDR5)"        = fem.all.up )

overlap.up.test <- fisher.test(matrix(c(19,223,458,nrow(fem.age.all)-19-223-458),2,2), alternative = "greater")
overlap.up.test$p.value ### 0.0003307205

pdf(paste0(Sys.Date(),"_Irf2_CUTRUN_direct_targets_5kb_and_aging_F_up.pdf"))
ggvenn(my.ups, 
       fill_color = c("peru", "deeppink4"),
       stroke_size = 0.1, set_name_size = 4,
       auto_scale = T, show_stats = "c") + ggtitle(paste0("p = ",scientific(signif(overlap.up.test$p.value,3))))
dev.off()

########################
my.downs <- list("Irf2 targets (5 kb)"      = Irf2_targets_5kb  ,
                 "Aging F down (FDR5)"       = fem.all.dwn)


overlap.dwn.test <- fisher.test(matrix(c(5,472,110, nrow(fem.age.all)-5-472-110),2,2), alternative = "greater")
overlap.dwn.test$p.value ### 0.3136529

pdf(paste0(Sys.Date(),"_Irf2_CUTRUN_direct_targets_5kb_and_aging_F_down.pdf"))
ggvenn(my.downs, 
       fill_color = c("peru", "deeppink"),
       stroke_size = 0.1, set_name_size = 4,
       auto_scale = T, show_stats = "c") + ggtitle(paste0("p = ",signif(overlap.dwn.test$p.value,3)))
dev.off()



#### plot heatmap
fem.irf2 <- intersect(c(fem.all.up,fem.all.dwn),Irf2_targets_5kb)

pdf(paste0(Sys.Date(),"_Female_Aging_FDR5_and_Irf2_CUTRUN_targets.pdf"), height = 6, width = 6)
pheatmap::pheatmap(tissue.cts[fem.irf2,], scale = "row", cluster_cols = F, 
                   colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
                   main = "Female age-regulation, Irf2 CUTRUN targets", cellheight = 12, cellwidth = 12)
dev.off()

###############################################################################################################################################

#######################
sink(file = paste(Sys.Date(),"_F_aging_Irf2_CUTandRUN_Venn_session_Info.txt", sep =""))
sessionInfo()
sink()


