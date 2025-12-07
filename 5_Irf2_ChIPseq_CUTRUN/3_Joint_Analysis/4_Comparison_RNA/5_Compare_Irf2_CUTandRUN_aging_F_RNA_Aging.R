setwd('/Users/berenice/Desktop/Irf2_CUTRUN/Comparison_RNAseq')
options(stringsAsFactors = F)

# load libraries for analysis
library(bitops)
library(ggvenn)
library(ggplot2)
library(scales)

# load libraries for analysis
library(DESeq2)
library(phenoTest)

# 2025-12-02
# Compare Irf2 CUT and RUN aging vs Female RNA aging

############################################################################################################################################
# Load TF targets peaks/genes
Irf2_peaks.fdr5   <- read.csv('/Users/berenice/Desktop/Irf2_CUTRUN/DiffBind/MSPC_4samples/2025-11-28_PeriMac_Irf2_stringent_Aging_DiffBind_DESeq2_FDR5_Peaks.txt' , header = T, sep = "\t")

# gene is called target if within 2kb of TSS
Irf2_difftargets_5kb     <- unique(Irf2_peaks.fdr5$Gene.Name  [abs(Irf2_peaks.fdr5  $Distance.to.TSS) < 5000])
Irf2_difftargets_5kb     <- Irf2_difftargets_5kb[!is.na(Irf2_difftargets_5kb)] # 167

# Load female all aging
fem.age.all      <- read.table('/Volumes/BB_Travel_2/Mph_Aging/RNAseq/DEseq2_analysis/DESeq2/2022-03-30_Peritoneal_Macrophages_AGING_Female_ALL_genes_statistics.txt')
fem.all.up       <- rownames(fem.age.all)[bitAnd(fem.age.all$log2FoldChange_F > 0, fem.age.all$padj_F < 0.05) >0]
fem.all.dwn      <- rownames(fem.age.all)[bitAnd(fem.age.all$log2FoldChange_F < 0, fem.age.all$padj_F < 0.05) >0]

#### load expression data
tissue.cts <- read.table('/Volumes/BB_Travel_2/Mph_Aging/RNAseq/DEseq2_analysis/DESeq2//2022-03-30_Peritoneal_Macrophages_DESeq2_model_with_AGE_SEX_log2_counts_matrix.txt', sep = "\t", header = T)


###############################################################################################################################################


####################################################################################################
##### Overlap analysis with female aging

intersect(fem.all.up,Irf2_difftargets_5kb)
# [1] "Dhx58"   "Mx2"     "H2-D1"   "B2m"     "Dnaaf3"  "Trim30a"

########################
my.ups <- list("Irf2 targets (FDR5)"      = Irf2_difftargets_5kb  ,
               "Aging F up (FDR5)"        = fem.all.up )

overlap.up.test <- fisher.test(matrix(c(5,160,236,nrow(fem.age.all)-5-160-236),2,2), alternative = "greater")
overlap.up.test$p.value ### 0.1328872

pdf(paste0(Sys.Date(),"_Irf2_CUTRUN_direct_targets_FDR5aging_and_aging_F_up.pdf"))
ggvenn(my.ups, 
       fill_color = c("peru", "deeppink4"),
       stroke_size = 0.1, set_name_size = 4,
       auto_scale = T, show_stats = "c") + ggtitle(paste0("p = ",scientific(signif(overlap.up.test$p.value,3))))
dev.off()
###############################################################################################################################################

###############################################################################################################################################
#### B. GSEA analysis

# Geneset enrichment
res.fem           <- data.frame(fem.age.all)
res.fem$GeneName  <- rownames(fem.age.all )

######################## B. Prepare GeneLists using DEseq2 t-statistic to rank genes ########################
res.fem.geneList         <- res.fem$stat
names(res.fem.geneList)  <- res.fem$GeneName
res.fem.geneList         <- sort(res.fem.geneList , decreasing = TRUE)


######################## C. Prep gene Sets ########################
# Load Irf target genes with binding change
my.TF.gs <- list( "Irf2_targets_FDR5"    = setdiff(Irf2_difftargets_5kb  ,"")  )


######################## C. Gene Set Enrichment Analysis ########################
# set seed to stabilize output
set.seed(1234567890)

# run phenotest GSEA
gsea.res.sex <- gsea( x         =  res.fem.geneList   ,
                      gsets     =  my.TF.gs           ,
                      mc.cores  =  1                  ,
                      logScale  =  FALSE              ,
                      B         =  10000              ,
                      minGenes  =  5                  ,
                      maxGenes  =  15000               )
my.summary <- data.frame(summary(gsea.res.sex))
gsea.res.sex$significance$summary
#                      n        es      nes pval.es pval.nes fdr
# Irf2_targets_FDR5 166 0.4191537 1.778311 8.358639e-08 8.288565e-08 8.288565e-08

pdf(paste(Sys.Date(), "Macrophage_Irf2_CUTandRUN_DirectTarget5kb_List_GSEA_plot_FEMALE_aging.pdf", sep = "_"))
plot.gseaData(gsea.res.sex, es.nes='nes', selGsets='Irf2_targets_FDR5', color = "purple")
dev.off()

###############################################################################################################################################

#######################
sink(file = paste(Sys.Date(),"_F_aging_Irf2_CUTandRUN_Venn_session_Info.txt", sep =""))
sessionInfo()
sink()


