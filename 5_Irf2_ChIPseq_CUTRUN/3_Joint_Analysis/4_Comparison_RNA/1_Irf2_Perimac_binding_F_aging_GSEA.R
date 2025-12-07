setwd('/Users/berenice/Desktop/Irf2_CUTRUN/Comparison_RNAseq/Comparison_F_Aging')
options(stringsAsFactors = F)

# load libraries for analysis
library(DESeq2)
library(Vennerable)
library(bitops)
library(phenoTest)

# 2025-12-01
# Compare Irf2 PeriMac CUTandRUN and female aging
# Use consistent peak seats from MSPC


# Load TF targets peaks/genes
Irf2_peaks   <- read.csv('/Users/berenice/Desktop/Irf2_CUTRUN/MSPC/HOMER_MSPC_Irf2_PeriMac_4samples_ConsensusPeaks.xls' , header = T, sep = "\t")

# gene is called target if within 2kb of TSS
Irf2_targets_500bp   <- unique(Irf2_peaks$Gene.Name  [abs(Irf2_peaks  $Distance.to.TSS) < 500 ]) # 239
Irf2_targets_1kb     <- unique(Irf2_peaks$Gene.Name  [abs(Irf2_peaks  $Distance.to.TSS) < 1000]) # 299
Irf2_targets_2kb     <- unique(Irf2_peaks$Gene.Name  [abs(Irf2_peaks  $Distance.to.TSS) < 2000]) # 366
Irf2_targets_5kb     <- unique(Irf2_peaks$Gene.Name  [abs(Irf2_peaks  $Distance.to.TSS) < 5000]) # 478

# Load female all aging
fem.age.all      <- read.table('/Volumes/BB_Travel_2/Mph_Aging/RNAseq/DEseq2_analysis/DESeq2/2022-03-30_Peritoneal_Macrophages_AGING_Female_ALL_genes_statistics.txt')
fem.all.up       <- rownames(fem.age.all)[bitAnd(fem.age.all$log2FoldChange_F > 0, fem.age.all$padj_F < 0.05) >0]
fem.all.dwn      <- rownames(fem.age.all)[bitAnd(fem.age.all$log2FoldChange_F < 0, fem.age.all$padj_F < 0.05) >0]
###############################################################################################################################################


###############################################################################################################################################
# Geneset enrichment
res.fem           <- data.frame(fem.age.all)
res.fem$GeneName  <- rownames(fem.age.all )

######################## B. Prepare GeneLists using DEseq2 t-statistic to rank genes ########################
res.fem.geneList         <- res.fem$stat
names(res.fem.geneList)  <- res.fem$GeneName
res.fem.geneList         <- sort(res.fem.geneList , decreasing = TRUE)


######################## C. Prep gene Sets ########################
# load gene lists

# Load TF target genes
my.TF.gs <- list( "Irf2_targets_500bp"   = setdiff(Irf2_targets_500bp  ,"") ,
                  "Irf2_targets_1kb"     = setdiff(Irf2_targets_1kb    ,"") ,
                  "Irf2_targets_2kb"     = setdiff(Irf2_targets_2kb    ,"") ,
                  "Irf2_targets_5kb"     = setdiff(Irf2_targets_5kb    ,"") )


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
# Irf2_targets_500bp 239 0.5191921 2.341849       0        0   0
# Irf2_targets_1kb   299 0.4812641 2.205479       0        0   0
# Irf2_targets_2kb   366 0.4573264 2.137614       0        0   0
# Irf2_targets_5kb   478 0.4289279 2.056152       0        0   0

pdf(paste(Sys.Date(), "Macrophage_Irf2_CUTandRUN_DirectTarget5kb_List_GSEA_plot_FEMALE_aging.pdf", sep = "_"))
plot.gseaData(gsea.res.sex, es.nes='nes', selGsets='Irf2_targets_5kb', color = "purple")
dev.off()



############################################################################################
# Make bubble chart summary
library(ggplot2) 
library(scales) 
theme_set(theme_bw())

############## Plot TFs
gsea.summary    <- data.frame(gsea.res.sex$significance$summary)
gsea.summary.v2 <- cbind(rownames(gsea.summary),gsea.summary)
colnames(gsea.summary.v2)[1] <- "GeneSet"

# get merged datafame for ggplot
gsea.summary.v2$minusLog10FDR <- -log10(gsea.summary.v2$fdr + 1e-30) ### to avoid -Inf
gsea.summary.v2$condition <- rep("F_aging",nrow(gsea.summary.v2))

my.max <- 3
my.min <- -3
my.values <- c(my.min,0.75*my.min,0.5*my.min,0.25*my.min,0,0.25*my.max,0.5*my.max,0.75*my.max,my.max)
my.scaled <- rescale(my.values, to = c(0, 1))
my.color.vector <- c("darkblue","dodgerblue4","dodgerblue3","dodgerblue1","white","lightcoral","brown1","firebrick2","firebrick4")

# to preserve the wanted order
gsea.summary.v2$condition <- factor(gsea.summary.v2$condition, levels = unique(gsea.summary.v2$condition))
gsea.summary.v2$GeneSet   <- factor(gsea.summary.v2$GeneSet, levels = rev(unique(gsea.summary.v2$GeneSet)))

pdf(paste0(Sys.Date(),"F_Aging_GSEA_Irf2_CUTandRUN_MSPCs_GeneSets.pdf"),height = 4, width=5)
my.plot <- ggplot(gsea.summary.v2,aes(x=condition,y=GeneSet,colour=nes,size=minusLog10FDR))+ theme_bw()+ geom_point(shape = 16)
my.plot <- ggplot(gsea.summary.v2,aes(x=condition,y=GeneSet,colour=nes,size=minusLog10FDR))+ theme(text = element_text(size=16))+ geom_point(shape = 16)
my.plot <- my.plot + ggtitle("GSEA") + labs(x = "F aging", y = "Irf2 ChIP target GeneSet")
my.plot <- my.plot + scale_colour_gradientn(colours = my.color.vector,space = "Lab", na.value = "grey50", guide = "colourbar", values = my.scaled, limits = c(my.min,my.max))
my.plot <- my.plot + scale_size_area(limits = c(1,30))
print(my.plot)
dev.off()


#######################
sink(file = paste(Sys.Date(),"_F_aging_Irf2_CUTandRUN_session_Info.txt", sep =""))
sessionInfo()
sink()


