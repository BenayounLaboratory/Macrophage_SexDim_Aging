setwd('/Volumes/BB_HQ_1/Immune_sex_dimorphism_Aging/Macrophages/Clean_Macrophage_Data/RNAseq/Functional_enrichment/Top_TF_Validation/Irf2_ChIPseq_BMDM/')
options(stringsAsFactors = F)

# load libraries for analysis
library(DESeq2)
library(Vennerable)
library(bitops)
library(phenoTest)

# 2024-06-18
# Compare Irf2 BMDM TF ChIP and female aging
# Use consistent peak seats from MSPC

# 2024-08-19
# get GSEA plot

# 2025-01-30
# rerun cehck

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
#                     n        es      nes pval.es     pval.nes          fdr
# Irf2_targets_500bp  82 0.6722671 2.643124       0 0.000000e+00 0.000000e+00
# Irf2_targets_1kb   104 0.6616832 2.638718       0 2.220446e-16 2.960595e-16
# Irf2_targets_2kb   116 0.6506689 2.641868       0 2.220446e-16 4.440892e-16
# Irf2_targets_5kb   146 0.5935683 2.482883       0 0.000000e+00 0.000000e+00

pdf(paste(Sys.Date(), "Macrophage_Irf2_DirectTarget5kb_List_GSEA_plot_FEMALE_aging.pdf", sep = "_"))
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

pdf(paste0(Sys.Date(),"F_Aging_GSEA_Irf2_MSPCs_ChIP_GeneSets.pdf"),height = 4, width=5)
my.plot <- ggplot(gsea.summary.v2,aes(x=condition,y=GeneSet,colour=nes,size=minusLog10FDR))+ theme_bw()+ geom_point(shape = 16)
my.plot <- ggplot(gsea.summary.v2,aes(x=condition,y=GeneSet,colour=nes,size=minusLog10FDR))+ theme(text = element_text(size=16))+ geom_point(shape = 16)
my.plot <- my.plot + ggtitle("GSEA") + labs(x = "F aging", y = "Irf2 ChIP target GeneSet")
my.plot <- my.plot + scale_colour_gradientn(colours = my.color.vector,space = "Lab", na.value = "grey50", guide = "colourbar", values = my.scaled, limits = c(my.min,my.max))
my.plot <- my.plot + scale_size_area(limits = c(1,30))
print(my.plot)
dev.off()


#######################
sink(file = paste(Sys.Date(),"_F_aging_Irf2_ChIP_session_Info.txt", sep =""))
sessionInfo()
sink()


