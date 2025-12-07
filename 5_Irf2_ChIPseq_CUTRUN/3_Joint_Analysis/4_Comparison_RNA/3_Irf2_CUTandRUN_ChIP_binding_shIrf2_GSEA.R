setwd('/Users/berenice/Desktop/Irf2_CUTRUN/Comparison_RNAseq')
options(stringsAsFactors = F)

library(DESeq2)
library(Vennerable)
library(bitops)
library(phenoTest)

# 2025-12-01
# Compare Irf2 BMDM ChIP/PeriMac CUTandRUN and shIrf2 RAW/J774
# Use consistent peak seats from MSPC


###############################################################################################################################################
# Load TF targets peaks/genes
# Load TF targets peaks/genes
Irf2_perimac   <- read.csv('/Users/berenice/Desktop/Irf2_CUTRUN/MSPC/HOMER_MSPC_Irf2_PeriMac_4samples_ConsensusPeaks.xls' , header = T, sep = "\t")
Irf2_bmdms     <- read.csv('/Users/berenice/Desktop/Irf2_CUTRUN/01_ChIPseq/MSPC/HOMERMSPC_Irf2_BMDM_ConsensusPeaks.xls' , header = T, sep = "\t")

# gene is called target if within 5kb of TSS
Irf2_5kb_perimac     <- unique(Irf2_perimac$Gene.Name  [abs(Irf2_perimac  $Distance.to.TSS) < 5000]) # 478
Irf2_5kb_bmdms       <- unique(Irf2_bmdms$Gene.Name    [abs(Irf2_bmdms    $Distance.to.TSS) < 5000]) # 189

# Load shRNA RNAseq
shIrf2.raw  <- read.table('/Volumes/BB_Travel_2/Mph_Aging/Intervention_Data/shRNA_analyses/RAW267_cells/DEseq2/2024-08-16_RAW_sh_Irf2_data_DESeq2_shIrf2_DEseq2_results_COMBINED_all_genes_statistics.txt', header = T)
shIrf2.j774 <- read.table('/Volumes/BB_Travel_2/Mph_Aging/Intervention_Data/shRNA_analyses/J774_cells/DESeq2/2025-01-30_J774_sh_Irf2_data_DESeq2_shIrf2_DEseq2_results_COMBINED_all_genes_statistics.txt' , header = T)

# get mean stat from both shRNAs
shIrf2.raw$av_stat  <- apply(shIrf2.raw[,c("stat.sh1", "stat.sh2")],1,mean)
shIrf2.j774$av_stat <- apply(shIrf2.j774[,c("stat.sh1", "stat.sh2")],1,mean)

######################## A. Prepare GeneLists using DEseq2 t-statistic to rank genes ########################
raw.irf2.geneList         <- shIrf2.raw$av_stat
names(raw.irf2.geneList)  <- shIrf2.raw$Row.names
raw.irf2.geneList         <- sort(raw.irf2.geneList , decreasing = TRUE)

j774.irf2.geneList         <- shIrf2.j774$av_stat
names(j774.irf2.geneList)  <- shIrf2.j774$Row.names
j774.irf2.geneList         <- sort(j774.irf2.geneList , decreasing = TRUE)

######################## B. Prep gene Sets ########################
# Load Irf2 target genes
my.TF.gs <- list( "Irf2_targets_5kb_CUTRUN"     = setdiff(Irf2_5kb_perimac   ,"") ,
                  "Irf2_targets_5kb_ChIP"       = setdiff(Irf2_5kb_bmdms     ,"") ,
                  "Irf2_targets_5kb_BOTH"       = setdiff(intersect(Irf2_5kb_perimac ,Irf2_5kb_bmdms)  ,"") )


######################## C. Gene Set Enrichment Analysis ########################
# set seed to stabilize output
set.seed(1234567890)

# run phenotest GSEA
gsea.raw <- gsea( x         =  raw.irf2.geneList   ,
                  gsets     =  my.TF.gs           ,
                  mc.cores  =  1                  ,
                  logScale  =  FALSE              ,
                  B         =  10000              ,
                  minGenes  =  5                  ,
                  maxGenes  =  15000               )
my.summary.raw <- data.frame(summary(gsea.raw))
gsea.raw$significance$summary
#                           n        es      nes pval.es pval.nes fdr
# Irf2_targets_5kb_CUTRUN 478 0.5553814 2.510939       0        0   0
# Irf2_targets_5kb_ChIP   189 0.6525489 2.655458       0        0   0
# Irf2_targets_5kb_BOTH   139 0.7057130 2.743788       0        0   0

gsea.j774 <- gsea( x         =  j774.irf2.geneList   ,
                   gsets     =  my.TF.gs           ,
                   mc.cores  =  1                  ,
                   logScale  =  FALSE              ,
                   B         =  10000              ,
                   minGenes  =  5                  ,
                   maxGenes  =  15000               )
my.summary.j774 <- data.frame(summary(gsea.j774))
gsea.j774$significance$summary
#                           n        es      nes      pval.es     pval.nes          fdr
# Irf2_targets_5kb_CUTRUN 478 0.3367518 1.633147 2.220446e-16 1.255721e-03 1.255721e-03
# Irf2_targets_5kb_ChIP   189 0.4254962 1.877231 1.976197e-14 2.220446e-16 3.330669e-16
# Irf2_targets_5kb_BOTH   139 0.4495599 1.889950 6.777456e-11 2.220446e-16 6.661338e-16

############################################################################################
# Make bubble chart summary
library(ggplot2) 
library(scales) 
theme_set(theme_bw())

############## Plot TFs
my.summary.raw $condition <- "shIrf2_J264"
my.summary.j774$condition <- "shIrf2_J774"
my.summary.raw $GeneSet <- rownames(my.summary.raw )
my.summary.j774$GeneSet <- rownames(my.summary.j774)

gsea.summary    <- rbind(my.summary.raw,my.summary.j774)

# get merged datafame for ggplot
gsea.summary$minusLog10FDR <- -log10(gsea.summary$fdr + 1e-30) ### to avoid -Inf

my.max <-  3
my.min <- -3
my.values <- c(my.min,0.75*my.min,0.5*my.min,0.25*my.min,0,0.25*my.max,0.5*my.max,0.75*my.max,my.max)
my.scaled <- rescale(my.values, to = c(0, 1))
my.color.vector <- c("peru","goldenrod3","goldenrod2","goldenrod1","white","seagreen1","seagreen2","seagreen3","seagreen")

# to preserve the wanted order
gsea.summary$condition <- factor(gsea.summary$condition, levels = unique(gsea.summary$condition))
gsea.summary$GeneSet   <- factor(gsea.summary$GeneSet, levels = rev(unique(gsea.summary$GeneSet)))

pdf(paste0(Sys.Date(),"shIRF2_GSEA_Irf2_MSPCs_CUTandRUN_ChIP_5kb_GeneSets.pdf"),height = 3.5, width=6)
my.plot <- ggplot(gsea.summary,aes(x=condition,y=GeneSet,colour=nes,size=minusLog10FDR))+ theme_bw()+ geom_point(shape = 16)
my.plot <- ggplot(gsea.summary,aes(x=condition,y=GeneSet,colour=nes,size=minusLog10FDR))+ theme(text = element_text(size=16))+ geom_point(shape = 16)
my.plot <- my.plot + ggtitle("GSEA") + labs(x = "", y = "Irf2 target GeneSet")
my.plot <- my.plot + scale_colour_gradientn(colours = my.color.vector,space = "Lab", na.value = "grey50", guide = "colourbar", values = my.scaled, limits = c(my.min,my.max))
my.plot <- my.plot + scale_size_area(limits = c(1,30))
print(my.plot)
dev.off()


#######################
sink(file = paste(Sys.Date(),"_shIrf2_Irf2_binding_session_Info.txt", sep =""))
sessionInfo()
sink()


