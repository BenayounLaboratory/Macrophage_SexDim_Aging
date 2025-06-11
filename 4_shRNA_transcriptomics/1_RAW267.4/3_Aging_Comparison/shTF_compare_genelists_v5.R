setwd('/Volumes/BB_HQ_1/Immune_sex_dimorphism_Aging/Macrophages/Clean_Macrophage_Data/RNAseq/Functional_enrichment/Top_TF_Validation/shTF_target_analysis')
options(stringsAsFactors = F)

# load libraries for analysis
library(DESeq2)
library(Vennerable)
library(bitops)
library(phenoTest)
library(ggplot2) 
library(scales) 

# 2023-06-30
# Compare shRNA and female aging
# Use consistent sets/combined p

# 2023-10-11
# unstranded run, add Tal1

# 2024-05-23
# run with new shTF data

# 2024-08-16
# run with all shTF data
# Each TF was processed on its own for reproducibility

# 2024-08-22
# use correct Mef2c data (previous code had bug on DE genes)

################################################################################
# 0. Load shTF responsive genes
irf.sig <- read.table('/Volumes/BB_HQ_1/Immune_sex_dimorphism_Aging/Macrophages/Clean_Macrophage_Data/Intervention_Data/shRNA_analyses/DEseq2/2024-08-16_RAW_sh_Irf2_data_DESeq2_shIrf2_DEseq2_results_COMBINED_FDR1e-4_genes_statistics.txt', header = T)
irf.up  <- irf.sig$Row.names[irf.sig$log2FoldChange.sh1  > 0]
irf.dwn <- irf.sig$Row.names[irf.sig$log2FoldChange.sh1 < 0]

tbl.sig <- read.table('/Volumes/BB_HQ_1/Immune_sex_dimorphism_Aging/Macrophages/Clean_Macrophage_Data/Intervention_Data/shRNA_analyses/DEseq2/2024-08-16_RAW_sh_Tbl1xr1_data_DESeq2_shTbl1xr1_DEseq2_results_COMBINED_FDR1e-4_genes_statistics.txt', header = T)
tbl.up  <- tbl.sig$Row.names[tbl.sig$log2FoldChange.sh1 > 0]
tbl.dwn <- tbl.sig$Row.names[tbl.sig$log2FoldChange.sh1 < 0]

tal.sig <- read.table('/Volumes/BB_HQ_1/Immune_sex_dimorphism_Aging/Macrophages/Clean_Macrophage_Data/Intervention_Data/shRNA_analyses/DEseq2/2024-08-16_RAW_sh_Tal1_data_DESeq2_shTal1_DEseq2_results_COMBINED_FDR1e-4_genes_statistics.txt', header = T)
tal.up  <- tal.sig$Row.names[tal.sig$log2FoldChange.sh1 > 0]
tal.dwn <- tal.sig$Row.names[tal.sig$log2FoldChange.sh1 < 0]

meis1.sig <- read.table('/Volumes/BB_HQ_1/Immune_sex_dimorphism_Aging/Macrophages/Clean_Macrophage_Data/Intervention_Data/shRNA_analyses/DEseq2/2024-08-16_RAW_sh_Meis1_data_DESeq2_shMeis1_DEseq2_results_COMBINED_FDR1e-4_genes_statistics.txt', header = T)
meis1.up  <- meis1.sig$Row.names[meis1.sig$log2FoldChange.sh1 > 0]
meis1.dwn <- meis1.sig$Row.names[meis1.sig$log2FoldChange.sh1 < 0]

mef2c.sig <- read.table('/Volumes/BB_HQ_1/Immune_sex_dimorphism_Aging/Macrophages/Clean_Macrophage_Data/Intervention_Data/shRNA_analyses/DEseq2/2024-08-22_RAW_sh_Mef2c_data_DESeq2_shMef2c_DEseq2_results_COMBINED_FDR1e-4_genes_statistics.txt', header = T)
mef2c.up  <- mef2c.sig$Row.names[mef2c.sig$log2FoldChange.sh1 > 0]
mef2c.dwn <- mef2c.sig$Row.names[mef2c.sig$log2FoldChange.sh1 < 0]

# Load female all aging
fem.age.all <- read.table('../../../DEseq2_analysis/DESeq2/2022-03-30_Peritoneal_Macrophages_AGING_Female_ALL_genes_statistics.txt')
fem.all.up       <- rownames(fem.age.all)[bitAnd(fem.age.all$log2FoldChange_F > 0, fem.age.all$padj_F < 0.05) >0]
fem.all.dwn      <- rownames(fem.age.all)[bitAnd(fem.age.all$log2FoldChange_F < 0, fem.age.all$padj_F < 0.05) >0]


################################################################################
# 1. Run Gene Set enrichment analysis

######################## A. format DESeq2 results  ########################
res.fem           <- data.frame(fem.age.all)
res.fem$GeneName  <- rownames(fem.age.all )

######################## B. Prepare GeneLists using DEseq2 t-statistic to rank genes ########################
res.fem.geneList         <- res.fem$stat
names(res.fem.geneList)  <- res.fem$GeneName
res.fem.geneList         <- sort(res.fem.geneList , decreasing = TRUE)

######################## C. Prep gene Sets ########################
# load gene lists

# Load IRF, TBL and TAL responsive genes
my.irf.gs <- list("shIrf2_up"          = irf.up ,
                  "shIrf2_dwn"         = irf.dwn)
   
my.tbl.gs <- list("shTbl1xr1_up"       = tbl.up ,
                  "shTbl1xr1_dwn"      = tbl.dwn)
   
my.tal1.gs <- list("shTal1_up"          = tal.up ,
                   "shTal1_dwn"         = tal.dwn)

my.meis1.gs <- list("shMeis1_up"       = meis1.up ,
                    "shMeis1_dwn"      = meis1.dwn )

my.mef2c.gs <- list("shMef2c_up"       = mef2c.up ,
                    "shMef2c_dwn"      = mef2c.dwn )

######################## D. Gene Set Enrichment Analysis ########################
# set seed to stabilize output
set.seed(1234567890)

# run phenotest GSEA
gsea.res.irf2    <- gsea( x =  res.fem.geneList  , gsets     =  my.irf.gs   , mc.cores  =  1 , logScale  =  FALSE , B =  10000 , minGenes  =  5 , maxGenes  =  5000 )
gsea.res.tbl1xr1 <- gsea( x =  res.fem.geneList  , gsets     =  my.tbl.gs   , mc.cores  =  1 , logScale  =  FALSE , B =  10000 , minGenes  =  5 , maxGenes  =  5000 )
gsea.res.tal1    <- gsea( x =  res.fem.geneList  , gsets     =  my.tal1.gs  , mc.cores  =  1 , logScale  =  FALSE , B =  10000 , minGenes  =  5 , maxGenes  =  5000 )
gsea.res.meis1   <- gsea( x =  res.fem.geneList  , gsets     =  my.meis1.gs , mc.cores  =  1 , logScale  =  FALSE , B =  10000 , minGenes  =  5 , maxGenes  =  5000 )
gsea.res.mef2c   <- gsea( x =  res.fem.geneList  , gsets     =  my.mef2c.gs , mc.cores  =  1 , logScale  =  FALSE , B =  10000 , minGenes  =  5 , maxGenes  =  5000 )

my.summary.irf2    <- data.frame(gsea.res.irf2   $significance$summary)
my.summary.tbl1xr1 <- data.frame(gsea.res.tbl1xr1$significance$summary)
my.summary.tal1    <- data.frame(gsea.res.tal1   $significance$summary)
my.summary.meis1   <- data.frame(gsea.res.meis1  $significance$summary)
my.summary.mef2c   <- data.frame(gsea.res.mef2c  $significance$summary)

# n         es       nes      pval.es     pval.nes          fdr
# shIrf2_up  1011  0.2868500  1.466009 6.931422e-11 3.501643e-13 3.501643e-13
# shIrf2_dwn 1046 -0.2359014 -1.279573 3.755922e-04 2.095153e-04 2.095153e-04
# n         es        nes      pval.es     pval.nes          fdr
# shTbl1xr1_up  564 -0.1898773 -0.9862234 5.100182e-01 5.183325e-01 5.183325e-01
# shTbl1xr1_dwn 696  0.3010994  1.5088945 1.780798e-13 2.735893e-09 2.735893e-09
# n         es       nes   pval.es  pval.nes       fdr
# shTal1_up  316  0.2206010  1.026784 0.3823784 0.3809383 0.3809383
# shTal1_dwn 421 -0.2040501 -1.030401 0.3481647 0.3503116 0.3503116
# n         es       nes      pval.es     pval.nes          fdr
# shMeis1_up  593 -0.2114267 -1.097197 1.400528e-01 1.376368e-01 1.376368e-01
# shMeis1_dwn 628 -0.3525580 -1.849613 2.220446e-16 2.220446e-16 4.440892e-16
# n         es       nes      pval.es     pval.nes          fdr
# shMef2c_up  1004 -0.2072920 -1.116918 0.06446197 0.05549102 0.05549102
# shMef2c_dwn 1155 -0.2939089 -1.602635 0.00000000 0.00000000 0.00000000

################################################################################
# 2.  Make bubble chart summary

theme_set(theme_bw())

############## Plot results
my.summary.irf2   $TF <- "shIrf2"
my.summary.tbl1xr1$TF <- "shTbl1xr1"
my.summary.tal1   $TF <- "shTal1"
my.summary.meis1  $TF <- "shMeis1"
my.summary.mef2c  $TF <- "shMef2c"

my.summary.irf2   $GeneSet <- c("shRNA_up","shRNA_down")
my.summary.tbl1xr1$GeneSet <- c("shRNA_up","shRNA_down")
my.summary.tal1   $GeneSet <- c("shRNA_up","shRNA_down")
my.summary.meis1  $GeneSet <- c("shRNA_up","shRNA_down")
my.summary.mef2c  $GeneSet <- c("shRNA_up","shRNA_down")

# get merged datafame for ggplot
my.merged.TF.screen <- rbind(my.summary.irf2   ,
                             my.summary.tbl1xr1,
                             my.summary.tal1   ,
                             my.summary.meis1  ,
                             my.summary.mef2c  )

# get merged datafame for ggplot
my.merged.TF.screen$minusLog10FDR <- -log10(my.merged.TF.screen$fdr + 1e-20)

my.max <- 2
my.min <- -2
my.values <- c(my.min,0.75*my.min,0.5*my.min,0.25*my.min,0,0.25*my.max,0.5*my.max,0.75*my.max,my.max)
my.scaled <- rescale(my.values, to = c(0, 1))
my.color.vector <- c("darkblue","dodgerblue4","dodgerblue3","dodgerblue1","white","lightcoral","brown1","firebrick2","firebrick4")

# to preserve the wanted order
my.merged.TF.screen$TF        <- factor(my.merged.TF.screen$TF , levels = c("shIrf2","shMef2c","shMeis1","shTal1", "shTbl1xr1"))
my.merged.TF.screen$GeneSet   <- factor(my.merged.TF.screen$GeneSet, levels = rev(unique(my.merged.TF.screen$GeneSet)))

pdf(paste0(Sys.Date(),"_F_Aging_GSEA_shTF_RNAseq_Screen_GeneSets_CLEAN.pdf"),height = 3, width=7)
my.plot <- ggplot(my.merged.TF.screen,aes(x=TF,y=GeneSet,colour=nes,size=minusLog10FDR))+ theme_bw()+ geom_point(shape = 16)
my.plot <- ggplot(my.merged.TF.screen,aes(x=TF,y=GeneSet,colour=nes,size=minusLog10FDR))+ theme(text = element_text(size=16))+ geom_point(shape = 16)
my.plot <- my.plot + ggtitle("GSEA (Female aging)") + labs(x = "F aging", y = "shRNA target GeneSet")
my.plot <- my.plot + scale_colour_gradientn(colours = my.color.vector,space = "Lab", na.value = "grey50", guide = "colourbar", values = my.scaled, limits = c(my.min,my.max))
my.plot <- my.plot + scale_size_area(limits = c(0,20)) + scale_x_discrete(guide = guide_axis(angle = 45)) + scale_size(range = c(2,8))
print(my.plot)
dev.off()


#######################
sink(file = paste(Sys.Date(),"_F_aging_shRNA_session_Info.txt", sep =""))
sessionInfo()
sink()


