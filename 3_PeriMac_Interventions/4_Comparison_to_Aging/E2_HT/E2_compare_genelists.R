setwd('/Volumes/BB_HQ_1/Immune_sex_dimorphism_Aging/Macrophages/Clean_Macrophage_Data/RNAseq/Functional_enrichment/E2_supplementation_comparison')
options(stringsAsFactors = F)

# load libraries for analysis
library(DESeq2)
library(Vennerable)
library(bitops)

# 2023-10-24
# Compare E2 HT and female aging

# Load OVX responsive genes
e2.sig <- read.table('/Volumes/BB_HQ_1/Immune_sex_dimorphism_Aging/Macrophages/Clean_Macrophage_Data/Intervention_Data/E2_supplementation/DESeq2/2023-10-24_Peritoneal_Macrophage_Aging_E2sup_DESeq2_Old_E2_FDR5_genes_statistics.txt')
e2.up  <- rownames(e2.sig)[e2.sig$log2FoldChange > 0]
e2.dwn <- rownames(e2.sig)[e2.sig$log2FoldChange < 0]

# Load female specific aging
fem.spe.age <- read.table('../../DEseq2_analysis/DESeq2/2022-03-30_Peritoneal_Macrophages_AGING_Female.NOT.Male_FDR5_genes_statistics.txt')
fem.spe.up  <- rownames(fem.spe.age)[fem.spe.age$log2FoldChange_F > 0]
fem.spe.dwn <- rownames(fem.spe.age)[fem.spe.age$log2FoldChange_F < 0]

# Load female all aging
fem.age.all <- read.table('../../DEseq2_analysis/DESeq2/2022-03-30_Peritoneal_Macrophages_AGING_Female_ALL_genes_statistics.txt')
fem.all.up       <- rownames(fem.age.all)[bitAnd(fem.age.all$log2FoldChange_F > 0, fem.age.all$padj_F < 0.05) >0]
fem.all.dwn      <- rownames(fem.age.all)[bitAnd(fem.age.all$log2FoldChange_F < 0, fem.age.all$padj_F < 0.05) >0]

###############################################################################################################################################
library(phenoTest)

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

# Load E2 responsive genes

my.E2.gs <- list("E2_up"       = e2.up,
                 "E2_dwn"      = e2.dwn)

######################## C. Gene Set Enrichment Analysis ########################
# set seed to stabilize output
set.seed(1234567890)

# run phenotest GSEA
gsea.res.sex <- gsea( x         =  res.fem.geneList ,
                      gsets     =  my.E2.gs,
                      mc.cores  =  1                 ,
                      logScale  =  FALSE             ,
                      B         =  10000              ,
                      minGenes  =  5                 ,
                      maxGenes  =  5000               )
my.summary <- data.frame(summary(gsea.res.sex))
gsea.res.sex$significance$summary
#          n         es       nes     pval.es   pval.nes        fdr
# E2_up   15 -0.5726292 -1.570947 0.033374425 0.02005031 0.02005031
# E2_dwn 129  0.3871436  1.564133 0.007907683 0.02200035 0.02200035

pdf(paste(Sys.Date(), "E2_responsive_GSEA_plot_Impact_F_aging.pdf", sep = "_"))
plot.gseaData(gsea.res.sex, es.nes='nes', selGsets='E2_up', color = "purple")
plot.gseaData(gsea.res.sex, es.nes='nes', selGsets='E2_dwn', color = "purple")
dev.off()



############################################################################################
# Make bubble chart summary
library(ggplot2) 
library(scales) 
theme_set(theme_bw())

############## Plot OVX
gsea.summary <- data.frame(gsea.res.sex$significance$summary)
my.e2 <- cbind(rownames(gsea.summary),gsea.summary)
colnames(my.e2)[1] <- "GeneSet"


# get merged datafame for ggplot
my.e2$minusLog10FDR <- -log10(my.e2$fdr)
my.e2$condition <- rep("F_aging",nrow(my.e2))

my.max <- 2.5
my.min <- -2
my.values <- c(my.min,0.75*my.min,0.5*my.min,0.25*my.min,0,0.25*my.max,0.5*my.max,0.75*my.max,my.max)
my.scaled <- rescale(my.values, to = c(0, 1))
my.color.vector <- c("darkblue","dodgerblue4","dodgerblue3","dodgerblue1","white","lightcoral","brown1","firebrick2","firebrick4")

# to preserve the wanted order
my.e2$condition <- factor(my.e2$condition, levels = unique(my.e2$condition))
my.e2$GeneSet   <- factor(my.e2$GeneSet, levels = rev(unique(my.e2$GeneSet)))

pdf(paste0(Sys.Date(),"F_Aging_GSEA_E2sup_GeneSets.pdf"),height = 4, width=5)
my.plot <- ggplot(my.e2,aes(x=condition,y=GeneSet,colour=nes,size=minusLog10FDR))+ theme_bw()+ geom_point(shape = 16)
my.plot <- ggplot(my.e2,aes(x=condition,y=GeneSet,colour=nes,size=minusLog10FDR))+ theme(text = element_text(size=16))+ geom_point(shape = 16)
my.plot <- my.plot + ggtitle("GSEA") + labs(x = "F aging", y = "E2 responsive GeneSet")
my.plot <- my.plot + scale_colour_gradientn(colours = my.color.vector,space = "Lab", na.value = "grey50", guide = "colourbar", values = my.scaled, limits = c(my.min,my.max))
my.plot <- my.plot + scale_size_area(limits = c(1,3))
print(my.plot)
dev.off()
#######################


#######################
sink(file = paste(Sys.Date(),"_F_aging_E2_supp_session_Info.txt", sep =""))
sessionInfo()
sink()


