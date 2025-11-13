setwd('/Volumes/BB_HQ_1/Immune_sex_dimorphism_Aging/Macrophages/Clean_Macrophage_Data/RNAseq/Functional_enrichment/Top_process_for_validation/Metabolism_genes/')
options(stringsAsFactors = F)

library(DESeq2)
library(phenoTest)


# 2024-08-09
# Analyze metabolism gene lists
# concentrate on curated HALLMARK lists

######################## A. Load DEseq2 results for analysis ######################## 
# load DEseq2 results
load('../../../DEseq2_analysis/DEseq2/2022-03-30_Peritoneal_Macrophages_RNA_seq_results.RData')
my.mph1.RNAseq.process
# [1] "Sex"             "Aging"           "Female_Aging"    "Male_Aging"      "FemaleMaleAging"

my.mph.aging.F           <- data.frame(my.mph1.RNAseq.process$Female_Aging)
my.mph.aging.F$GeneName  <- rownames(my.mph.aging.F )

my.mph.aging.M           <- data.frame(my.mph1.RNAseq.process$Male_Aging)
my.mph.aging.M$GeneName  <- rownames(my.mph.aging.M )

######################## B. Prepare GeneLists using DEseq2 t-statistic to rank genes ######################## 
mph.aging.F.geneList         <- my.mph.aging.F$stat
names(mph.aging.F.geneList)  <- my.mph.aging.F$GeneName
mph.aging.F.geneList         <- sort(mph.aging.F.geneList , decreasing = TRUE)

mph.aging.M.geneList         <- my.mph.aging.M$stat
names(mph.aging.M.geneList)  <- my.mph.aging.M$GeneName
mph.aging.M.geneList         <- sort(mph.aging.M.geneList , decreasing = TRUE)


######################## C. Prep gene Sets ########################

###### HALLMARK
hall.glyco <- read.gmt('HALLMARK_GLYCOLYSIS.v2023.2.Mm.gmt')
hall.oxphos <- read.gmt('HALLMARK_OXIDATIVE_PHOSPHORYLATION.v2023.2.Mm.gmt')

# load gene lists
my.glyc.gs <- list("Hallmark_Glycolysis"                 = unique(hall.glyco$HALLMARK_GLYCOLYSIS) )
my.oxphos.gs <- list("Hallmark_Oxidative_Phosphorylation"         = unique(hall.oxphos$HALLMARK_OXIDATIVE_PHOSPHORYLATION) )

metab.gs <- list("Hallmark_Glycolysis"                 = unique(hall.glyco$HALLMARK_GLYCOLYSIS) ,
                 "Hallmark_Oxidative_Phosphorylation"         = unique(hall.oxphos$HALLMARK_OXIDATIVE_PHOSPHORYLATION))

######################## C. Gene Set Enrichment Analysis ######################## 

# set seed to stabilize output
set.seed(123456789)

# run phenotest GSEA
gsea.data.F <- gsea( x         =  mph.aging.F.geneList , 
                     gsets     =  metab.gs, 
                     mc.cores  =  2                 , 
                     logScale  =  FALSE             , 
                     B         =  10000              , 
                     minGenes  =  5                 , 
                     maxGenes  =  5000               )
my.summary.F <- data.frame(summary(gsea.data.F))
my.summary.F
#                                      n         es        nes     pval.es    pval.nes         fdr
# Hallmark_Glycolysis                200 -0.3217400 -1.4728047 0.002461097 0.002625896 0.005251791
# Hallmark_Oxidative_Phosphorylation 195 -0.1699869 -0.7938202 0.949592676 0.949566620 0.949566620

gsea.data.M <- gsea( x         =  mph.aging.M.geneList , 
                     gsets     =  metab.gs, 
                     mc.cores  =  2                 , 
                     logScale  =  FALSE             , 
                     B         =  10000              , 
                     minGenes  =  5                 , 
                     maxGenes  =  5000               )
my.summary.M <- data.frame(summary(gsea.data.M))
my.summary.M
#                                      n        es      nes pval.es     pval.nes          fdr
# Hallmark_Glycolysis                200 0.4878508 2.224792       0 0.000000e+00 0.000000e+00
# Hallmark_Oxidative_Phosphorylation 195 0.6214242 2.865363       0 2.220446e-16 4.440892e-16

# write results to file
my.outfile <- paste(Sys.Date(), "Macrophage_Metabolism_HALLMARK_Female_GSEA_Analysis_table.txt", sep = "_")
write.table(my.summary.F, file = my.outfile, quote = F, sep = "\t")

my.outfile <- paste(Sys.Date(), "Macrophage_Metabolism_HALLMARK_Male_GSEA_Analysis_table.txt", sep = "_")
write.table(my.summary.M, file = my.outfile, quote = F, sep = "\t")

############################################################################################


############################################################################################
# Make bubble chart summary
library(ggplot2) 
library(scales) 
theme_set(theme_bw())

my.f.metab <- cbind(rownames(my.summary.F),my.summary.F)
my.m.metab <- cbind(rownames(my.summary.M),my.summary.M)

colnames(my.f.metab)[1] <- "GeneSet"
colnames(my.m.metab)[1] <- "GeneSet"

# get merged datafame for ggplot
merged.gsea <- rbind(my.f.metab,my.m.metab)
merged.gsea$minusLog10FDR <- -log10(merged.gsea$fdr + 1e-20)
merged.gsea$condition <- c(rep("Females",nrow(my.f.metab)),rep("Males",nrow(my.m.metab)))

my.max <- 3
my.min <- -3
my.values <- c(my.min,0.75*my.min,0.5*my.min,0.25*my.min,0,0.25*my.max,0.5*my.max,0.75*my.max,my.max)
my.scaled <- rescale(my.values, to = c(0, 1))
my.color.vector <- c("darkblue","dodgerblue4","dodgerblue3","dodgerblue1","white","lightcoral","brown1","firebrick2","firebrick4")

# to preserve the wanted order
merged.gsea$condition <- factor(merged.gsea$condition, levels = unique(merged.gsea$condition))
merged.gsea$GeneSet  <- factor(merged.gsea$GeneSet, levels = rev(unique(merged.gsea$GeneSet)))

pdf(paste0(Sys.Date(),"HALLMARK_Metabolism_GSEA_Male_Female_Aging_Perimac_CLEAN.pdf"),height = 3, width=7)
my.plot <- ggplot(merged.gsea,aes(x=condition,y=GeneSet,colour=nes,size=minusLog10FDR))+ theme_bw()+ geom_point(shape = 16)
my.plot <- ggplot(merged.gsea,aes(x=condition,y=GeneSet,colour=nes,size=minusLog10FDR))+ theme(text = element_text(size=16))+ geom_point(shape = 16)
my.plot <- my.plot + ggtitle("GSEA") + labs(x = "PeriMac Aging", y = "Phagocytosis regulator GeneSet")
my.plot <- my.plot + scale_colour_gradientn(colours = my.color.vector,space = "Lab", na.value = "grey50", guide = "colourbar", values = my.scaled, limits = c(-3, 3))
my.plot <- my.plot + scale_size_area(limits = c(0,20))
print(my.plot)
dev.off()  




#######################
sink(file = paste(Sys.Date(),"_session_Info.txt", sep =""))
sessionInfo()
sink()



