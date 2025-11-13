setwd('/Volumes/BB_HQ_1/Immune_sex_dimorphism_Aging/Macrophages/Clean_Macrophage_Data/RNAseq/Functional_enrichment/Top_process_for_validation/Polarization/')
options(stringsAsFactors = F)

library(DESeq2)
library(phenoTest)
library(qusage)

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


######################## C. Prep Neutrophil gene Sets ######################## 
# load gene lists polarization
Polar.gsets  <- read.gmt("/Volumes/BB_HQ_1/Immune_sex_dimorphism_Aging/Macrophages/Macrophage_datasets_for_Comparison/M1_M2_macrophages/GSE103958_Bulk/DEseq2/2021-11-18_GSE103958_BMDM_polarization_genes_lists.gmt"              )

# format for GSEA
my.mph.curated.gs <- c(Polar.gsets)


######################## C. Gene Set Enrichment Analysis ######################## 

# set seed to stabilize output 
set.seed(987654321)

# run phenotest GSEA
gsea.data.F <- gsea( x         =  mph.aging.F.geneList , 
                     gsets     =  my.mph.curated.gs, 
                     mc.cores  =  2                 , 
                     logScale  =  FALSE             , 
                     B         =  10000              , 
                     minGenes  =  5                 , 
                     maxGenes  =  5000               )

my.summary.F <- data.frame(summary(gsea.data.F))
my.summary.F
#                     n         es       nes      pval.es     pval.nes          fdr
# BMDM_M1_enriched 1392  0.3522367  1.836844 2.220446e-16 2.220446e-16 2.220446e-16
# BMDM_M2_enriched  760 -0.2224965 -1.180919 2.807181e-02 1.809302e-02 1.809302e-02

gsea.data.M <- gsea( x         =  mph.aging.M.geneList , 
                     gsets     =  my.mph.curated.gs, 
                     mc.cores  =  2                 , 
                     logScale  =  FALSE             , 
                     B         =  10000              , 
                     minGenes  =  5                 , 
                     maxGenes  =  5000               )

my.summary.M <- data.frame(summary(gsea.data.M))
my.summary.M
#                     n        es      nes      pval.es     pval.nes          fdr
# BMDM_M1_enriched 1392 0.4766482 2.536113 0.000000e+00 2.220446e-16 4.440892e-16
# BMDM_M2_enriched  760 0.4511345 2.339842 2.220446e-16 0.000000e+00 0.000000e+00

# write results to file
my.outfile <- paste(Sys.Date(), "Macrophage_Polarization_Female_GSEA_Analysis_table.txt", sep = "_")
write.table(my.summary.F, file = my.outfile, quote = F, sep = "\t")

my.outfile <- paste(Sys.Date(), "Macrophage_Polarization_Male_GSEA_Analysis_table.txt", sep = "_")
write.table(my.summary.M, file = my.outfile, quote = F, sep = "\t")

pdf(paste(Sys.Date(), "Polarization_List_GSEA_plot_FEMALE_aging.pdf", sep = "_"))
plot.gseaData(gsea.data.F, es.nes='nes', selGsets='BMDM_M1_enriched', color = "purple")
plot.gseaData(gsea.data.F, es.nes='nes', selGsets='BMDM_M2_enriched', color = "purple")
dev.off()

pdf(paste(Sys.Date(), "Polarization_List_GSEA_plot_MALE_aging.pdf", sep = "_"))
plot.gseaData(gsea.data.M, es.nes='nes', selGsets='BMDM_M1_enriched', color = "purple")
plot.gseaData(gsea.data.M, es.nes='nes', selGsets='BMDM_M2_enriched', color = "purple")
dev.off()
############################################################################################


############################################################################################
# Make bubble chart summary
library(ggplot2) 
library(scales) 
theme_set(theme_bw())


my.f.polar <- cbind(rownames(my.summary.F),my.summary.F)
my.m.polar <- cbind(rownames(my.summary.M),my.summary.M)

colnames(my.f.polar)[1] <- "GeneSet"
colnames(my.m.polar)[1] <- "GeneSet"

# get merged datafame for ggplot
merged.gsea <- rbind(my.f.polar,my.m.polar)
merged.gsea$minusLog10FDR <- -log10(merged.gsea$fdr + 1e-20)
merged.gsea$condition <- c(rep("Females",nrow(my.f.polar)),rep("Males",nrow(my.m.polar)))

my.max <- 2
my.min <- -2
my.values <- c(my.min,0.75*my.min,0.5*my.min,0.25*my.min,0,0.25*my.max,0.5*my.max,0.75*my.max,my.max)
my.scaled <- rescale(my.values, to = c(0, 1))
my.color.vector <- c("darkblue","dodgerblue4","dodgerblue3","dodgerblue1","white","lightcoral","brown1","firebrick2","firebrick4")

# to preserve the wanted order
merged.gsea$condition <- factor(merged.gsea$condition, levels = unique(merged.gsea$condition))
merged.gsea$GeneSet  <- factor(merged.gsea$GeneSet, levels = rev(unique(merged.gsea$GeneSet)))

pdf(paste0(Sys.Date(),"Polarization_GSEA_Male_Female_Aging_Perimac_CLEAN.pdf"),height = 3, width=6)
my.plot <- ggplot(merged.gsea,aes(x=condition,y=GeneSet,colour=nes,size=minusLog10FDR))+ theme_bw()+ geom_point(shape = 16)
my.plot <- ggplot(merged.gsea,aes(x=condition,y=GeneSet,colour=nes,size=minusLog10FDR))+ theme(text = element_text(size=16))+ geom_point(shape = 16)
my.plot <- my.plot + ggtitle("GSEA") + labs(x = "PeriMac Aging", y = "Polarization GeneSets")
my.plot <- my.plot + scale_colour_gradientn(colours = my.color.vector,space = "Lab", na.value = "grey50", guide = "colourbar", values = my.scaled, limits = c(-3, 3))
my.plot <- my.plot + scale_size_area(limits = c(0,20))
print(my.plot)
dev.off()  



#######################
sink(file = paste(Sys.Date(),"_session_Info.txt", sep =""))
sessionInfo()
sink()



