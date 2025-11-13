setwd('/Volumes/BB_HQ_1/Immune_sex_dimorphism_Aging/Macrophages/Clean_Macrophage_Data/RNAseq/Functional_enrichment/Top_process_for_validation/Phagocytosis_genes/noSVA')
options(stringsAsFactors = F)

library(DESeq2)
library(phenoTest)

######################## A. Load DEseq2 results for analysis ######################## 
# load DEseq2 results
load('../../../../DEseq2_analysis/DEseq2/2022-03-30_Peritoneal_Macrophages_RNA_seq_results.RData')
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
# load gene lists
my.bassik.phago  <- read.table('/Volumes/BB_HQ_1//Immune_sex_dimorphism_Aging/Macrophages/Macrophage_datasets_for_Comparison/Mike_Bassik_Phagocytosis/2021-09-09_Mouse_homologs_Ens104.txt', sep = "\t", header = T)
my.wyss.phago    <- read.table('/Volumes/BB_HQ_1/Immune_sex_dimorphism_Aging/Macrophages/Macrophage_datasets_for_Comparison/Wyss_Corray_Screen/2022-01-07_Biomart_ens105_Name_conversion_FDR5_hits_BV2_screen_286.txt', sep = "\t", header = T)

my.mph.curated.gs <- list("Bassik_Screen"       = unique(my.bassik.phago$Mouse.gene.name),
                          "Wyss_Screen"         = unique(my.wyss.phago$Gene.name),
                          "combined_CRISPR"     = unique(c(my.bassik.phago$Mouse.gene.name,my.wyss.phago$Gene.name))  )

######################## C. Gene Set Enrichment Analysis ######################## 

# set seed to stabilize output (add 2020/5/20)
set.seed(123456789)

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
#                    n         es       nes     pval.es    pval.nes         fdr
# Bassik_Screen   315 -0.2753039 -1.350637 0.005280782 0.007162912 0.010744368
# Wyss_Screen     286 -0.2772437 -1.333924 0.017301554 0.009134388 0.009134388
# combined_CRISPR 580 -0.2716229 -1.397647 0.001829401 0.003845479 0.011536436


gsea.data.M <- gsea( x         =  mph.aging.M.geneList , 
                     gsets     =  my.mph.curated.gs, 
                     mc.cores  =  2                 , 
                     logScale  =  FALSE             , 
                     B         =  10000              , 
                     minGenes  =  5                 , 
                     maxGenes  =  5000               )

my.summary.M <- data.frame(summary(gsea.data.M))
my.summary.M
#                   n        es      nes     pval.es   pval.nes        fdr
# Bassik_Screen   315 0.2816259 1.357785 0.014613372 0.01159193 0.01738789
# Wyss_Screen     286 0.2870096 1.357025 0.018194950 0.01173207 0.01173207
# combined_CRISPR 580 0.2744515 1.393783 0.001203164 0.00644685 0.01934055


# write results to file
my.outfile <- paste(Sys.Date(), "Macrophage_phagocytosis_Female_GSEA_Analysis_table.txt", sep = "_")
write.table(my.summary.F, file = my.outfile, quote = F, sep = "\t")

my.outfile <- paste(Sys.Date(), "Macrophage_phagocytosis_Male_GSEA_Analysis_table.txt", sep = "_")
write.table(my.summary.M, file = my.outfile, quote = F, sep = "\t")

pdf(paste(Sys.Date(), "Macrophage_Bassik_List_GSEA_plot_FEMALE_aging.pdf", sep = "_"))
plot.gseaData(gsea.data.F, es.nes='nes', selGsets='Bassik_Screen', color = "purple")
dev.off()
pdf(paste(Sys.Date(), "Macrophage_Wyss_List_GSEA_plot_FEMALE_aging.pdf", sep = "_"))
plot.gseaData(gsea.data.F, es.nes='nes', selGsets='Wyss_Screen', color = "purple")
dev.off()
pdf(paste(Sys.Date(), "Macrophage_combinedCRISPR_List_GSEA_plot_FEMALE_aging.pdf", sep = "_"))
plot.gseaData(gsea.data.F, es.nes='nes', selGsets='combined_CRISPR', color = "purple")
dev.off()

pdf(paste(Sys.Date(), "Macrophage_Bassik_List_GSEA_plot_MALE_aging.pdf", sep = "_"))
plot.gseaData(gsea.data.M, es.nes='nes', selGsets='Bassik_Screen', color = "purple")
dev.off()
pdf(paste(Sys.Date(), "Macrophage_Wyss_List_GSEA_plot_MALE_aging.pdf", sep = "_"))
plot.gseaData(gsea.data.M, es.nes='nes', selGsets='Wyss_Screen', color = "purple")
dev.off()
pdf(paste(Sys.Date(), "Macrophage_combinedCRISPR_List_GSEA_plot_MALE_aging.pdf", sep = "_"))
plot.gseaData(gsea.data.M, es.nes='nes', selGsets='combined_CRISPR', color = "purple")
dev.off()
############################################################################################


############################################################################################
# Make bubble chart summary
library(ggplot2) 
library(scales) 
theme_set(theme_bw())


my.f.phago <- cbind(rownames(my.summary.F),my.summary.F)
my.m.phago <- cbind(rownames(my.summary.M),my.summary.M)

colnames(my.f.phago)[1] <- "GeneSet"
colnames(my.m.phago)[1] <- "GeneSet"

# get merged datafame for ggplot
merged.gsea <- rbind(my.f.phago,my.m.phago)
merged.gsea$minusLog10FDR <- -log10(merged.gsea$fdr)
merged.gsea$condition <- c(rep("Females",nrow(my.f.phago)),rep("Males",nrow(my.m.phago)))

my.max <- 2
my.min <- -2
my.values <- c(my.min,0.75*my.min,0.5*my.min,0.25*my.min,0,0.25*my.max,0.5*my.max,0.75*my.max,my.max)
my.scaled <- rescale(my.values, to = c(0, 1))
my.color.vector <- c("darkblue","dodgerblue4","dodgerblue3","dodgerblue1","white","lightcoral","brown1","firebrick2","firebrick4")

# to preserve the wanted order
merged.gsea$condition <- factor(merged.gsea$condition, levels = unique(merged.gsea$condition))
merged.gsea$GeneSet  <- factor(merged.gsea$GeneSet, levels = rev(unique(merged.gsea$GeneSet)))

pdf(paste0(Sys.Date(),"CRISPR_Screen_GSEA_Male_Female_Aging_Perimac_CLEAN.pdf"),height = 4, width=6)
my.plot <- ggplot(merged.gsea,aes(x=condition,y=GeneSet,colour=nes,size=minusLog10FDR))+ theme_bw()+ geom_point(shape = 16)
my.plot <- ggplot(merged.gsea,aes(x=condition,y=GeneSet,colour=nes,size=minusLog10FDR))+ theme(text = element_text(size=16))+ geom_point(shape = 16)
my.plot <- my.plot + ggtitle("GSEA") + labs(x = "PeriMac Aging", y = "Phagocytosis regulator GeneSet")
my.plot <- my.plot + scale_colour_gradientn(colours = my.color.vector,space = "Lab", na.value = "grey50", guide = "colourbar", values = my.scaled, limits = c(-3, 3))
my.plot <- my.plot + scale_size_area(limits = c(0,5))
print(my.plot)
dev.off()  




#######################
sink(file = paste(Sys.Date(),"_session_Info.txt", sep =""))
sessionInfo()
sink()



