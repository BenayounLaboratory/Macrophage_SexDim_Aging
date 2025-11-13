setwd('/Volumes/BB_HQ_1/Immune_sex_dimorphism_Aging/Macrophages/Clean_Macrophage_Data/RNAseq/Functional_enrichment/XCIescape_Xist')
options(stringsAsFactors = F)

library(DESeq2)
library(phenoTest)

######################## A. Load DEseq2 results for analysis ######################## 
# load DEseq2 results
load('../../DEseq2_analysis/DEseq2/2022-03-30_Peritoneal_Macrophages_RNA_seq_results.RData')
my.mph1.RNAseq.process
# [1] "Sex"             "Aging"           "Female_Aging"    "Male_Aging"      "FemaleMaleAging"

my.mph.aging.F           <- data.frame(my.mph1.RNAseq.process$Female_Aging)
my.mph.aging.F$GeneName  <- rownames(my.mph.aging.F )


######################## B. Prepare GeneLists using DEseq2 t-statistic to rank genes ######################## 
mph.aging.F.geneList         <- my.mph.aging.F$stat
names(mph.aging.F.geneList)  <- my.mph.aging.F$GeneName
mph.aging.F.geneList         <- sort(mph.aging.F.geneList , decreasing = TRUE)


######################## C. Prep X/y-link gene Sets ######################## 
# load gene lists
my.XCI.escape      <- read.table('/Volumes/BB_HQ_1/PATHWAY_ANNOT/X_inactivation_escape/2023-11-01_updated_known_mouse_XCI_genes.txt', sep = "\t", header = F)
my.mm10.genes      <- read.csv('../../OmicsCircos/2019-06-12_HOMER_annotated_TSS_mm10.xls', sep = "\t", header = T)
my.x.genes         <- unique(my.mm10.genes$Gene.Name[my.mm10.genes$Chr == "chrX"])

my.mph.curated.gs <- list("XCI_escapees"     = unique(my.XCI.escape$V1),
                          "X_linked_genes"   = my.x.genes)

######################## D. Gene Set Enrichment Analysis ######################## 

# set seed to stabilize output
set.seed(123456789)

# run phenotest GSEA
######################################
gsea.data.F <- gsea( x         =  mph.aging.F.geneList , 
                     gsets     =  my.mph.curated.gs, 
                     mc.cores  =  2                 , 
                     logScale  =  FALSE             , 
                     B         =  10000              , 
                     minGenes  =  5                 , 
                     maxGenes  =  5000               )

my.summary.F <- data.frame(summary(gsea.data.F))
my.summary.F
#                   n         es       nes     pval.es    pval.nes         fdr
# XCI_escapees    250 -0.3040105 -1.424572 0.003568083 0.001795507 0.003591014
# X_linked_genes 1069 -0.2538431 -1.299561 0.003231949 0.014667380 0.014667380

# write results to file
my.outfile <- paste(Sys.Date(), "Macrophage_Xlinked_genelists_Female_GSEA_Analysis_table.txt", sep = "_")
write.table(my.summary.F, file = my.outfile, quote = F, sep = "\t")
############################################################################################


############################################################################################
# Make bubble chart summary
library(ggplot2) 
library(scales) 
theme_set(theme_bw())

############## 
my.xlist <- cbind(rownames(my.summary.F),my.summary.F)
colnames(my.xlist)[1] <- "GeneSet"


# get merged datafame for ggplot
my.xlist$minusLog10FDR <- -log10(my.xlist$fdr)
my.xlist$condition <- rep("F_aging",nrow(my.xlist))

my.max <- 2.5
my.min <- -2
my.values <- c(my.min,0.75*my.min,0.5*my.min,0.25*my.min,0,0.25*my.max,0.5*my.max,0.75*my.max,my.max)
my.scaled <- rescale(my.values, to = c(0, 1))
my.color.vector <- c("darkblue","dodgerblue4","dodgerblue3","dodgerblue1","white","lightcoral","brown1","firebrick2","firebrick4")

# to preserve the wanted order
my.xlist$condition <- factor(my.xlist$condition, levels = unique(my.xlist$condition))
my.xlist$GeneSet   <- factor(my.xlist$GeneSet, levels = rev(unique(my.xlist$GeneSet)))

pdf(paste0(Sys.Date(),"F_Aging_GSEA_X_related_GeneSets.pdf"),height = 4, width=5)
my.plot <- ggplot(my.xlist,aes(x=condition,y=GeneSet,colour=nes,size=minusLog10FDR))+ theme_bw()+ geom_point(shape = 16)
my.plot <- ggplot(my.xlist,aes(x=condition,y=GeneSet,colour=nes,size=minusLog10FDR))+ theme(text = element_text(size=16))+ geom_point(shape = 16)
my.plot <- my.plot + ggtitle("GSEA") + labs(x = "F aging", y = "X related GeneSet")
my.plot <- my.plot + scale_colour_gradientn(colours = my.color.vector,space = "Lab", na.value = "grey50", guide = "colourbar", values = my.scaled, limits = c(my.min,my.max))
my.plot <- my.plot + scale_size_area(limits = c(1,4))
print(my.plot)
dev.off()


############################################################################################


#######################
sink(file = paste(Sys.Date(),"_session_Info.txt", sep =""))
sessionInfo()
sink()



