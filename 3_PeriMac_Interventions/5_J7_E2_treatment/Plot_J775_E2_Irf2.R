# setwd('/Volumes/BB_HQ_1/Immune_sex_dimorphism_Aging/Macrophages/E2_treaments/RNAseq/DEseq2')
setwd('/Volumes/BB_Travel_2/J774_E2/DESeq2')
options(stringsAsFactors = F)

# load libraries for analysis
library(DESeq2)
library(Vennerable)
library(bitops)
library(beeswarm)

# 2025-11-24
# Plot Irf2 expression

####################################################################################
# Load E2 response
e2.sig <- read.table('/Volumes/BB_Travel_2/J774_E2/DESeq2/2025-11-24_J774A.1_E2_treatment_DESeq2_MM10_E2_all_genes_statistics.txt')

e2.cts <- read.table('/Volumes/BB_Travel_2/J774_E2/DESeq2/2025-11-24_J774A.1_E2_treatment_DESeq2_MM10_log2_counts_matrix.txt')

irf2.j <- list("Vehicle" = as.numeric(e2.cts["Irf2",1:3]),
               "E2"      = as.numeric(e2.cts["Irf2",4:7]))

pdf(paste0(Sys.Date(),"_Irf2_boxplot_JE2.pdf"), height = 5, width = 3)
boxplot(irf2.j, 
        col = c("deeppink","mediumorchid1"), 
        outline = F, las = 1,
        ylim = c(8.5,10),
        ylab = "DESeq2 VST−normalized log2 counts"
        )
beeswarm(irf2.j, add = T, pch = 16, cex = 1)
text(1.5,10, signif(e2.sig["Irf2",]$padj,2) )
dev.off()

#######################
sink(file = paste(Sys.Date(),"J774A.1_E2_RNAseq_analysis_session_Info.txt", sep =""))
sessionInfo()
sink()

