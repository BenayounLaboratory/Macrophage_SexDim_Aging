setwd('/Users/berenice/Desktop/LTE2/RNA-seq/Cohort6/DESeq2')
options(stringsAsFactors = F)

# load libraries for analysis
library(DESeq2)
library(Vennerable)
library(bitops)
library(beeswarm)

# 2025-11-26
# Plot Irf2 expression

####################################################################################
# Load E2 response
e2.sig  <- read.table('2025-11-26_Peritoneal_Macrophage_Aging_LTE2c6_DESeq2_Old_E2_all_genes_statistics.txt')
age.sig <- read.table('2025-11-26_Peritoneal_Macrophage_Aging_LTE2c6_DESeq2_Aging_CTL_all_genes_statistics.txt')
e2.cts  <- read.table('2025-11-26_Peritoneal_Macrophage_Aging_LTE2c6_DESeq2_log2_counts_matrix.txt')

irf2.j <- list("Young"   = as.numeric(e2.cts["Irf2",1:5]),
               "Old_Veh" = as.numeric(e2.cts["Irf2",6:9]),
               "Old_E2"  = as.numeric(e2.cts["Irf2",10:13]) )

pdf(paste0(Sys.Date(),"_Irf2_boxplot_LTE2.pdf"), height = 5, width = 4)
boxplot(irf2.j, 
        col = c("deeppink","deeppink4","magenta"), 
        outline = F, las = 1,
        ylim = c(10,14),
        ylab = "DESeq2 VST−normalized log2 counts"
        )
beeswarm(irf2.j, add = T, pch = 16, cex = 1)
text(1.5,14, signif(age.sig["Irf2",]$padj,2) )
text(2.5,14, signif(e2.sig["Irf2",]$padj,2) )
dev.off()

#######################
sink(file = paste(Sys.Date(),"J774A.1_E2_RNAseq_analysis_session_Info.txt", sep =""))
sessionInfo()
sink()

