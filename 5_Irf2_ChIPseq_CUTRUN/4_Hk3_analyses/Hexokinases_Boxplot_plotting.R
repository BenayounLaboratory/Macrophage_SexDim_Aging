setwd('/Volumes/BB_Travel_2/Mph_Aging/RNAseq/DEseq2_analysis/DESeq2/')
options(stringsAsFactors = FALSE)

# load libraries for analysis
library(DESeq2)
library(pheatmap)
library('bitops')
library('limma')
library(RColorBrewer)
library(fields)


# 2025-11-27
# plot Hexokinases expressions

###### read results
my.cts    <- read.table('/Volumes/BB_Travel_2/Mph_Aging/RNAseq/DEseq2_analysis/DESeq2/2022-03-30_Peritoneal_Macrophages_DESeq2_model_with_AGE_SEX_log2_counts_matrix.txt')
M.age.spe <- read.table('/Volumes/BB_Travel_2/Mph_Aging/RNAseq/DEseq2_analysis/DESeq2/2022-03-30_Peritoneal_Macrophages_AGING_Male.NOT.Female_FDR5_genes_statistics.txt')
F.age.spe <- read.table('/Volumes/BB_Travel_2/Mph_Aging/RNAseq/DEseq2_analysis/DESeq2/2022-03-30_Peritoneal_Macrophages_AGING_Female.NOT.Male_FDR5_genes_statistics.txt')
fem.mal.age.comp <- read.table('/Volumes/BB_Travel_2/Mph_Aging/RNAseq/DEseq2_analysis/DESeq2/2022-03-30_Peritoneal_Macrophages_AGING_Separated_Sex_Merged_Table_ALL_genes_statistics.txt')

# get ideal range
# quantile(my.cts["Hk3",])
# quantile(my.cts["Hk1",])
# quantile(my.cts["Hk2",])

get_sex_age_bplots <- function(my.gene, my.cts, my.sig.res, my.ylim = c(0,20), y.pval = 19) {
  # extract data frame of expression
  my.gene.exp <- data.frame("YF" = as.numeric(my.cts[my.gene,1:5]  ) ,
                            "OF" = as.numeric(my.cts[my.gene,6:10]  ) ,
                            "YM" = as.numeric(my.cts[my.gene,11:15] ) ,
                            "OM" = as.numeric(my.cts[my.gene,16:20]) )
  
  pdfname <- paste0(Sys.Date(),"_",my.gene,"_Gene_Expression_boxplot_DEseq2.pdf")
  pdf(pdfname, height = 5, width = 3)
  boxplot(my.gene.exp,
          col = c("deeppink", "deeppink4","deepskyblue", "deepskyblue4"),
          ylab = "DESeq2 VST-normalized log2 counts",
          ylim = my.ylim, outline = F, main = my.gene,
          las = 1)
  beeswarm::beeswarm(my.gene.exp, add = T, pch = 16, cex = 0.75)
  text(1.5, y.pval, signif(my.sig.res[my.gene,]$padj_F,2))
  text(3.5, y.pval, signif(my.sig.res[my.gene,]$padj_M,2))
  dev.off()
  
}


# ### hexokinase genes
get_sex_age_bplots("Hk3"    , my.cts, fem.mal.age.comp,  c(9, 12), 12)
get_sex_age_bplots("Hk1"    , my.cts, fem.mal.age.comp,  c(10, 13), 13)
get_sex_age_bplots("Hk2"    , my.cts, fem.mal.age.comp,  c(8, 11), 11)


#######################
sink(file = paste(Sys.Date(),"Macrophage_aging_sex_RNAseq_HEATMAPS_session_Info.txt", sep =""))
sessionInfo()
sink()

