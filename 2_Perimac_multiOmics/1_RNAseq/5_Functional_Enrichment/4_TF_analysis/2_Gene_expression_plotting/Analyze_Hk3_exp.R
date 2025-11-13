setwd('/Volumes/BB_HQ_1//Immune_sex_dimorphism_Aging/Macrophages/Clean_Macrophage_Data/RNAseq/Functional_enrichment/Top_TF_Validation/TF_analysis/noSVA/')
options(stringsAsFactors = F)

# 2024-09-26
# Plot Xist expression

#### load expression data
tissue.cts <- read.table('../../../../DEseq2_analysis/DEseq2/2022-03-30_Peritoneal_Macrophages_DESeq2_model_with_AGE_SEX_log2_counts_matrix.txt', sep = "\t", header = T)

# get Sex-specific changes
fem.mal.age.comp <- read.table('../../../../DEseq2_analysis/DEseq2/2022-03-30_Peritoneal_Macrophages_AGING_Separated_Sex_Merged_Table_ALL_genes_statistics.txt')
fem.spe          <- read.table('../../../../DEseq2_analysis/DEseq2/2022-03-30_Peritoneal_Macrophages_AGING_Female.NOT.Male_FDR5_genes_statistics.txt')
mal.spe          <- read.table('../../../../DEseq2_analysis/DEseq2/2022-03-30_Peritoneal_Macrophages_AGING_Male.NOT.Female_FDR5_genes_statistics.txt')


########################################################################################################
#### DEBUG
# my.gene    <- "Rps6ka3"
# my.cts     <- tissue.cts
# my.sig.res <- fem.mal.age.comp
# my.ylim    <- c(8,14)
# y.pval     <- 13

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
          ylim = my.ylim, outline = F, main = my.gene)
  beeswarm::beeswarm(my.gene.exp, add = T, pch = 16, cex = 0.75)
  text(1.5, y.pval, signif(my.sig.res[my.gene,]$padj_F,2))
  text(3.5, y.pval, signif(my.sig.res[my.gene,]$padj_M,2))
  dev.off()
  
}


# ### Female specific TF with changes in ATAC accessibility
get_sex_age_bplots("Hk3"    , tissue.cts, fem.mal.age.comp,  c(8, 12), 12)

#######################
sink(file = paste(Sys.Date(),"_session_Info.txt", sep =""))
sessionInfo()
sink()




