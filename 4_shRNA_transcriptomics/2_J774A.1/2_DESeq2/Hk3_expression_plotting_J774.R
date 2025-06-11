setwd('/Volumes/BB_HQ_1/Immune_sex_dimorphism_Aging/Macrophages/Clean_Macrophage_Data/Intervention_Data/shRNA_analyses/J774_cells/DESeq2')
options(stringsAsFactors = F)

library(beeswarm)

# 2025-01-30
# analyze Hk3 expression from RNAseq

##################################################################
# 1. Read RAW results
cts.irf2     <- read.table('2025-01-30_J774_sh_Irf2_data_DESeq2_log2_counts_matrix.txt', sep = "\t", header = T)
deg.irf2     <- read.table('2025-01-30_J774_sh_Irf2_data_DESeq2_shIrf2_DEseq2_results_COMBINED_FDR1e-4_genes_statistics.txt', sep = "\t", header = T)


shNon = 1:3
shLuc = 4:6
shTF1 = 7:9
shTF2 = 10:12

my.colors       <-  c(rep("orange",length(shNon)), rep("peru",length(shLuc)), rep("seagreen",length(shTF1)), rep("seagreen3",length(shTF2)))

Hk.exp <- list("shControl"   = as.numeric(cts.irf2["Hk3", c(shNon, shLuc) ]),
               "shIrf2"        = as.numeric(cts.irf2["Hk3",  c(shTF1, shTF2)]))

pdf(paste0(Sys.Date(),"_shIrf2_J774A.1_Normalized_counts_boxplot_Hk3.pdf"), height = 5, width = 3)
boxplot(Hk.exp, outline = F, las = 1, ylim = c(12,15) , main = "Hk3 (J774A.1)", ylab = "DEseq2 VST normalized log2 counts", cex.axis = 0.75)
beeswarm(Hk.exp, add = T, pch = 15, pwcol = my.colors)
text(1.5, 15, signif(deg.irf2[deg.irf2$Row.names == "Hk3",]$p_combined,2))
dev.off()




# extract data frame of expression
my.gene.exp <- data.frame("YF" = as.numeric(cts.irf2[my.gene,1:5]  ) ,
                          "OF" = as.numeric(cts.irf2[my.gene,6:10]  ) ,
                          "YM" = as.numeric(cts.irf2[my.gene,11:15] ) ,
                          "OM" = as.numeric(cts.irf2[my.gene,16:20]) )

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




#######################
sink(file = paste(Sys.Date(),"shTF_RAWcells_RNAseq_session_Info.txt", sep =""))
sessionInfo()
sink()




