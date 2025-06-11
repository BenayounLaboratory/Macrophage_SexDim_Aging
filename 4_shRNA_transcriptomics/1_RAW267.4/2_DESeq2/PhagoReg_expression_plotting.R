setwd('/Volumes/BB_HQ_1/Immune_sex_dimorphism_Aging/Macrophages/Clean_Macrophage_Data/Intervention_Data/shRNA_analyses/RAW267_cells/DEseq2/')
options(stringsAsFactors = F)

library(beeswarm)


# 2025-04-29
# analyze Rps6ka3, Magt1 expression from RNAseq

##################################################################
# 1. Read RAW results
cts.meis1     <- read.table('2024-08-16_RAW_sh_Meis1_data_DESeq2_log2_counts_matrix.txt', sep = "\t", header = T)
deg.meis1     <- read.table('2024-08-16_RAW_sh_Meis1_data_DESeq2_shMeis1_DEseq2_results_COMBINED_all_genes_statistics.txt', sep = "\t", header = T)


shNon = 1:3
shLuc = 4:6
shTF1 = 7:9
shTF2 = 10:12

my.colors       <-  c(rep("orange",length(shNon)), rep("peru",length(shLuc)), rep("seagreen",length(shTF1)), rep("seagreen3",length(shTF2)))

Rps.exp <- list("shControl"   = as.numeric(cts.meis1["Rps6ka3", c(shNon, shLuc) ]),
                "shMeis1"     = as.numeric(cts.meis1["Rps6ka3",  c(shTF1, shTF2)]))

pdf(paste0(Sys.Date(),"_shIrf2_RAW267.4_Normalized_counts_boxplot_Rps6ka3.pdf"), height = 5, width = 3)
boxplot(Rps.exp, outline = F, las = 1, ylim = c(11,12) , main = "Rps6ka3 (RAW267.4)", ylab = "DEseq2 VST normalized log2 counts", cex.axis = 0.75)
beeswarm(Rps.exp, add = T, pch = 16, pwcol = my.colors)
text(1.5, 12, signif(deg.meis1[deg.meis1$Row.names == "Rps6ka3",]$p_combined,2))
dev.off()



#####
Mag.exp <- list("shControl"   = as.numeric(cts.meis1["Magt1", c(shNon, shLuc) ]),
                "shMeis1"     = as.numeric(cts.meis1["Magt1",  c(shTF1, shTF2)]))

pdf(paste0(Sys.Date(),"_shIrf2_RAW267.4_Normalized_counts_boxplot_Magt1.pdf"), height = 5, width = 3)
boxplot(Mag.exp, outline = F, las = 1, ylim = c(11,12) , main = "Magt1 (RAW267.4)", ylab = "DEseq2 VST normalized log2 counts", cex.axis = 0.75)
beeswarm(Mag.exp, add = T, pch = 16, pwcol = my.colors)
text(1.5, 12, signif(deg.meis1[deg.meis1$Row.names == "Magt1",]$p_combined,2))
dev.off()





# extract data frame of expression
my.gene.exp <- data.frame("YF" = as.numeric(cts.meis1[my.gene,1:5]  ) ,
                          "OF" = as.numeric(cts.meis1[my.gene,6:10]  ) ,
                          "YM" = as.numeric(cts.meis1[my.gene,11:15] ) ,
                          "OM" = as.numeric(cts.meis1[my.gene,16:20]) )

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




