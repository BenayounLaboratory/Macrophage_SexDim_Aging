setwd('/Volumes/BB_Travel_2/New_Analyses/Hk_expression')
options(stringsAsFactors = FALSE)

# load libraries for analysis
library(DESeq2)
library(pheatmap)
library('bitops')
library('limma')
library(RColorBrewer)
library(fields)
library(beeswarm)


# 2025-11-27
# plot Hexokinases expressions

##################################################################
# 1. Read RAW/J7 results
raw.cts.irf2     <- read.table('/Volumes/BB_Travel_2/Mph_Aging/Intervention_Data/shRNA_analyses/RAW267_cells/DEseq2/2024-08-16_RAW_sh_Irf2_data_DESeq2_log2_counts_matrix.txt', sep = "\t", header = T)
raw.deg.irf2     <- read.table('/Volumes/BB_Travel_2/Mph_Aging/Intervention_Data/shRNA_analyses/RAW267_cells/DEseq2/2024-08-16_RAW_sh_Irf2_data_DESeq2_shIrf2_DEseq2_results_COMBINED_all_genes_statistics.txt', sep = "\t", header = T)

j77.cts.irf2     <- read.table('/Volumes/BB_Travel_2/Mph_Aging/Intervention_Data/shRNA_analyses/J774_cells/DESeq2/2025-01-30_J774_sh_Irf2_data_DESeq2_log2_counts_matrix.txt', sep = "\t", header = T)
j77.deg.irf2     <- read.table('/Volumes/BB_Travel_2/Mph_Aging/Intervention_Data/shRNA_analyses/J774_cells/DESeq2/2025-01-30_J774_sh_Irf2_data_DESeq2_shIrf2_DEseq2_results_COMBINED_all_genes_statistics.txt', sep = "\t", header = T)

shNon = 1:3
shLuc = 4:6
shTF1 = 7:9
shTF2 = 10:12

my.colors       <-  c(rep("orange",length(shNon)), rep("peru",length(shLuc)), rep("seagreen",length(shTF1)), rep("seagreen3",length(shTF2)))

Hk2.exp.raw <- list("RAW_shControl"      = as.numeric(raw.cts.irf2["Hk2", c(shNon, shLuc) ]),
                    "RAW_shIrf2"         = as.numeric(raw.cts.irf2["Hk2",  c(shTF1, shTF2)]),
                    "J774_shControl"     = as.numeric(j77.cts.irf2["Hk2", c(shNon, shLuc) ]),
                    "J774_shIrf2"        = as.numeric(j77.cts.irf2["Hk2",  c(shTF1, shTF2)]))

Hk1.exp.raw <- list("RAW_shControl"      = as.numeric(raw.cts.irf2["Hk1", c(shNon, shLuc) ]),
                    "RAW_shIrf2"         = as.numeric(raw.cts.irf2["Hk1",  c(shTF1, shTF2)]),
                    "J774_shControl"     = as.numeric(j77.cts.irf2["Hk1", c(shNon, shLuc) ]),
                    "J774_shIrf2"        = as.numeric(j77.cts.irf2["Hk1",  c(shTF1, shTF2)]))


pdf(paste0(Sys.Date(),"_shIrf2_RAW267.4_J774A.1_Normalized_counts_boxplot_Hk2.pdf"), height = 5, width = 5)
boxplot(Hk2.exp.raw, outline = F, las = 1, ylim = c(10,14) ,
        main = "Hk2", ylab = "DEseq2 VST normalized log2 counts", cex.axis = 0.75)
beeswarm(Hk2.exp.raw, add = T, pwcol = rep(my.colors,2), pwpch = c(rep(16, 12),rep(15, 12)))
text(1.5, 14, signif(raw.deg.irf2[raw.deg.irf2$Row.names == "Hk2",]$p_combined,2))
text(3.5, 14, signif(j77.deg.irf2[j77.deg.irf2$Row.names == "Hk2",]$p_combined,2))
abline(v = 2.5, col = "grey")
dev.off()

pdf(paste0(Sys.Date(),"_shIrf2_RAW267.4_J774A.1_Normalized_counts_boxplot_Hk1.pdf"), height = 5, width = 5)
boxplot(Hk1.exp.raw, outline = F, las = 1, ylim = c(10,14) ,
        main = "Hk1", ylab = "DEseq2 VST normalized log2 counts", cex.axis = 0.75)
beeswarm(Hk1.exp.raw, add = T, pwcol = rep(my.colors,2), pwpch = c(rep(16, 12),rep(15, 12)))
text(1.5, 14, signif(raw.deg.irf2[raw.deg.irf2$Row.names == "Hk1",]$p_combined,2))
text(3.5, 14, signif(j77.deg.irf2[j77.deg.irf2$Row.names == "Hk1",]$p_combined,2))
abline(v = 2.5, col = "grey")
dev.off()


#######################
sink(file = paste(Sys.Date(),"Macrophage_aging_sex_RNAseq_HEATMAPS_session_Info.txt", sep =""))
sessionInfo()
sink()

