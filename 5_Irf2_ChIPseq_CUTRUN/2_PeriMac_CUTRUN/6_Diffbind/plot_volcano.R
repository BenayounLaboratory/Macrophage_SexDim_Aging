setwd('/Users/berenice/Desktop/Irf2_CUTRUN/DiffBind/MSPC_4samples')
options(stringsAsFactors = F)

library("scales")

# 2025-12-02
# use diffbind to analyze Irf cut and run
# use stringent peak definition to improve statistical power (peak in all samples)
################################################################################

diffbind.res <- read.csv('2025-11-28_PeriMac_Irf2_stringent_Aging_DiffBind_DESeq2_ALL_Peaks.txt', header = T, sep = "\t")

my.sig <- diffbind.res$FDR<0.05
sum(my.sig) # 541

pdf(paste0(Sys.Date(),"_Diffbind_DESeq2_Irf2_CUTRUN_Volcano.pdf"), height = 4.25, width = 3)
plot(diffbind.res$Fold, -log10(diffbind.res$FDR), 
     pch  = 16, cex = 0.5, 
     col  = alpha("grey50", alpha = 0.3),
     xlab = "log2 FC",
     ylab = "-log10(FDR)",
     xlim = c(-2,2),
     ylim = c(0,2) ,
     las  = 1)
points(diffbind.res$Fold[my.sig], -log10(diffbind.res$FDR)[my.sig], 
       pch = 1, cex = 0.5, col = alpha("deeppink", alpha = 0.3))
abline(h = -log10(0.05), col = "red", lty ="dashed")
abline(v = 0, col = "grey", lty ="dashed")
dev.off()

pdf(paste0(Sys.Date(),"_Diffbind_DESeq2_Irf2_CUTRUN_Volcano_with_Hk3.pdf"), height = 4.25, width = 3)
plot(diffbind.res$Fold, -log10(diffbind.res$FDR), 
     pch  = 16, cex = 0.5, 
     col  = alpha("grey50", alpha = 0.3),
     xlab = "log2 FC",
     ylab = "-log10(FDR)",
     xlim = c(-2,2),
     ylim = c(0,2) ,
     las  = 1)
points(diffbind.res$Fold[my.sig], -log10(diffbind.res$FDR)[my.sig], 
       pch = 1, cex = 0.5, col = alpha("deeppink", alpha = 0.3))
abline(h = -log10(0.05), col = "red", lty ="dashed")
abline(v = 0, col = "grey", lty ="dashed")
points(diffbind.res$Fold[diffbind.res$Gene.Name == "Hk3"],
       -log10(diffbind.res$FDR)[diffbind.res$Gene.Name == "Hk3"],
       pch = 1, cex = 0.5, col = "black")
text(diffbind.res$Fold[diffbind.res$Gene.Name == "Hk3"],
     -log10(diffbind.res$FDR)[diffbind.res$Gene.Name == "Hk3"],
     "Hk3", pos = 4)
dev.off()

################################################################################

#######################
sink(file = paste(Sys.Date(),"Diffbind_Irf2_CUTandRUN_session_Info.txt", sep =""))
sessionInfo()
sink()

