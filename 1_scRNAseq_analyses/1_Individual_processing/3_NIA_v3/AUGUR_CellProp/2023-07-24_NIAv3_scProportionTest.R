setwd('/Volumes/BB_HQ_1/Immune_sex_dimorphism_Aging/Macrophages/Clean_Macrophage_Data/Peritoneal_scRNAseq/v3_NIA/scProportionTest')
options(stringsAsFactors = F)

library(Seurat)
library(scProportionTest)
library(ggplot2)


######################################################################
# load Seurat analysis and SingleR annotation
load('../Seurat/2023-07-21_Seurat_NIAv3_10x_peritoneal_Singlets_SingleR_ConsistentID_SingleR.RData')
NIA.singlets.consistent
# An object of class Seurat 
# 48285 features across 10758 samples within 3 assays 
# Active assay: SCT (15990 features, 2000 variable features)
# 2 other assays present: RNA, HTO
# 2 dimensional reductions calculated: pca, umap

######################################################################
# create prop test object
perit.NIAv3.prop_test <- sc_utils(NIA.singlets.consistent)


# Once the object is created, the permutation testing and bootstrapping can be run.
perit.NIAv3.prop_test.F <- permutation_test(perit.NIAv3.prop_test, 
                                            cluster_identity = "SingleR_ImmGen",
                                            sample_1 = "YF", 
                                            sample_2 = "OF",
                                            sample_identity = "Condition")

perit.NIAv3.prop_test.M <- permutation_test(perit.NIAv3.prop_test, 
                                            cluster_identity = "SingleR_ImmGen",
                                            sample_1 = "YM", 
                                            sample_2 = "OM",
                                            sample_identity = "Condition")

## Modify function
permutation_plot_mod <- function (sc_utils_obj, FDR_threshold = 0.05, cols_vals = c("deeppink", "deeppink4","grey"), my.title) {
  plot_data <- sc_utils_obj@results$permutation
  
  # fix order alphabetically
  plot_data$clusters <- factor(plot_data$clusters, levels = rev(plot_data$clusters))
  
  # get significance and color scale
  plot_data$significance <- "n.s."
  plot_data$significance[plot_data$FDR < FDR_threshold & plot_data$boot_mean_log2FD > 0] <- paste("FDR <", FDR_threshold, "(Increased)")
  plot_data$significance[plot_data$FDR < FDR_threshold & plot_data$boot_mean_log2FD < 0] <- paste("FDR <", FDR_threshold, "(Decreased)")
  plot_data$significance <- factor(plot_data$significance, levels = c(paste("FDR <", FDR_threshold, "(Decreased)"), paste("FDR <", FDR_threshold, "(Increased)"),  "n.s."))
  
  # fix range for symmetry
  max_range <- max(abs(plot_data[,7:8]))
  
  # plot
  p <- ggplot(plot_data, aes(x = clusters, y = obs_log2FD)) +   theme_bw() +
         geom_pointrange(aes(ymin = boot_CI_2.5, ymax = boot_CI_97.5, color = significance)) +
         geom_hline(yintercept = 0) + scale_color_manual(values = cols_vals) + ylim(c(-max_range,max_range)) + 
         coord_flip() + ggtitle(my.title)
  # p
  return(p)
}

# A point-range plot of the results can then be created.
perm.f <- permutation_plot_mod(perit.NIAv3.prop_test.F, cols_vals = c("deeppink"   , "deeppink4","grey")   , my.title = "Female Aging (NIAv3)")
perm.m <- permutation_plot_mod(perit.NIAv3.prop_test.M, cols_vals = c("deepskyblue", "deepskyblue4","grey"), my.title = "Male Aging (NIAv3)")

pdf(paste(Sys.Date(),"scProportionTest_SingleR_CONSISTENT_NIAv3Cohort.pdf",sep = "_"), height = 4, width = 10)
perm.f + perm.m
dev.off()
############################################################################################################

#######################
sink(file = paste(Sys.Date(),"_scProportionTest_NIAv3Cohort_analysis_session_Info.txt", sep =""))
sessionInfo()
sink()
