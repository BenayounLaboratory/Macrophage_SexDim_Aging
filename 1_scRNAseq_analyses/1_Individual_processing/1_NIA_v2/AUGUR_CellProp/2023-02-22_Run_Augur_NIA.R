setwd('/Volumes/BB_HQ_1/Immune_sex_dimorphism_Aging/Macrophages/Clean_Macrophage_Data/Peritoneal_scRNAseq/v2_set/Augur')
options(stringsAsFactors = F)


library(Seurat)
library(Augur)
library(viridis)


######################################################################
# load Seurat analysis and SingleR annotation
load('../Seurat_Doublet_Identification_clustering/2023-02-22_Seurat_10x_peritoneal_Singlets_SingleR_ConsistentID_SingleR.RData')
my.perit.singlets.consistent
# An object of class Seurat 
# 48214 features across 18805 samples within 2 assays 
# Active assay: SCT (15929 features, 2000 variable features)
#  1 other assay present: RNA
#  2 dimensional reductions calculated: pca, umap


######################################################################
# Run AUGUR based on singleR annotation (consistent with both annotation libraries)
augur.perit <-  calculate_auc(as.matrix(my.perit.singlets.consistent@assays$SCT@data),
                              my.perit.singlets.consistent@meta.data, 
                              cell_type_col = "SingleR_ImmGen", 
                              label_col = "Condition",
                              n_threads = 1)
augur.perit
# $AUC
# # A tibble: 4 x 2
# cell_type     auc
# <chr>       <dbl>
# 1 Macrophages 0.592
# 2 B cells     0.538
# 3 Monocytes   0.531
# 4 T cells     0.515

save(augur.perit, file = paste0(Sys.Date(),"_Augur_peritoneal_cells_v2NIA_object_Consistent_labels.RData"))

pdf(paste0(Sys.Date(),"_Augur_peritoneal_cells_UMAP_v2NIA_Consistent_labels.pdf"), width = 3, height = 3)
plot_umap(augur.perit,my.perit.singlets.consistent, cell_type_col = "SingleR_RNAseq")
dev.off()

pdf(paste0(Sys.Date(),"_Augur_peritoneal_cells_UMAP_Red_Blue_v2NIA_Consistent_labels.pdf"), width = 3, height = 3)
plot_umap(augur.perit,my.perit.singlets.consistent, cell_type_col = "SingleR_RNAseq", palette = colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50))
dev.off()

pdf(paste0(Sys.Date(),"_Augur_peritoneal_cells_Lollipop_v2NIA_Consistent_labels.pdf"), width = 3, height = 3)
plot_lollipop(augur.perit)
dev.off()

#######################
sink(file = paste(Sys.Date(),"_AUGUR_scRNAseq_peritoneal_v2NIA_analysis_session_Info.txt", sep =""))
sessionInfo()
sink()
