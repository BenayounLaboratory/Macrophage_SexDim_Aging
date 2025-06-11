setwd('/Volumes/BB_HQ_1/Immune_sex_dimorphism_Aging/Macrophages/Clean_Macrophage_Data/Peritoneal_scRNAseq/v3_set/Augur')
options(stringsAsFactors = F)

library(Seurat)
library(Augur)
library(viridis)


######################################################################
# load Seurat analysis and SingleR annotation
load('../Seurat/2023-02-22_Seurat_JAX_10x_peritoneal_Singlets_SingleR_ConsistentID_SingleR.RData')
JAX.singlets.consistent
# An object of class Seurat 
# 49037 features across 14595 samples within 3 assays 
# Active assay: SCT (16740 features, 2000 variable features)
#  2 other assays present: RNA, HTO
#  2 dimensional reductions calculated: pca, umap


######################################################################
# Run AUGUR based on singleR annotation (consistent with both annotation libraries)
augur.perit <-  calculate_auc(as.matrix(JAX.singlets.consistent@assays$SCT@data),
                              JAX.singlets.consistent@meta.data, 
                              cell_type_col = "SingleR_ImmGen", 
                              label_col = "Condition",
                              n_threads = 1)
augur.perit
# $AUC
# # A tibble: 4 x 2
# cell_type     auc
# <chr>       <dbl>
# 1 Macrophages 0.587
# 2 B cells     0.533
# 3 Monocytes   0.529
# 4 T cells     0.511

save(augur.perit, file = paste0(Sys.Date(),"_JAX_Augur_peritoneal_cells_object_Consistent_labels.RData"))

# for some reason, barcode column trips Augur up
JAX.singlets.consistent@meta.data <- JAX.singlets.consistent@meta.data[,colnames(JAX.singlets.consistent@meta.data)  != "barcode"]

pdf(paste0(Sys.Date(),"_Augur_peritoneal_cells_UMAP_Consistent_labels.pdf"), width = 3, height = 3)
plot_umap(augur.perit,JAX.singlets.consistent, cell_type_col = "SingleR_RNAseq")
dev.off()

pdf(paste0(Sys.Date(),"_Augur_peritoneal_cells_UMAP_Red_Blue_Consistent_labels.pdf"), width = 3, height = 3)
plot_umap(augur.perit,JAX.singlets.consistent, cell_type_col = "SingleR_RNAseq", palette = colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50))
dev.off()

pdf(paste0(Sys.Date(),"_Augur_peritoneal_cells_Lollipop_Consistent_labels.pdf"), width = 3, height = 3)
plot_lollipop(augur.perit)
dev.off()
###################################################################################################################

###################################################################################################################
sink(file = paste(Sys.Date(),"_AUGUR_scRNAseq_peritoneal_analysis_session_Info.txt", sep =""))
sessionInfo()
sink()
