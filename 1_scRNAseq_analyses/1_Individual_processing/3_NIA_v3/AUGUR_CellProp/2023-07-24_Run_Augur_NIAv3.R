setwd('/Volumes/BB_HQ_1/Immune_sex_dimorphism_Aging/Macrophages/Clean_Macrophage_Data/Peritoneal_scRNAseq/v3_NIA/AUGUR')
options(stringsAsFactors = F)

library(Seurat)
library(Augur)
library(viridis)


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
# Run AUGUR based on singleR annotation (consistent with both annotation libraries)
augur.perit <-  calculate_auc(as.matrix(NIA.singlets.consistent@assays$SCT@data),
                              NIA.singlets.consistent@meta.data, 
                              cell_type_col = "SingleR_ImmGen", 
                              label_col = "Condition",
                              n_threads = 1)
augur.perit
# $AUC
# # A tibble: 3 x 2
# cell_type     auc
# <chr>       <dbl>
# 1 Macrophages 0.567
# 2 B cells     0.545
# 3 T cells     0.517

save(augur.perit, file = paste0(Sys.Date(),"_NIAv3_Augur_peritoneal_cells_object_Consistent_labels.RData"))

# for some reason, barcode column trips Augur up
NIA.singlets.consistent@meta.data <- NIA.singlets.consistent@meta.data[,colnames(NIA.singlets.consistent@meta.data)  != "barcode"]

pdf(paste0(Sys.Date(),"_Augur_peritoneal_cells_NIAv3_UMAP_Consistent_labels.pdf"), width = 3, height = 3)
plot_umap(augur.perit,NIA.singlets.consistent, cell_type_col = "SingleR_RNAseq")
dev.off()

pdf(paste0(Sys.Date(),"_Augur_peritoneal_cells_NIAv3_UMAP_Red_Blue_Consistent_labels.pdf"), width = 3, height = 3)
plot_umap(augur.perit,NIA.singlets.consistent, cell_type_col = "SingleR_RNAseq", palette = colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50))
dev.off()

pdf(paste0(Sys.Date(),"_Augur_peritoneal_cells_NIAv3_Lollipop_Consistent_labels.pdf"), width = 3, height = 3)
plot_lollipop(augur.perit)
dev.off()
###################################################################################################################

###################################################################################################################
sink(file = paste(Sys.Date(),"_AUGUR_scRNAseq_NIAv3_peritoneal_analysis_session_Info.txt", sep =""))
sessionInfo()
sink()
