setwd('/Volumes/BB_HQ_1/Immune_sex_dimorphism_Aging/Macrophages/Clean_Macrophage_Data/Peritoneal_scRNAseq/Combined_Analysis/Macrophage_Analysis')
options(stringsAsFactors=FALSE)

library(Seurat)
library(beeswarm)
library(ggplot2)

########################################################
# 0. Load macrophages from all experiments 

# Load all 3 cohort-wise SingleR Seurat Objects, only consistent macrophages

#%%%%%%%%%%%%%%% NIA v2 %%%%%%%%%%%%%%%
load('../../v2_set/Seurat_Doublet_Identification_clustering/2023-02-22_Seurat_10x_peritoneal_Singlets_SingleR_Consistent_Macrophages_NIA.RData')
my.perit.singlets.Mph.v3
# An object of class Seurat 
# 48214 features across 4935 samples within 2 assays 
# Active assay: SCT (15929 features, 2000 variable features)
# 1 other assay present: RNA
# 2 dimensional reductions calculated: pca, umap

#%%%%%%%%%%%%%%% JAX v3 %%%%%%%%%%%%%%%
load('../../v3_set/Seurat/2023-02-22_Seurat_10x_peritoneal_Singlets_SingleR_Consistent_Macrophages.RData')
JAX.singlets.Mph.v3
# An object of class Seurat 
# 49037 features across 5475 samples within 3 assays 
# Active assay: SCT (16740 features, 2000 variable features)
# 2 other assays present: RNA, HTO
# 2 dimensional reductions calculated: pca, umap

# https://github.com/satijalab/seurat/issues/4633
# remove HTO for SCE to not bug out
JAX.Mph.diet <- DietSeurat(JAX.singlets.Mph.v3, assays = c("SCT","RNA"), dimreducs = c("pca","umap"))
JAX.Mph.diet
# An object of class Seurat 
# 49025 features across 5475 samples within 2 assays 
# Active assay: SCT (16740 features, 2000 variable features)
# 1 other assay present: RNA
# 2 dimensional reductions calculated: pca, umap

# subset metadata
JAX.Mph.diet@meta.data <- JAX.Mph.diet@meta.data[,c("nCount_RNA"       ,
                                                    "nFeature_RNA"     ,
                                                    "percent.mito"     ,
                                                    "S.Score"          ,
                                                    "G2M.Score"        ,
                                                    "Phase"            ,
                                                    "nCount_SCT"       ,
                                                    "nFeature_SCT"     ,
                                                    "Condition"        ,
                                                    "Sample_ID"        ,
                                                    "SingleR_ImmGen"   ,
                                                    "SingleR_RNAseq"   )]
JAX.Mph.diet@meta.data$source <- "JAX_v3"


#%%%%%%%%%%%%%%% NIA v3 %%%%%%%%%%%%%%%
load('../../v3_NIA/Seurat/2023-07-21_Seurat_10x_peritoneal_Singlets_SingleR_Consistent_Macrophages.RData')
NIA.singlets.Mph.v3
# An object of class Seurat 
# 48285 features across 5573 samples within 3 assays 
# Active assay: SCT (15990 features, 2000 variable features)
# 2 other assays present: RNA, HTO
# 2 dimensional reductions calculated: pca, umap

# https://github.com/satijalab/seurat/issues/4633
# remove HTO for SCE to not bug out
NIA.mph.diet <- DietSeurat(NIA.singlets.Mph.v3, assays = c("SCT","RNA"), dimreducs = c("pca","umap"))
NIA.mph.diet
# An object of class Seurat 
# 48275 features across 5573 samples within 2 assays 
# Active assay: SCT (15990 features, 2000 variable features)
# 1 other assay present: RNA
# 2 dimensional reductions calculated: pca, umap

# subset metadata
NIA.mph.diet@meta.data <- NIA.mph.diet@meta.data[,c("nCount_RNA"       ,
                                                    "nFeature_RNA"     ,
                                                    "percent.mito"     ,
                                                    "S.Score"          ,
                                                    "G2M.Score"        ,
                                                    "Phase"            ,
                                                    "nCount_SCT"       ,
                                                    "nFeature_SCT"     ,
                                                    "Condition"        ,
                                                    "Sample_ID"        ,
                                                    "SingleR_ImmGen"   ,
                                                    "SingleR_RNAseq"   )]
NIA.mph.diet@meta.data$source <- "NIA_v3"

################
##### Clean factor
JAX.Mph.diet$Condition <- factor(JAX.Mph.diet$Condition, levels = c("YF", "OF", "YM", "OM"))
NIA.mph.diet$Condition <- factor(NIA.mph.diet$Condition, levels = c("YF", "OF", "YM", "OM"))
my.perit.singlets.Mph.v3$Condition <- factor(my.perit.singlets.Mph.v3$Condition, levels = c("YF", "OF", "YM", "OM"))

library(ggpubr)

##### Irf2
p.irf.v2  <- VlnPlot(my.perit.singlets.Mph.v3, features = "Irf2", group.by = "Condition", col = c("deeppink","deeppink4","deepskyblue","deepskyblue4"))
p.irf.v3N <- VlnPlot(NIA.mph.diet            , features = "Irf2", group.by = "Condition", col = c("deeppink","deeppink4","deepskyblue","deepskyblue4"))
p.irf.v3J <- VlnPlot(JAX.Mph.diet            , features = "Irf2", group.by = "Condition", col = c("deeppink","deeppink4","deepskyblue","deepskyblue4"))
p.irf <- ggarrange(p.irf.v2, p.irf.v3N, p.irf.v3J, ncol = 3, nrow = 1, labels = c("NIA v2", "NIA v3", "JAX v3"),  common.legend = T, legend = "right")
p.irf

pdf(paste0(Sys.Date(),"_Irf2_Exp_violins_3_cohorts.pdf"), height = 3.5, width = 7)
plot(p.irf)
dev.off()

##### Esr1
p.esr.v2  <- VlnPlot(my.perit.singlets.Mph.v3, features = "Esr1", group.by = "Condition", col = c("deeppink","deeppink4","deepskyblue","deepskyblue4"))
p.esr.v3N <- VlnPlot(NIA.mph.diet            , features = "Esr1", group.by = "Condition", col = c("deeppink","deeppink4","deepskyblue","deepskyblue4"))
p.esr.v3J <- VlnPlot(JAX.Mph.diet            , features = "Esr1", group.by = "Condition", col = c("deeppink","deeppink4","deepskyblue","deepskyblue4"))
p.esr <- ggarrange(p.esr.v2, p.esr.v3N, p.esr.v3J, ncol = 3, nrow = 1, labels = c("NIA v2", "NIA v3", "JAX v3"),  common.legend = T, legend = "right")
p.esr

pdf(paste0(Sys.Date(),"_Esr1_Exp_violins_3_cohorts.pdf"), height = 3.5, width = 7)
plot(p.esr)
dev.off()

##### Xist
p.Xist.v2  <- VlnPlot(my.perit.singlets.Mph.v3, features = "Xist", group.by = "Condition", col = c("deeppink","deeppink4","deepskyblue","deepskyblue4"))
p.Xist.v3N <- VlnPlot(NIA.mph.diet            , features = "Xist", group.by = "Condition", col = c("deeppink","deeppink4","deepskyblue","deepskyblue4"))
p.Xist.v3J <- VlnPlot(JAX.Mph.diet            , features = "Xist", group.by = "Condition", col = c("deeppink","deeppink4","deepskyblue","deepskyblue4"))
p.Xist     <- ggarrange(p.Xist.v2, p.Xist.v3N, p.Xist.v3J, ncol = 3, nrow = 1, labels = c("NIA v2", "NIA v3", "JAX v3"),  common.legend = T, legend = "right")
p.Xist

pdf(paste0(Sys.Date(),"_Xist_Exp_violins_3_cohorts.pdf"), height = 3.5, width = 7)
plot(p.Xist)
dev.off()


################################################
# get tests on violins

Idents(my.perit.singlets.Mph.v3) <- "Condition"
Idents(NIA.mph.diet            ) <- "Condition"
Idents(JAX.Mph.diet            ) <- "Condition"


#
get_test_scRNA <- function (my.feature, my.obj.list = list(my.perit.singlets.Mph.v3, NIA.mph.diet, JAX.Mph.diet)) {
  
  # calculate wilcoxon through Seurat for feature of interest
  F.v2  <- FindMarkers(my.obj.list[[1]], ident.1 = "YF", ident.2 = "OF", features = my.feature, logfc.threshold = 0, min.pct = 0)
  F.v3N <- FindMarkers(my.obj.list[[2]], ident.1 = "YF", ident.2 = "OF", features = my.feature, logfc.threshold = 0, min.pct = 0)
  F.v3J <- FindMarkers(my.obj.list[[3]], ident.1 = "YF", ident.2 = "OF", features = my.feature, logfc.threshold = 0, min.pct = 0)
  
  M.v2  <- FindMarkers(my.obj.list[[1]], ident.1 = "YM", ident.2 = "OM", features = my.feature, logfc.threshold = 0, min.pct = 0)
  M.v3N <- FindMarkers(my.obj.list[[2]], ident.1 = "YM", ident.2 = "OM", features = my.feature, logfc.threshold = 0, min.pct = 0)
  M.v3J <- FindMarkers(my.obj.list[[3]], ident.1 = "YM", ident.2 = "OM", features = my.feature, logfc.threshold = 0, min.pct = 0)
  
  # collate for return
  my.wilcox.pvals <- c(F.v2 $p_val, M.v2 $p_val, F.v3N$p_val, M.v3N$p_val, F.v3J$p_val, M.v3J$p_val)
  names(my.wilcox.pvals) <- c("F_v2"  , "M_v2", "F_v3N" , "M_v3N", "F_v3J" , "M_v3J")
  
  my.wilcox.pvals
}

Irf2.pvals <- get_test_scRNA("Irf2")
Esr1.pvals <- get_test_scRNA("Esr1")
Xist.pvals <- get_test_scRNA("Xist")

my.res <- rbind(Irf2.pvals,
                Esr1.pvals,
                Xist.pvals)
#                   F_v2       M_v2        F_v3N        M_v3N        F_v3J        M_v3J
# Irf2.pvals 0.401205884 0.02016652 1.698433e-04 3.962175e-04 4.265559e-20 7.334037e-08
# Esr1.pvals 0.000340508 0.41383242 3.193705e-01 1.265224e-01 4.309920e-01 1.575560e-03
# Xist.pvals 0.002258074 0.03778703 9.240861e-10 1.642168e-58 2.764026e-19 2.247102e-48

write.table(my.res , file = "Esr1_Irf2_Xist_Wilcox_scRNA_perCohort.txt", col.names = T, row.names = T, quote = F, sep = "\t")

############################################################################################################

#######################
sink(file = paste(Sys.Date(),"_scRNAseq_macrophage_v3Cohorts_COMBINED_TopGene_Analysis.txt", sep =""))
sessionInfo()
sink()

