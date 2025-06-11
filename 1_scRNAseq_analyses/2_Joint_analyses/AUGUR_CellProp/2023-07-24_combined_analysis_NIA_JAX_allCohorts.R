setwd('/Volumes/BB_HQ_1/Immune_sex_dimorphism_Aging/Macrophages/Clean_Macrophage_Data/Peritoneal_scRNAseq/Combined_Analysis/Cell_Proportions_AUGUR/')
options(stringsAsFactors = F)

# load libraries for analysis
library('Seurat')
library(bitops)
library('ggplot2')
library('dplyr')
library(Augur)
library(viridis)
library(beeswarm)
library(scProportionTest)
library(ggplot2)



################################################################################################################
### 1. load annotated datasets from NIA (v2) and JAX (v3)

# Load all 3 cohort-wise SingleR Seurat Objects, only consistent cells
load('../../v2_set/Seurat_Doublet_Identification_clustering/2023-02-22_Seurat_10x_peritoneal_Singlets_SingleR_ConsistentID_SingleR.RData')
load('../../v3_set/Seurat/2023-02-22_Seurat_JAX_10x_peritoneal_Singlets_SingleR_ConsistentID_SingleR.RData')
load('../../v3_NIA/Seurat/2023-07-21_Seurat_NIAv3_10x_peritoneal_Singlets_SingleR_ConsistentID_SingleR.RData')

###########################
### NIA data (v2)
my.perit.singlets.consistent
# An object of class Seurat 
# 48214 features across 18805 samples within 2 assays 
# Active assay: SCT (15929 features, 2000 variable features)
#  1 other assay present: RNA
#  2 dimensional reductions calculated: pca, umap

# add provenance
my.perit.singlets.consistent@meta.data$source <- "NIA_v2"


###########################
### JAX data
JAX.singlets.consistent
# An object of class Seurat 
# 49037 features across 14595 samples within 3 assays 
# Active assay: SCT (16740 features, 2000 variable features)
#  2 other assays present: RNA, HTO
#  2 dimensional reductions calculated: pca, umap

# https://github.com/satijalab/seurat/issues/4633
# remove HTO for SCE to not bug out
JAX.singlets.diet <- DietSeurat(JAX.singlets.consistent, assays = c("SCT","RNA"), dimreducs = c("pca","umap"))
JAX.singlets.diet
# An object of class Seurat 
# 49025 features across 14595 samples within 2 assays 
# Active assay: SCT (16740 features, 2000 variable features)
#  1 other assay present: RNA
#  2 dimensional reductions calculated: pca, umap

# subset metadata
JAX.singlets.diet@meta.data <- JAX.singlets.diet@meta.data[,c("nCount_RNA"       ,
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
JAX.singlets.diet@meta.data$source <- "JAX_v3"


###########################
### NIA data  (v3)
NIA.singlets.consistent
# An object of class Seurat 
# 48285 features across 10758 samples within 3 assays 
# Active assay: SCT (15990 features, 2000 variable features)
#  2 other assays present: RNA, HTO
#  2 dimensional reductions calculated: pca, umap

# https://github.com/satijalab/seurat/issues/4633
# remove HTO for SCE to not bug out
NIA.singlets.diet <- DietSeurat(NIA.singlets.consistent, assays = c("SCT","RNA"), dimreducs = c("pca","umap"))
NIA.singlets.diet
# An object of class Seurat 
# 48275 features across 10758 samples within 2 assays 
# Active assay: SCT (15990 features, 2000 variable features)
# 1 other assay present: RNA
# 2 dimensional reductions calculated: pca, umap

# subset metadata
NIA.singlets.diet@meta.data <- NIA.singlets.diet@meta.data[,c("nCount_RNA"       ,
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
NIA.singlets.diet@meta.data$source <- "NIA_v3"



# make object list for processing
perit.list <- list("NIA_v2" = my.perit.singlets.consistent,
                   "JAX_v3" = JAX.singlets.diet,
                   "NIA_v3" = NIA.singlets.diet)


################################################################################################################
### 2. Integration with RPCA

# normalize and identify variable features for each dataset independently
perit.list <- lapply(X = perit.list,
                     FUN = function(x) {
                       x <- NormalizeData(x)
                       x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = perit.list, nfeatures = 3000)

# run PCA on each dataset using these features
perit.list <- lapply(X = perit.list,
                     FUN = function(x) {
                       x <- ScaleData(x, features = features, verbose = FALSE)
                       x <- RunPCA(x, features = features, verbose = FALSE)})

# Perform integration
perit.anchors <- FindIntegrationAnchors(object.list = perit.list, anchor.features = features, reduction = "rpca", k.anchor = 10)

# this command creates an 'integrated' data assay
perit.combined <- IntegrateData(anchorset = perit.anchors)
perit.combined
# An object of class Seurat 
# 52639 features across 44158 samples within 3 assays 
# Active assay: integrated (3000 features, 3000 variable features)
# 2 other assays present: RNA, SCT

# Now we can run a single integrated analysis on all cells!
# specify that we will perform downstream analysis on the corrected data
# note that the original unmodified data still resides in the 'RNA' assay
DefaultAssay(perit.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
perit.combined <- ScaleData(perit.combined, verbose = FALSE)
perit.combined <- RunPCA(perit.combined, npcs = 50, verbose = FALSE)

# Determine the ‘dimensionality’ of the dataset
pdf(paste0(Sys.Date(), "_Peritoneal_Lavage_COMBINED_elbowplot.pdf"))
ElbowPlot(perit.combined, ndims = 50)
dev.off()


################################################################################
# https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html
# To give us an idea of the number of PCs needed to be included:
# We can calculate where the principal components start to elbow by taking the larger value of:
#    - The point where the principal components only contribute 5% of standard deviation
#    - The principal components cumulatively contribute 90% of the standard deviation.
#    - The point where the percent change in variation between the consecutive PCs is less than 0.1%.

# Determine percent of variation associated with each PC
pct <- perit.combined[["pca"]]@stdev / sum(perit.combined[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
co1 # 41

# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

# last point where change of % of variation is more than 0.1%.
co2 # 15

# Minimum of the two calculation
pcs <- min(co1, co2)
pcs # 15

# Based on these metrics, first 16 PCs to generate the clusters.
# We can plot the elbow plot again and overlay the information determined using our metrics:

# Create a dataframe with values
plot_df <- data.frame(pct  = pct,
                      cumu = cumu,
                      rank = 1:length(pct))

# Elbow plot to visualize
pdf(paste0(Sys.Date(), "_Peritoneal_Lavage_COMBINED_elbowplot_threshold_analysis.pdf"), height = 5, width= 6)
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) +
  geom_text() +
  geom_vline(xintercept = 90, color = "grey") +
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()
dev.off()
###############################################################################


# Run UMAP and clustering based on dimensionality
perit.combined <- RunUMAP(perit.combined, reduction = "pca", dims = 1:pcs)
perit.combined <- FindNeighbors(perit.combined, reduction = "pca", dims = 1:pcs)
perit.combined <- FindClusters(perit.combined, resolution = 0.3)

# Visualization
p1 <- DimPlot(perit.combined, reduction = "umap", group.by = "source", shuffle = T)
p2 <- DimPlot(perit.combined, reduction = "umap", group.by = "SingleR_ImmGen", label = TRUE, repel = TRUE)

pdf(paste0(Sys.Date(),"PostIntegration_COMBINED_Peritoneal_aging_UMAP_SOURCE_and_SingleR_ImmGen.pdf"), width = 12, height = 5)
p1 + p2
dev.off()


table(perit.combined$Sample_ID)
#   JAX_OF1   JAX_OF2   JAX_OF3   JAX_OM1   JAX_OM2   JAX_OM3   JAX_YF1   JAX_YF2   JAX_YF3   JAX_YM1   JAX_YM2   JAX_YM3   NIA_OF1   NIA_OF2   NIA_OM1   NIA_OM2 
#      1221      1081      1319      1154      1157      1233      1042      1382      1303      1026      1363      1314       619      1293      1009      1317 
#   NIA_YF1   NIA_YF2   NIA_YF3   NIA_YM1   NIA_YM2   NIA_YM3 OFCohort1 OFCohort2 OMCohort1 OMCohort2 YFCohort1 YFCohort2 YMCohort1 YMCohort2 
#       670      1620       980      1601       979       670      3342      2119      3538      1663      3879      1777      1046      1441 



################################################################################################################
### 3. Integrated analysis on all cells for cell proportion (library-wise)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(perit.combined) <- "integrated"

# keep only cells with consistent labeling
pdf(paste0(Sys.Date(),"SingleR_consistent_COMBINED_Peritoneal_aging_UMAP_SingleR_ImmGen_CLEAN.pdf"), width = 7, height = 5)
DimPlot(perit.combined, reduction = "umap", group.by = "SingleR_ImmGen")
dev.off()

pdf(paste0(Sys.Date(),"SingleR_consistent_COMBINED_Peritoneal_aging_UMAP_SingleR_ImmGen_split_by_sample.pdf"), width = 20, height = 25)
DimPlot(perit.combined, reduction = "umap", group.by = "SingleR_ImmGen", split.by = "Sample_ID", ncol = 4)
dev.off()


# get cell type proportions
my.freq.table <- prop.table(x = table(perit.combined@meta.data$SingleR_RNAseq, perit.combined@meta.data$Sample_ID), margin = 2)

my.yf <- grep("YF", colnames(my.freq.table))
my.of <- grep("OF", colnames(my.freq.table))
my.ym <- grep("YM", colnames(my.freq.table))
my.om <- grep("OM", colnames(my.freq.table))

# reorder by bio group
# YF/OF/YM/OM,  NIA then JAX
my.freq.table <- my.freq.table[,c(my.yf, my.of, my.ym, my.om)]

my.freq.table.av <- apply(my.freq.table,1,mean)
my.freq.table.av.sort <- sort(my.freq.table.av, decreasing = T,index.return = T)

my.Bcell.freqs <- list("Bcell_YF" = my.freq.table["B cells",1:8  ],
                       "Bcell_OF" = my.freq.table["B cells",9:15 ],
                       "Bcell_YM" = my.freq.table["B cells",16:23],
                       "Bcell_OM" = my.freq.table["B cells",24:30])

my.Mph.freqs <- list("Mph_YF" = my.freq.table["Macrophages",1:8  ],
                     "Mph_OF" = my.freq.table["Macrophages",9:15 ],
                     "Mph_YM" = my.freq.table["Macrophages",16:23],
                     "Mph_OM" = my.freq.table["Macrophages",24:30])


my.Tcell.freqs <- list("Tcell_YF" = my.freq.table["T cells",1:8  ],
                       "Tcell_OF" = my.freq.table["T cells",9:15 ],
                       "Tcell_YM" = my.freq.table["T cells",16:23],
                       "Tcell_OM" = my.freq.table["T cells",24:30])


my.mph.age.F <- wilcox.test(my.Mph.freqs$Mph_YF, my.Mph.freqs$Mph_OF)
my.mph.age.M <- wilcox.test(my.Mph.freqs$Mph_YM, my.Mph.freqs$Mph_OM)

my.b.age.F <- wilcox.test(my.Bcell.freqs$Bcell_YF, my.Bcell.freqs$Bcell_OF)
my.b.age.M <- wilcox.test(my.Bcell.freqs$Bcell_YM, my.Bcell.freqs$Bcell_OM)

my.t.age.F <- wilcox.test(my.Tcell.freqs$Tcell_YF, my.Tcell.freqs$Tcell_OF)
my.t.age.M <- wilcox.test(my.Tcell.freqs$Tcell_YM, my.Tcell.freqs$Tcell_OM)

### for beeswarm
my.pch.ptwise <- c(rep(15,3),rep(16,5) ,  # YF
                   rep(15,3),rep(16,4) ,  # OF
                   rep(15,3),rep(16,5) ,  # YM
                   rep(15,3),rep(16,4) )  # OM

my.col.ptwise <- c(rep("black",3),rep("gray33",3),rep("black",2) ,  # YF
                   rep("black",3),rep("gray33",2),rep("black",2) ,  # OF
                   rep("black",3),rep("gray33",3),rep("black",2) ,  # YM
                   rep("black",3),rep("gray33",2),rep("black",2) )  # OM


pdf(paste(Sys.Date(),"Bcell_Tcell_Mph_freq_boxplot_SingleR_CONSISTENT_COMBINED_NIA_JAX_ALL_Cohorts.pdf",sep = "_"), height = 4, width = 8)
par(mfrow=c(1,3))
boxplot(my.Bcell.freqs,
        col = c("deeppink","deeppink4","deepskyblue","deepskyblue4"),
        las = 2, ylab = "B cell proportion (Fraction)",
        main = "B cell proportion (scRNAseq)",
        ylim = c(0,1), outline = F)
beeswarm(my.Bcell.freqs, add = T, cex = 1.5, pwpch = my.pch.ptwise, pwcol = my.col.ptwise )
text(1.5,1,signif(my.b.age.F$p.value,2))
text(3.5,1,signif(my.b.age.M$p.value,2))

boxplot(my.Mph.freqs,
        col = c("deeppink","deeppink4","deepskyblue","deepskyblue4"),
        las = 2, ylab = "Macrophage proportion (Fraction)",
        main = "Macrophage proportion (scRNAseq)",
        ylim = c(0,1), outline = F)
beeswarm(my.Mph.freqs, add = T, cex = 1.5, pwpch = my.pch.ptwise, pwcol = my.col.ptwise  )
text(1.5,1,signif(my.mph.age.F$p.value,2))
text(3.5,1,signif(my.mph.age.M$p.value,2))

boxplot(my.Tcell.freqs,
        col = c("deeppink","deeppink4","deepskyblue","deepskyblue4"),
        las = 2, ylab = "T cell  proportion (Fraction)",
        main = "T cell proportion (scRNAseq)",
        ylim = c(0,1), outline = F)
beeswarm(my.Tcell.freqs, add = T, cex = 1.5, pwpch = my.pch.ptwise, pwcol = my.col.ptwise  )
text(1.5,1,signif(my.t.age.F$p.value,2))
text(3.5,1,signif(my.t.age.M$p.value,2))

par(mfrow=c(1,1))
dev.off()

# subset and save macrophages
All.perit.singlets.Mph <- subset(perit.combined, subset = SingleR_ImmGen %in% "Macrophages")    # Macrophages, 15983 cells
save(All.perit.singlets.Mph, file = paste(Sys.Date(),"Seurat_COMBINED_10x_peritoneal_Singlets_SingleR_Consistent_Macrophages.RData",sep = "_"))


################################################################################################################
### 4. scProportionTest (global/integrated)

# create prop test object
perit.comb.prop_test <- sc_utils(perit.combined)


# Once the object is created, the permutation testing and bootstrapping can be run.
perit.comb.prop_test.F <- permutation_test(perit.comb.prop_test, 
                                           cluster_identity = "SingleR_ImmGen",
                                           sample_1 = "YF", 
                                           sample_2 = "OF",
                                           sample_identity = "Condition")

perit.comb.prop_test.M <- permutation_test(perit.comb.prop_test, 
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
perm.f <- permutation_plot_mod(perit.comb.prop_test.F, cols_vals = c("deeppink"   , "deeppink4","grey")   , my.title = "Female Aging (All cohorts)")
perm.m <- permutation_plot_mod(perit.comb.prop_test.M, cols_vals = c("deepskyblue", "deepskyblue4","grey"), my.title = "Male Aging (All cohorts)")

pdf(paste(Sys.Date(),"scProportionTest_SingleR_CONSISTENT_All_Cohort.pdf",sep = "_"), height = 4, width = 10)
perm.f + perm.m
dev.off()



################################################################################################################
### 5.  Cell Cycle analysis

# get cell type proportions
cc.freq.table <- prop.table(x = table(All.perit.singlets.Mph@meta.data$Phase, All.perit.singlets.Mph@meta.data$Sample_ID), margin = 2)

# reorder by bio group
cc.freq.table <- cc.freq.table[,c(my.yf, my.of, my.ym, my.om)]


my.G1.freqs <- list("G1_YF" = cc.freq.table["G1",1:8  ],
                    "G1_OF" = cc.freq.table["G1",9:15 ],
                    "G1_YM" = cc.freq.table["G1",16:23],
                    "G1_OM" = cc.freq.table["G1",24:30])

my.S.freqs <- list("S_YF" = cc.freq.table["S",1:8  ],
                   "S_OF" = cc.freq.table["S",9:15 ],
                   "S_YM" = cc.freq.table["S",16:23],
                   "S_OM" = cc.freq.table["S",24:30])

my.G2M.freqs <- list("G2M_YF" = cc.freq.table["G2M",1:8  ],
                     "G2M_OF" = cc.freq.table["G2M",9:15 ],
                     "G2M_YM" = cc.freq.table["G2M",16:23],
                     "G2M_OM" = cc.freq.table["G2M",24:30])

my.prolif.freqs <- list("SG2M_YF" = cc.freq.table["S",1:8  ] + cc.freq.table["G2M",1:8  ],
                        "SG2M_OF" = cc.freq.table["S",9:15 ] + cc.freq.table["G2M",9:15 ],
                        "SG2M_YM" = cc.freq.table["S",16:23] + cc.freq.table["G2M",16:23],
                        "SG2M_OM" = cc.freq.table["S",24:30] + cc.freq.table["G2M",24:30])


my.G1.age.F <- wilcox.test(my.G1.freqs$G1_YF, my.G1.freqs$G1_OF)
my.G1.age.M <- wilcox.test(my.G1.freqs$G1_YM, my.G1.freqs$G1_OM)

my.S.age.F <- wilcox.test(my.S.freqs$S_YF, my.S.freqs$S_OF)
my.S.age.M <- wilcox.test(my.S.freqs$S_YM, my.S.freqs$S_OM)

my.G2M.age.F <- wilcox.test(my.G2M.freqs$G2M_YF, my.G2M.freqs$G2M_OF)
my.G2M.age.M <- wilcox.test(my.G2M.freqs$G2M_YM, my.G2M.freqs$G2M_OM)

my.pro.age.F <- wilcox.test(my.prolif.freqs$SG2M_YF, my.prolif.freqs$SG2M_OF)
my.pro.age.M <- wilcox.test(my.prolif.freqs$SG2M_YM, my.prolif.freqs$SG2M_OM)


pdf(paste(Sys.Date(),"G1_S_G2M_freq_boxplot_MACROPHAGES_CONSISTENT_COMBINED_NIA_JAX.pdf",sep = "_"), height = 4, width = 8)
par(mfrow=c(1,3))
boxplot(my.G1.freqs,
        col = c("deeppink","deeppink4","deepskyblue","deepskyblue4"),
        las = 2, ylab = "G1 cell proportion (Fraction)",
        main = "G1 cell proportion (Mph)",
        ylim = c(0,1), outline = F)
beeswarm(my.G1.freqs, add = T, cex = 1.5, pwpch = my.pch.ptwise, pwcol = my.col.ptwise )
text(1.5,1,signif(my.G1.age.F$p.value,2))
text(3.5,1,signif(my.G1.age.M$p.value,2))

boxplot(my.S.freqs,
        col = c("deeppink","deeppink4","deepskyblue","deepskyblue4"),
        las = 2, ylab = "S proportion (Fraction)",
        main = "S proportion (Mph)",
        ylim = c(0,1), outline = F)
beeswarm(my.S.freqs, add = T, cex = 1.5, pwpch = my.pch.ptwise, pwcol = my.col.ptwise)
text(1.5,1,signif(my.S.age.F$p.value,2))
text(3.5,1,signif(my.S.age.M$p.value,2))

boxplot(my.G2M.freqs,
        col = c("deeppink","deeppink4","deepskyblue","deepskyblue4"),
        las = 2, ylab = "G2M cell  proportion (Fraction)",
        main = "G2M cell proportion (scRNAseq)",
        ylim = c(0,1), outline = F)
beeswarm(my.G2M.freqs, add = T, cex = 1.5, pwpch = my.pch.ptwise, pwcol = my.col.ptwise )
text(1.5,1,signif(my.G2M.age.F$p.value,2))
text(3.5,1,signif(my.G2M.age.M$p.value,2))

par(mfrow=c(1,1))
dev.off()


pdf(paste(Sys.Date(),"G1_versus_SG2M_freq_boxplot_MACROPHAGES_CONSISTENT_COMBINED_NIA_JAX.pdf",sep = "_"), height = 4, width = 8)
par(mfrow=c(1,3))
boxplot(my.G1.freqs,
        col = c("deeppink","deeppink4","deepskyblue","deepskyblue4"),
        las = 2, ylab = "G1 cell proportion (Fraction)",
        main = "G1 cell proportion (Mph)",
        ylim = c(0,1), outline = F)
beeswarm(my.G1.freqs, add = T, cex = 1.5, pwpch = my.pch.ptwise, pwcol = my.col.ptwise  )
text(1.5,1,signif(my.G1.age.F$p.value,2))
text(3.5,1,signif(my.G1.age.M$p.value,2))

boxplot(my.prolif.freqs,
        col = c("deeppink","deeppink4","deepskyblue","deepskyblue4"),
        las = 2, ylab = "S/G2M cell  proportion (Fraction)",
        main = "S/G2M cell proportion (scRNAseq)",
        ylim = c(0,1), outline = F)
beeswarm(my.prolif.freqs, add = T, cex = 1.5, pwpch = my.pch.ptwise, pwcol = my.col.ptwise )
text(1.5,1,signif(my.pro.age.F$p.value,2))
text(3.5,1,signif(my.pro.age.M$p.value,2))
par(mfrow=c(1,1))
dev.off()


################################################################################################################
### 6.  Augur
augur.perit <-  calculate_auc(as.matrix(perit.combined@assays$integrated@data),
                              perit.combined@meta.data,
                              cell_type_col = "SingleR_ImmGen",
                              label_col = "Condition",
                              n_threads = 1)
augur.perit
# $AUC
# # A tibble: 4 x 2
# cell_type     auc
# <chr>       <dbl>
# 1 Macrophages 0.632
# 2 Monocytes   0.624
# 3 T cells     0.590
# 4 B cells     0.561
save(augur.perit, file = paste0(Sys.Date(),"_Augur_3Cohorts_COMBINED_peritoneal_cells_object.RData"))

# for some reason, barcode column trips Augur up
perit.combined@meta.data <- perit.combined@meta.data[,colnames(perit.combined@meta.data)  != "barcode"]

pdf(paste0(Sys.Date(),"_Augur_3Cohorts_COMBINED_peritoneal_cells_UMAP.pdf"), width = 3, height = 3)
plot_umap(augur.perit,perit.combined, cell_type_col = "SingleR_ImmGen")
dev.off()

pdf(paste0(Sys.Date(),"_Augur_3Cohorts_COMBINED_peritoneal_cells_UMAP_Red_Blue.pdf"), width = 3, height = 3)
plot_umap(augur.perit,perit.combined, cell_type_col = "SingleR_ImmGen", palette = colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50))
dev.off()

pdf(paste0(Sys.Date(),"_Augur_3Cohorts_COMBINED_peritoneal_cells_Lollipop.pdf"), width = 3, height = 3)
plot_lollipop(augur.perit)
dev.off()
############################################################################################################

#######################
sink(file = paste(Sys.Date(),"_scRNAseq_peritoneal_3Cohorts_COMBINED_analysis_session_Info.txt", sep =""))
sessionInfo()
sink()

