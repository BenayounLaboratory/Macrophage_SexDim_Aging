setwd('/Volumes/BB_HQ_1/Immune_sex_dimorphism_Aging/Macrophages/Clean_Macrophage_Data/Peritoneal_scRNAseq/v3_NIA/Seurat')
options(stringsAsFactors = F)

# General use packages
library('Seurat')
library(bitops)
library(sctransform)
library(SingleR)
library(pheatmap)
library(dplyr)
library('ggplot2')

# Doublet finder
library(DoubletFinder)

# removal of ambient RNA
# https://github.com/constantAmateur/SoupX
library('SoupX')

# cell type annotation
library(SingleR)  # SingleR_1.0.1


########################################################################################################################################################
#### 0. Assumed doublet information/to calculate %age for prediction
# Targeted Cell Recovery  # of Cells Loaded	Barcodes Detected	Singlets Multiplets	Multiplet Rate
#  3,000	           4,950           ~3,000        ~2,900	   ~80	    ~2.4%
#  4,000	           6,600           ~3,900	       ~3,800	   ~140     ~3.2%
#  5,000	           8,250           ~4,800        ~4,600	   ~210     ~4.0%
#  6,000            9,900           ~5,700	       ~5,400	   ~300	    ~4.8%
#  7,000            11,550	         ~6,600	       ~6,200	   ~400	    ~5.6%
#  8,000            13,200	         ~7,500	       ~7,000	   ~510	    ~6.4%
#  9,000            14,850	         ~8,400	       ~7,700	   ~640	    ~7.2%
# 10,000            16,500	         ~9,200	       ~8,400	   ~780	    ~8.0%
# 12,000            19,800	         ~10,900	     ~9,800	   ~1,100   ~9.6%
# 14,000            23,100	         ~12,500	     ~11,000	 ~1,500   ~11.2%
# 16,000            26,400	         ~14,000	     ~12,100	 ~1,900   ~12.8%
# 18,000            29,700	         ~15,500	     ~13,100	 ~2,300   ~14.4%
# 20,000            33,000	         ~16,900	     ~14,100	 ~2,800   ~16.0%

pred.10x.dblt <- data.frame( "cell_number" = c(3000,4000,5000,6000,7000,8000,9000, 10000, 12000, 14000, 16000, 18000, 20000),
                             "dblt_rate"   = c(2.4 ,3.2 ,4.0 ,4.8 ,5.6 ,6.4 ,7.2 , 8.0  , 9.6, 11.2, 12.8, 14.4, 16.0))

pred_dblt_lm <- lm(dblt_rate ~ cell_number, data = pred.10x.dblt)

pdf(paste0(Sys.Date(), "_10x_cell_number_vs_doublet_rate.pdf"))
plot(dblt_rate ~ cell_number, data = pred.10x.dblt)
abline(pred_dblt_lm, col = "red", lty = "dashed")
dev.off()
########################################################################################################################################################


########################################################################################################################################################
######## 1. Read data from CellRanger, perform background removal

# Calculate and clean the contribution of ambient RNA with SoupX

# read 10x libraries cell ranger gene barcode matrices for SoupX
cellRanger.NIA_PL.1  <- load10X("/Volumes/BB_HQ_1/Immune_sex_dimorphism_Aging/Macrophages/Clean_Macrophage_Data/Peritoneal_scRNAseq/v3_NIA/Cell_Ranger/NIA_perit_p1/outs/")

# estimate RNA soup from bottom 10% droplets
cellRanger.NIA_PL.1 <- autoEstCont(cellRanger.NIA_PL.1) # Estimated global rho of 0.03

# adjust counts based on RNA soup (and round)
out.NIA_PL.1 <- adjustCounts(cellRanger.NIA_PL.1, roundToInt=TRUE)

# get seurat objects
seurat.NIA_PL.1 <- CreateSeuratObject( out.NIA_PL.1 )

# Rename cells for merge cleaned seurat objects
my.NIA.perit <- RenameCells(seurat.NIA_PL.1, add.cell.id = "NIA_pool1")
my.NIA.perit
# An object of class Seurat 
# 32285 features across 16024 samples within 1 assay 
# Active assay: RNA (32285 features, 0 variable features)

# create batch label
Batch <- rep("NA", length(colnames(my.NIA.perit@assays$RNA)))
Batch[ grep("NIA_pool1"  , colnames(my.NIA.perit@assays$RNA)) ]   <- "NIA_v3"
Batch <- data.frame(Batch)
rownames(Batch) <- colnames(my.NIA.perit@assays$RNA)

# update Seurat with metadata
my.NIA.perit <- AddMetaData(object = my.NIA.perit, metadata = as.vector(Batch)    , col.name = "Batch"       )

################################################################################################################################################################
#### 2. QC on mitochondrial reads and depth
# The number of genes and UMIs (nGene and nUMI) are automatically calculated for every object by Seurat.  
# The % of UMI mapping to MT-genes is a common scRNA-seq QC metric.
### QC on mitochondrial reads
my.NIA.perit[["percent.mito"]] <- PercentageFeatureSet(my.NIA.perit, pattern = "^mt-")

pdf(paste(Sys.Date(),"Peritoneal_cells_aging_violinPlots_QC_gene_UMI_mito.pdf", sep = "_"), height = 5, width = 10)
VlnPlot(object = my.NIA.perit, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, group.by = "Batch")
dev.off()

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(my.NIA.perit, feature1 = "nCount_RNA", feature2 = "percent.mito")
plot2 <- FeatureScatter(my.NIA.perit, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

pdf(paste(Sys.Date(),"Peritoneal_cells_aging_QC_scatter.pdf", sep = "_"), height = 5, width = 10)
plot1 + plot2
dev.off()

# We filter out cells that have unique gene counts or less than 300  
my.NIA.perit <- subset(my.NIA.perit, subset = nFeature_RNA > 300 & percent.mito < 25 & nFeature_RNA < 6000)
my.NIA.perit
# An object of class Seurat 
# 32285 features across 15449 samples within 1 assay 
# Active assay: RNA (32285 features, 0 variable features)

plot1 <- FeatureScatter(my.NIA.perit, feature1 = "nCount_RNA", feature2 = "percent.mito")
plot2 <- FeatureScatter(my.NIA.perit, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

pdf(paste(Sys.Date(),"Peritoneal_cells_aging_QC_scatter_post_filter.pdf", sep = "_"), height = 5, width = 10)
plot1 + plot2
dev.off()

################################################################################################################################################################
#### 3. Normalizing the data
# global-scaling normalization method ???LogNormalize??? that normalizes the gene expression measurements for each cell by the total expression, 
# multiplies this by a scale factor (10,000 by default), and log-transforms the result.
my.NIA.perit <- NormalizeData(object = my.NIA.perit, normalization.method = "LogNormalize",  scale.factor = 10000)

################################################################################################################################################################
#### 4. Cell cycle regression
# Read in a list of cell cycle markers, from Tirosh et al, 2015
cc.genes <- readLines(con = "regev_lab_cell_cycle_genes.txt")

# make into mouse gene names
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

cc.genes.mouse <- firstup(tolower(cc.genes))


# We can segregate this list into markers of G2/M phase and markers of S
# phase
s.genes <- cc.genes.mouse[1:43]
g2m.genes <- cc.genes.mouse[44:97]

# Assign Cell-Cycle Scores
my.NIA.perit <- CellCycleScoring(object = my.NIA.perit, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# view cell cycle scores and phase assignments
head(x = my.NIA.perit@meta.data)
#                                 orig.ident nCount_RNA nFeature_RNA  Batch percent.mito      S.Score    G2M.Score Phase     old.ident
# NIA_pool1_AAACCCAAGAAGCCTG-1 SeuratProject       4640         1342 NIA_v3     3.491379  0.049625054  0.02250210     S SeuratProject
# NIA_pool1_AAACCCAAGTGCCGAA-1 SeuratProject       8783         2614 NIA_v3     3.324604  0.002481749 -0.06771357     S SeuratProject
# NIA_pool1_AAACCCAAGTTCCGTA-1 SeuratProject      20954         3779 NIA_v3     2.624797 -0.059331138 -0.06347658    G1 SeuratProject
# NIA_pool1_AAACCCACACAACCGC-1 SeuratProject      22852         3991 NIA_v3     1.325923 -0.030191638 -0.08483664    G1 SeuratProject
# NIA_pool1_AAACCCACACGATAGG-1 SeuratProject       2278         1051 NIA_v3     3.116769 -0.035940704 -0.06762239    G1 SeuratProject
# NIA_pool1_AAACCCACATACACCA-1 SeuratProject       3456         1297 NIA_v3     1.880787  0.074786698 -0.03368829     S SeuratProject

################################################################################################################################################################
##### 5.Find and remove doublets 

# Normalize data for further processing
my.NIA.perit <- SCTransform(object = my.NIA.perit, vars.to.regress = c("nFeature_RNA", "percent.mito"))
my.NIA.perit <- RunPCA(my.NIA.perit)

# Determine the ‘dimensionality’ of the dataset
pdf(paste0(Sys.Date(), "_Peritoneal_Lavage_NIAv3_elbowplot.pdf"))
ElbowPlot(my.NIA.perit, ndims = 50)
dev.off()

################################################################################
# https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html
# To give us an idea of the number of PCs needed to be included:
# We can calculate where the principal components start to elbow by taking the larger value of:
#    - The point where the principal components only contribute 5% of standard deviation
#    - The principal components cumulatively contribute 90% of the standard deviation.
#    - The point where the percent change in variation between the consecutive PCs is less than 0.1%.

# Determine percent of variation associated with each PC
pct <- my.NIA.perit[["pca"]]@stdev / sum(my.NIA.perit[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
co1 # 41

# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

# last point where change of % of variation is more than 0.1%.
co2 # 21

# Minimum of the two calculation
pcs <- min(co1, co2)
pcs # 21

# Based on these metrics, first 16 PCs to generate the clusters.
# We can plot the elbow plot again and overlay the information determined using our metrics:

# Create a dataframe with values
plot_df <- data.frame(pct  = pct,
                      cumu = cumu,
                      rank = 1:length(pct))

# Elbow plot to visualize
pdf(paste0(Sys.Date(), "_Peritoneal_Lavage_NIA_elbowplot_threshold_analysis.pdf"), height = 5, width= 6)
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) +
  geom_text() +
  geom_vline(xintercept = 90, color = "grey") +
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()
dev.off()
###############################################################################


my.NIA.perit <- RunUMAP(my.NIA.perit, dims = 1:pcs)
my.NIA.perit <- FindNeighbors(my.NIA.perit, dims = 1:pcs)
my.NIA.perit <- FindClusters(object = my.NIA.perit)


#############################################################################################################
################################  A. Doublet Finder [based on gene expression]  ################################ 
#############################################################################################################

#### need to split by original 10x sample to make sure to identify real doublets when running doublet finding

## Assume doublet rate based on 10x information (+ 20% cells to be safe)
pred.dblt.rate <- predict(pred_dblt_lm, data.frame("cell_number" = 1.2 * ncol(my.NIA.perit)))/100
pred.dblt.rate
# 1 
# 0.1483104 

## pK Identification (no ground-truth)
sweep.res.list_perit <- paramSweep_v3(my.NIA.perit, PCs = 1:30, sct = TRUE, num.cores	 = 4)
sweep.stats_perit    <- summarizeSweep(sweep.res.list_perit , GT = FALSE)
bcmvn_perit          <- find.pK(sweep.stats_perit)

# need some R gymnastics since the Pk is stored as a factor for some reason
# to get the pK number, need to first convert to character and THEN to numeric
# numeric first yeild row number
pk.perit <- as.numeric(as.character(bcmvn_perit[as.numeric(bcmvn_perit$pK[bcmvn_perit$BCmetric == max(bcmvn_perit$BCmetric)]),"pK"]))

## Homotypic Doublet Proportion Estimate
homotypic.prop <- modelHomotypic(my.NIA.perit@meta.data$seurat_clusters)     ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi       <- round(pred.dblt.rate*nrow(my.NIA.perit@meta.data))        ## Assuming doublet formation rate based on 10x manual

## Run DoubletFinder with varying classification stringencies
my.NIA.perit <- doubletFinder_v3(my.NIA.perit, PCs = 1:30, pN = 0.25, pK = pk.perit, nExp = nExp_poi,     reuse.pANN = FALSE, sct = T)

# rename column to enable subsetting
colnames(my.NIA.perit@meta.data)[grep("DF.classifications_0.25",colnames(my.NIA.perit@meta.data))] <- "DoubletFinder"


# plot UMAP plots 
# get classification name
pdf(paste(Sys.Date(),"NIAv3_Peritoneal_NIAv3Doublet_Finder_UMAP.pdf", sep = "_"), height = 5, width = 5)
print(DimPlot(my.NIA.perit, reduction = "umap", group.by = "DoubletFinder"))
dev.off()


table(my.NIA.perit@meta.data$DoubletFinder)
# Doublet Singlet 
#     2291   13158


save(my.NIA.perit, file = paste0(Sys.Date(),"_NIAv3_Peritoneal_cells_aging_Seurat_object_DoubletFinderList.RData"))


#############################################################################################################
##############################  A. HTO analysis [based on TotalSeq antibodies]  ############################# 
#############################################################################################################

#################################################################
############### a. read in cite-seq tools results ###############
#################################################################

################################ read HTO libraries 1
HTO_p1 <- Read10X('/Volumes/BB_HQ_1/Immune_sex_dimorphism_Aging/Macrophages/Clean_Macrophage_Data/Peritoneal_scRNAseq/v3_NIA/HTO_processing/NIA_perit_p1_HTO_parsed/umi_count', gene.column=1)
HTO_p1 <- HTO_p1[1:10,] # remove unmapped

# need to clean up cells with no HTOs
my.zero <- which(colSums(as.matrix(HTO_p1)) < 100)
HTO_p1 <- HTO_p1[,-my.zero] # remove cells with no hashtag
dim(HTO_p1) #  10 20390


################################################################
###############  b. Filter and demultiplex data  ###############
################################################################

# extract barcode information
my.NIA.perit <- AddMetaData(object = my.NIA.perit, metadata = as.vector(unlist(lapply(strsplit(rownames(my.NIA.perit@meta.data),"_|-"),'[[',3))), col.name = "barcode")

################################ process Pool #1 ################################
# Select cell barcodes detected by both RNA and HTO
joint.bcs.1 <- intersect(colnames(HTO_p1), my.NIA.perit@meta.data$barcode) # 5090

# Subset RNA Suerat object by joint cell barcodes
NIA.p1.umis <- my.NIA.perit[, my.NIA.perit@meta.data$barcode  %in% joint.bcs.1]

# Subset HTO counts by joint cell barcodes (correct order from subsetted seurat above)
NIA.p1.htos <- as.matrix(HTO_p1[, NIA.p1.umis@meta.data$barcode])

# make names that will correspond to Seurat cell names
colnames(NIA.p1.htos) <- rownames(NIA.p1.umis@meta.data)

# Confirm that the HTO have the correct names
rownames(NIA.p1.htos) 
# [1] "HTO_01_F_4m_800-ACCCACCAGTAAGACC"  "HTO_02_M_4m_810-GGTCGAGAGCATTCAA"  "HTO_03_F_20m_805-CTTGCCGCATGTCATT" "HTO_04_M_20m_815-AAAGCATTCTTCACGG"
# [5] "HTO_05_M_20m_816-CTTTGTCTTTGTGAGG" "HTO_06_F_20m_806-TATGCTGCCACGGTAA" "HTO_07_M_4m_811-GAGTCTGCCAGTATCC"  "HTO_08_F_4m_801-TATAGAACGCCAGGCC" 
# [9] "HTO_09_F_4m_802-TGCCTATGAAACAAGG"  "HTO_10_M_4m_812-CCGATTGTAACAGACC" 

# create a Seurat object copy
NIA.p1.hashtag <- NIA.p1.umis

# Normalize RNA data with log normalization
NIA.p1.hashtag <- NormalizeData(NIA.p1.hashtag, normalization.method = "LogNormalize",  scale.factor = 10000)

# Find and scale variable features
NIA.p1.hashtag <- FindVariableFeatures(NIA.p1.hashtag, selection.method = "mean.var.plot")
NIA.p1.hashtag <- ScaleData(NIA.p1.hashtag, features = VariableFeatures(NIA.p1.hashtag))

# Add HTO data as a new assay independent from RNA
NIA.p1.hashtag[["HTO"]] <- CreateAssayObject(counts = NIA.p1.htos)

# Normalize HTO data, here we use centered log-ratio (CLR) transformation
NIA.p1.hashtag <- NormalizeData(NIA.p1.hashtag, assay = "HTO", normalization.method = "CLR")

# Demultiplex cells based on HTO enrichment
NIA.p1.hashtag <- MULTIseqDemux(NIA.p1.hashtag, assay = "HTO") # https://github.com/satijalab/seurat/issues/2549

# Global classification results
table(NIA.p1.hashtag$MULTI_ID)
# Doublet  HTO-01-F-4m-800-ACCCACCAGTAAGACC  HTO-02-M-4m-810-GGTCGAGAGCATTCAA HTO-03-F-20m-805-CTTGCCGCATGTCATT 
#    2471                               794                              1842                               695 
# HTO-04-M-20m-815-AAAGCATTCTTCACGG HTO-05-M-20m-816-CTTTGTCTTTGTGAGG HTO-06-F-20m-806-TATGCTGCCACGGTAA  HTO-07-M-4m-811-GAGTCTGCCAGTATCC 
#                              1170                              1530                              1393                              1160 
# HTO-08-F-4m-801-TATAGAACGCCAGGCC  HTO-09-F-4m-802-TGCCTATGAAACAAGG  HTO-10-M-4m-812-CCGATTGTAACAGACC                          Negative 
#                              1915                              1165                               808                               329 

# Visualize enrichment for selected HTOs with ridge plots, and group cells based on the max HTO signal
Idents(NIA.p1.hashtag) <- "MULTI_ID"

pdf(paste0(Sys.Date(),"_NIA_Pool1_HTO_id_ridge_plot.pdf"), height = 10, width = 15)
RidgePlot(NIA.p1.hashtag, assay = "HTO", features = rownames(NIA.p1.hashtag[["HTO"]]), ncol = 2)
dev.off()

pdf(paste0(Sys.Date(),"_NIA_Pool1_HTO1_3_feature_Scatter_plot.pdf"), height = 5, width = 10)
FeatureScatter(NIA.p1.hashtag, feature1 = "HTO-01-F-4m-800-ACCCACCAGTAAGACC", feature2 = "HTO-03-F-20m-805-CTTGCCGCATGTCATT")
dev.off()

# Compare number of UMIs for singlets, doublets and negative cells
pdf(paste0(Sys.Date(),"_NIA_Pool1_nCount_RNA_violin_based_on_HTO_Demux.pdf"), height = 10, width = 15)
VlnPlot(NIA.p1.hashtag, features = "nCount_RNA", pt.size = 0.1, log = TRUE)
dev.off()

# Generate a two dimensional UMAP embedding for HTOs (grouping cells by singlets and doublets for simplicity)
# First, we will remove HTO-negative cells from the object (since they can't be attributed anyway)
NIA.p1.hashtag.subset <- subset(NIA.p1.hashtag, idents = "Negative", invert = TRUE)

# Calculate a UMAP embedding of the HTO data
DefaultAssay(NIA.p1.hashtag.subset) <- "HTO"
NIA.p1.hashtag.subset <- ScaleData(NIA.p1.hashtag.subset, features = rownames(NIA.p1.hashtag.subset), verbose = FALSE)
NIA.p1.hashtag.subset <- RunPCA(NIA.p1.hashtag.subset, features = rownames(NIA.p1.hashtag.subset), approx = FALSE)
NIA.p1.hashtag.subset <- RunUMAP(NIA.p1.hashtag.subset, dims = 1:10, perplexity = 100, check_duplicates = FALSE)

pdf(paste0(Sys.Date(),"_NIA_Pool1_UMAP_hashtag.pdf"), height = 6, width = 10)
DimPlot(NIA.p1.hashtag.subset)
dev.off()

# create metadata column to plot HTO heatmap
NIA.p1.hashtag@meta.data$HTO.global <- rep("Singlet",dim(NIA.p1.hashtag@meta.data)[1])
NIA.p1.hashtag@meta.data$HTO.global[NIA.p1.hashtag@meta.data$MULTI_ID %in% "Doublet"] <- "Doublet"
NIA.p1.hashtag@meta.data$HTO.global[NIA.p1.hashtag@meta.data$MULTI_ID %in% "Negative"] <- "Negative"

# Create an HTO heatmap, based on Figure 1C in the Cell Hashing paper.
pdf(paste0(Sys.Date(),"_NIA_Pool1_HTO_heatmap.pdf"), height = 5, width = 9)
HTOHeatmap(NIA.p1.hashtag, assay = "HTO",
           classification = "MULTI_ID",
           global.classification = "HTO.global")
dev.off()

save(NIA.p1.hashtag, file = paste0(Sys.Date(),"_NIA_Pool1_post_HTO_Demux.RData"))


#############################################################################
###############            c. Compare doublets              ###############
#############################################################################

# combined the 3 objects
NIA.combined.calls <- NIA.p1.hashtag
NIA.combined.calls
# An object of class Seurat 
# 48925 features across 15272 samples within 3 assays 
# Active assay: SCT (16630 features, 564 variable features)
# 2 other assays present: RNA, HTO
# 2 dimensional reductions calculated: pca, umap

# Compare HTO and scds calls
table(NIA.combined.calls@meta.data$DoubletFinder,NIA.combined.calls@meta.data$HTO.global)
#          Doublet Negative Singlet
# Doublet     838       68    1311
# Singlet    1633      261   11161

# Keep cells marked as singlets from both methods
NIA.singlets <- NIA.combined.calls[, bitAnd(NIA.combined.calls@meta.data$DoubletFinder == "Singlet", NIA.combined.calls@meta.data$HTO.global == "Singlet")>0]
NIA.singlets
# An object of class Seurat 
# 48925 features across 11161 samples within 3 assays 
# Active assay: SCT (16630 features, 564 variable features)
# 2 other assays present: RNA, HTO
# 2 dimensional reductions calculated: pca, umap

# table(NIA.singlets@meta.data$DoubletFinder,NIA.singlets@meta.data$HTO.global) ### ALL singlets, sanity check
save(NIA.singlets, file = paste0(Sys.Date(),"_NIA_Combined_Singlets_ONLY_Seurat_object.RData"))


################################################################################################################################################################
##### 6.Scaling the data and removing unwanted sources of variation

# Run SCT normalization, regress potential problems
NIA.singlets <- SCTransform(object = NIA.singlets, vars.to.regress = c("nFeature_RNA", "percent.mito", "S.Score", "G2M.Score"))

# find variable genes
NIA.singlets <- FindVariableFeatures(NIA.singlets, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(NIA.singlets), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(NIA.singlets)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

pdf(paste(Sys.Date(),"NIAv3_Peritoneal_cells_aging_Singlets_Variable genes.pdf", sep = "_"), height = 5, width = 10)
plot1 + plot2
dev.off()


################################################################################################################################################################
##### 7. Update metadata information based on HTO and output information on cell cycling

# drop doublet levels
NIA.singlets$MULTI_classification <- droplevels(NIA.singlets$MULTI_classification)

table(NIA.singlets$MULTI_classification)
# HTO-01-F-4m-800-ACCCACCAGTAAGACC  HTO-02-M-4m-810-GGTCGAGAGCATTCAA HTO-03-F-20m-805-CTTGCCGCATGTCATT HTO-04-M-20m-815-AAAGCATTCTTCACGG 
# 686                              1633                               647                              1097 
# HTO-05-M-20m-816-CTTTGTCTTTGTGAGG HTO-06-F-20m-806-TATGCTGCCACGGTAA  HTO-07-M-4m-811-GAGTCTGCCAGTATCC  HTO-08-F-4m-801-TATAGAACGCCAGGCC 
# 1378                              1341                              1012                              1670 
# HTO-09-F-4m-802-TGCCTATGAAACAAGG  HTO-10-M-4m-812-CCGATTGTAACAGACC 
# 1010                               687 

#######################
# create condition label
my.YF <- grep("F-4m" , NIA.singlets$MULTI_classification)
my.OF <- grep("F-20m", NIA.singlets$MULTI_classification)
my.YM <- grep("M-4m" , NIA.singlets$MULTI_classification)
my.OM <- grep("M-20m", NIA.singlets$MULTI_classification)

Condition <- rep("NA", length(colnames(NIA.singlets@assays$RNA)))
Condition[ my.YF  ]   <- "YF"
Condition[ my.OF  ]   <- "OF"
Condition[ my.YM  ]   <- "YM"
Condition[ my.OM  ]   <- "OM"

Condition <- data.frame(Condition)
rownames(Condition) <- colnames(NIA.singlets@assays$RNA)

NIA.singlets <- AddMetaData(object = NIA.singlets, metadata = as.vector(Condition), col.name = "Condition")


#######################
# create sample label
my.YF1 <- grep("F-4m-800" , NIA.singlets$MULTI_classification)
my.YF2 <- grep("F-4m-801" , NIA.singlets$MULTI_classification)
my.YF3 <- grep("F-4m-802" , NIA.singlets$MULTI_classification)
my.OF1 <- grep("F-20m-805", NIA.singlets$MULTI_classification)
my.OF2 <- grep("F-20m-806", NIA.singlets$MULTI_classification)
my.YM1 <- grep("M-4m-810" , NIA.singlets$MULTI_classification)
my.YM2 <- grep("M-4m-811" , NIA.singlets$MULTI_classification)
my.YM3 <- grep("M-4m-812" , NIA.singlets$MULTI_classification)
my.OM1 <- grep("M-20m-815", NIA.singlets$MULTI_classification)
my.OM2 <- grep("M-20m-816", NIA.singlets$MULTI_classification)

Sample <- rep("NA", length(colnames(NIA.singlets@assays$RNA)))
Sample[ my.YF1  ]   <- "NIA_YF1"
Sample[ my.YF2  ]   <- "NIA_YF2"
Sample[ my.YF3  ]   <- "NIA_YF3"
Sample[ my.OF1  ]   <- "NIA_OF1"
Sample[ my.OF2  ]   <- "NIA_OF2"
Sample[ my.YM1  ]   <- "NIA_YM1"
Sample[ my.YM2  ]   <- "NIA_YM2"
Sample[ my.YM3  ]   <- "NIA_YM3"
Sample[ my.OM1  ]   <- "NIA_OM1"
Sample[ my.OM2  ]   <- "NIA_OM2"

Sample <- data.frame(Sample)
rownames(Sample) <- colnames(NIA.singlets@assays$RNA)

NIA.singlets <- AddMetaData(object = NIA.singlets, metadata = as.vector(Sample), col.name = "Sample_ID")

#### check
head(NIA.singlets@meta.data)


#### Parse Cell Cyle information

# write predictions to file
write.table(NIA.singlets@meta.data, file = paste(Sys.Date(),"Peritoneal_cells_aging_v3NIA_CellCycle_predictions.txt", sep = "_"), sep = "\t", quote = F)


pdf(paste0(Sys.Date(), "Peritoneal_aging_v3NIA_cell_cycle_pie_charts.pdf"), height = 13, width = 8)
par(mfrow=c(2,2))
par(oma=c(0.1,0.1,0.1,0.1))
pie(table(NIA.singlets@meta.data$Phase[NIA.singlets@meta.data$Condition == "YF"]) , main = "YF")
pie(table(NIA.singlets@meta.data$Phase[NIA.singlets@meta.data$Condition == "YM"]) , main = "YM")
pie(table(NIA.singlets@meta.data$Phase[NIA.singlets@meta.data$Condition == "OF"]) , main = "OF")
pie(table(NIA.singlets@meta.data$Phase[NIA.singlets@meta.data$Condition == "OM"]) , main = "OM")
par(mfrow=c(1,1))
dev.off()


################################################################################################################################################################
##### 8. Cluster the cells

# rerun dimensionality reduction on merged/clean object
NIA.singlets <- RunPCA(NIA.singlets)
NIA.singlets <- RunUMAP(NIA.singlets, dims = 1:21)

# look at samples
pdf(paste(Sys.Date(),"Peritoneal_cells_aging_v3NIA_Singlets_UMAP_by_Group.pdf", sep = "_"), height = 5, width = 5)
DimPlot(NIA.singlets, reduction = "umap", group.by = "Condition")
DimPlot(NIA.singlets, reduction = "umap", group.by = "Sample_ID")
dev.off()

# perform clustering again now that only singlets remain
NIA.singlets <- FindNeighbors(NIA.singlets, 
                              dims = 1:21)
NIA.singlets <- FindClusters(object = NIA.singlets,
                             reduction = "pca",
                             k.param = 100,
                             dims = 1:21,
                             resolution = 0.3,
                             n.start = 100)
# Number of communities: 12

pdf(paste(Sys.Date(),"NIAv3_Peritoneal_cells_aging_UMAP_clusters_0.3res_Singlets.pdf", sep = "_"), height = 5, width = 6.5)
DimPlot(NIA.singlets, label = TRUE)
dev.off()

write.table(NIA.singlets@meta.data, file = paste(Sys.Date(),"NIAv3_Peritoneal_cells_aging_Seurat_PCA_SNN_clusters_and_Metadata_Singlets.txt", sep = "_"), sep = "\t", quote = F, col.names = F, row.names = T)

pdf(paste(Sys.Date(),"NIAv3_Peritoneal_cells_aging_Feature_plot_Singlets_UMI_counts.pdf", sep = "_"),  height = 5, width = 6.5)
FeaturePlot(object = NIA.singlets, features = c("nFeature_RNA")  )
dev.off()

################################################################################################################################################################
##### 9. Visualization of potential marker gene expression

### key markers
pdf(paste(Sys.Date(),"NIAv3_Peritoneal_cells_aging_Feature_plot_Singlets.pdf", sep = "_"), width = 15, height = 15)
FeaturePlot(object = NIA.singlets, features = c("Cd79a","Adgre1","Itgam","Bhlhe40","Cd14","Cd3e","Spi1","Cd209a")  )
FeaturePlot(object = NIA.singlets, features = c("Cd86", "Cd68", "Cd80","Tlr2", "Tlr4")  ) ### M1 "Cx3cr1","Nos2","Il1r", "Il1b"
FeaturePlot(object = NIA.singlets, features = c("Arg1", "Cd163")  ) # m2 "Cd206","Fizz1" Ym1/2
dev.off()

pdf(paste(Sys.Date(),"NIAv3_Peritoneal_cells_aging_Xist_ridge_plot_Singlets.pdf", sep = "_"))
RidgePlot(object = NIA.singlets, features = "Xist", ncol = 1, group.by = "Condition")
dev.off()

pdf(paste(Sys.Date(),"NIAv3_Peritoneal_cells_aging_Y_genes_ridge_plot_Singlets.pdf", sep = "_"), height = 5, width = 10)
RidgePlot(object = NIA.singlets, features = c("Ddx3y", "Eif2s3y"), ncol = 2, group.by = "Condition")
dev.off()

save(NIA.singlets, file = paste(Sys.Date(),"NIAv3_Peritoneal_10X_Singlets.Seurat.object.RData",sep = "_"))


################################################################################################################################################################
##### 10. Identify cell types using SingleR

################ run singleR and save object
my.singler = CreateSinglerObject(counts = NIA.singlets@assays$SCT@counts,
                                 annot = NIA.singlets@meta.data$Sample_ID,
                                 min.genes = 250, technology = "10X",
                                 species = "Mouse", project.name = "10X_Peritoneal")

my.singler$seurat               = NIA.singlets # (optional)
my.singler$meta.data$orig.ident = NIA.singlets@meta.data$orig.ident # the original identities, if not supplied in 'annot'

## if using Seurat v3.0 and over use:
my.singler$meta.data$xy       =    NIA.singlets@reductions$umap@cell.embeddings # the Umap coordinates
my.singler$meta.data$clusters =    NIA.singlets@active.ident                    # the Seurat clusters (if 'clusters' not provided)
my.singler$meta.data$Sample   =    NIA.singlets@meta.data$Sample                # the Seurat clusters (if 'clusters' not provided)

save(my.singler, file = paste(Sys.Date(),"NIAv3_Peritoneal_10X_SingleR_object_SCT.RData",sep = "_"))

#######################
### Main cell types
pdf(paste(Sys.Date(),"NIAv3_Peritoneal_10X_SingleR_object.main.cell.types.pdf",sep = "_"), height = 5, width = 7)
SingleR.DrawHeatmap(my.singler$singler[[1]]$SingleR.single.main, top.n = 50, clusters = my.singler$meta.data$Sample)
dev.off()

### All cell types
pdf(paste(Sys.Date(),"NIAv3_Peritoneal_10X_SingleR_object.all.cell.types.pdf",sep = "_"), height = 5, width = 7)
SingleR.DrawHeatmap(my.singler$singler[[1]]$SingleR.single, top.n = 50, clusters = my.singler$meta.data$Sample)
dev.off()

### Main cell types
pdf(paste(Sys.Date(),"NIAv3_Peritoneal_10X_SingleR_object.main.cell.types_MouseRNAseq.pdf",sep = "_"), height = 5, width = 7)
SingleR.DrawHeatmap(my.singler$singler[[2]]$SingleR.single.main, top.n = Inf,clusters = my.singler$meta.data$Sample)
dev.off()

### All cell types
pdf(paste(Sys.Date(),"NIAv3_Peritoneal_10X_SingleR_object.all.cell.types_MouseRNAseq.pdf",sep = "_"), height = 5, width = 7)
SingleR.DrawHeatmap(my.singler$singler[[2]]$SingleR.single, top.n = 50, clusters = my.singler$meta.data$Sample)
dev.off()

#Next, we can use the fine-tuned labels to color the t-SNE plot:
out = SingleR.PlotTsne(my.singler$singler[[1]]$SingleR.single.main,
                       my.singler$meta.data$xy,
                       do.label=FALSE,
                       do.letters = F,
                       labels = my.singler$singler[[1]]$SingleR.single.main$labels,
                       label.size = 4, dot.size = 1,
                       colors = c("blue3",
                                  "dodgerblue",
                                  "cyan3",
                                  "chartreuse3",
                                  "deeppink2",
                                  "darkgoldenrod1",
                                  "coral1",
                                  "aquamarine2",
                                  "firebrick1",
                                  "darkviolet",
                                  "darkorchid1",
                                  "deepskyblue",
                                  "darkgreen",
                                  "yellow",
                                  "darkgoldenrod",
                                  "firebrick4",
                                  "dodgerblue4"
                       ))

pdf(paste(Sys.Date(),"NIAv3_Peritoneal_10X_SingleR_UMAP_annotated_IMMGEN.pdf",sep = "_"), height = 10, width = 12)
out$p
dev.off()


out = SingleR.PlotTsne(my.singler$singler[[1]]$SingleR.single.main,
                       my.singler$meta.data$xy,
                       do.label=FALSE,
                       do.letters = F,
                       labels = my.singler$singler[[2]]$SingleR.single.main$labels,
                       label.size = 4, dot.size = 1,
                       colors = c("blue3",
                                  "dodgerblue",
                                  "cyan3",
                                  "chartreuse3",
                                  "deeppink2",
                                  "darkgoldenrod1",
                                  "coral1",
                                  "aquamarine2",
                                  "firebrick1",
                                  "darkviolet",
                                  "darkorchid1",
                                  "deepskyblue",
                                  "darkgreen",
                                  "yellow",
                                  "darkgoldenrod",
                                  "firebrick4",
                                  "dodgerblue4"
                       ))

pdf(paste(Sys.Date(),"NIAv3_Peritoneal_10X_SingleR_UMAP_annotated_MouseRNAseq.pdf",sep = "_"), height = 10, width = 12)
out$p
dev.off()

write.table(my.singler$singler[[1]]$SingleR.single.main$labels,
            file = paste(Sys.Date(),"NIAv3_Peritoneal_10X_SingleR_annotated_cells_Immgen_Main.txt",sep = "_"), sep = "\t", quote = F)
write.table(my.singler$singler[[1]]$SingleR.single$labels,
            file = paste(Sys.Date(),"NIAv3_Peritoneal_10X_SingleR_annotated_cells_Immgen_All.txt",sep = "_"), sep = "\t", quote = F)


write.table(my.singler$singler[[2]]$SingleR.single.main$labels,
            file = paste(Sys.Date(),"NIAv3_Peritoneal_10X_SingleR_annotated_cells_MouseRNAseq_Main.txt",sep = "_"), sep = "\t", quote = F)
write.table(my.singler$singler[[2]]$SingleR.single$labels,
            file = paste(Sys.Date(),"NIAv3_Peritoneal_10X_SingleR_annotated_cells_MouseRNAseq_All.txt",sep = "_"), sep = "\t", quote = F)

# Transfer SingleR annotations to Seurat object
sum(rownames(NIA.singlets@meta.data) == rownames(my.singler$singler[[1]]$SingleR.single$labels)) # 11161
dim(NIA.singlets@meta.data) # 11161    24

NIA.singlets@meta.data$SingleR_ImmGen          <- my.singler$singler[[1]]$SingleR.single.main$labels
NIA.singlets@meta.data$SingleR_ImmGen_detailed <- my.singler$singler[[1]]$SingleR.single$labels
NIA.singlets@meta.data$SingleR_RNAseq          <- my.singler$singler[[2]]$SingleR.single.main$labels
NIA.singlets@meta.data$SingleR_RNAseq_detailed <- my.singler$singler[[2]]$SingleR.single$labels

save(NIA.singlets, file = paste(Sys.Date(),"NIAv3_Peritoneal_10X_Singlets.Seurat_SingleR_annotated.object.RData",sep = "_"))

# How does cluster membership vary by replicate?
table(NIA.singlets@meta.data$SingleR_ImmGen, NIA.singlets@meta.data$Sample)
#              NIA_OF1 NIA_OF2 NIA_OM1 NIA_OM2 NIA_YF1 NIA_YF2 NIA_YF3 NIA_YM1 NIA_YM2 NIA_YM3
# B cells         398     859     324     662     269     636     418     368     371     225
# Basophils         0       0       0       0       1       0       0       0       0       0
# DC                3       3       8      11       4      16       4       4       6       2
# ILC               1       1       6       8       3       5       3       5       3       1
# Macrophages     210     395     537     582     386     935     523    1201     530     410
# Mast cells        0       0       0       8       0       2       2       1       0       5
# Monocytes         4       4       6      12       5      23      15      15      14       9
# NK cells          2       3       4       1       6       9       1       7      10       5
# NKT               0       0       0       0       0       3       2       0       1       0
# T cells          25      63     181      81      10      36      40      25      73      25
# Tgd               4      13      31      13       2       5       2       7       4       5

my.freq.table <- prop.table(x = table(NIA.singlets@meta.data$SingleR_ImmGen, NIA.singlets@meta.data$Sample), margin = 2)
my.freq.table <- my.freq.table[,c(5:7,1:2,8:10,3:4)]
my.freq.table.av <- apply(my.freq.table,1,mean)
my.freq.table.av.sort <- sort(my.freq.table.av, decreasing = T,index.return = T)

my.cols <- c("blue3",
             "dodgerblue",
             "cyan3",
             "chartreuse3",
             "deeppink2",
             "darkgoldenrod1",
             "coral1",
             "aquamarine2",
             "firebrick1",
             "darkviolet",
             "darkorchid1",
             "deepskyblue",
             "darkgreen",
             "yellow",
             "darkgoldenrod",
             "firebrick4",
             "dodgerblue4")

pdf(paste(Sys.Date(),"NIAv3_Peritoneal_10X_SingleR_IMMGEN_calls_barplot_sorted_Singlets.pdf",sep = "_"), height = 10, width = 15)
barplot(my.freq.table[my.freq.table.av.sort$ix,], las = 2,
        legend.text = rownames(my.freq.table)[my.freq.table.av.sort$ix],
        col = my.cols[my.freq.table.av.sort$ix],
        ylab = "Cell frequency in 10xGenomics single cell library")
box()
barplot(my.freq.table[my.freq.table.av.sort$ix,], las = 2,
        col = my.cols[my.freq.table.av.sort$ix],
        ylab = "Cell frequency in 10xGenomics single cell library")
box()
dev.off()

my.Bcell.freqs <- list("Bcell_YF" = my.freq.table["B cells",1:3  ],
                       "Bcell_OF" = my.freq.table["B cells",4:5  ],
                       "Bcell_YM" = my.freq.table["B cells",6:8  ],
                       "Bcell_OM" = my.freq.table["B cells",9:10 ])

my.Mph.freqs <- list("Mph_YF" = my.freq.table["Macrophages",1:3  ],
                     "Mph_OF" = my.freq.table["Macrophages",4:5  ],
                     "Mph_YM" = my.freq.table["Macrophages",6:8  ],
                     "Mph_OM" = my.freq.table["Macrophages",9:10 ])


my.Tcell.freqs <- list("Tcell_YF" = my.freq.table["T cells",1:3  ],
                       "Tcell_OF" = my.freq.table["T cells",4:5  ],
                       "Tcell_YM" = my.freq.table["T cells",6:8  ],
                       "Tcell_OM" = my.freq.table["T cells",9:10 ])

pdf(paste(Sys.Date(),"NIAv3_Peritoneal_10X_Bcell_Tcell_Mph_freq_boxplot.pdf",sep = "_"), height = 4, width = 8)
par(mfrow=c(1,3))
boxplot(my.Bcell.freqs,
        col = c("deeppink","deeppink4","deepskyblue","deepskyblue4"),
        las = 2, ylab = "B cell proportion (Ratio)",
        main = "B cell proportion (scRNAseq)",
        ylim = c(0,1))
beeswarm::beeswarm(my.Bcell.freqs, add = T, pch = 16, cex= 1)
boxplot(my.Mph.freqs,
        col = c("deeppink","deeppink4","deepskyblue","deepskyblue4"),
        las = 2, ylab = "Macrophage proportion (Ratio)",
        main = "Macrophage proportion (scRNAseq)",
        ylim = c(0,1))
beeswarm::beeswarm(my.Mph.freqs, add = T, pch = 16, cex= 1)
boxplot(my.Tcell.freqs,
        col = c("deeppink","deeppink4","deepskyblue","deepskyblue4"),
        las = 2, ylab = "T cell  proportion (Ratio)",
        main = "T cell proportion (scRNAseq)",
        ylim = c(0,1))
beeswarm::beeswarm(my.Tcell.freqs, add = T, pch = 16, cex= 1)
par(mfrow=c(1,1))
dev.off()

# save data for each cluster separately
NIA.singlets.Mph       <- subset(NIA.singlets, subset = SingleR_ImmGen %in% "Macrophages")    # Macrophages, 5709 cells
NIA.singlets.Bcells    <- subset(NIA.singlets, subset = SingleR_ImmGen %in% "B cells")        # B-cells, 4530 cells
NIA.singlets.Tcells    <- subset(NIA.singlets, subset = SingleR_ImmGen %in% "T cells")        # T-cells, 559 cells
NIA.singlets.Monocytes <- subset(NIA.singlets, subset = SingleR_ImmGen %in% "Monocytes")      # Monocytes, 107 cells

save(NIA.singlets.Mph,      file = paste(Sys.Date(),"Seurat_10x_peritoneal_Singlets_SingleR_IMMGEN_Macrophages.RData",sep = "_"))
save(NIA.singlets.Bcells,   file = paste(Sys.Date(),"Seurat_10x_peritoneal_Singlets_SingleR_IMMGEN_Bcells.RData",sep = "_"))
save(NIA.singlets.Tcells,   file = paste(Sys.Date(),"Seurat_10x_peritoneal_Singlets_SingleR_IMMGEN_Tcells.RData",sep = "_"))
save(NIA.singlets.Monocytes, file = paste(Sys.Date(),"Seurat_10x_peritoneal_Singlets_SingleR_IMMGEN_Monocytes.RData",sep = "_"))


################## Repeat with RNA-seq annotation
# How does cluster membership vary by replicate?
table(NIA.singlets@meta.data$SingleR_RNAseq, NIA.singlets@meta.data$Sample)
#               NIA_OF1 NIA_OF2 NIA_OM1 NIA_OM2 NIA_YF1 NIA_YF2 NIA_YF3 NIA_YM1 NIA_YM2 NIA_YM3
# B cells          410     868     339     673     270     642     420     367     376     227
# Erythrocytes       0       2       0       0       0       0       0       0       0       0
# Macrophages      195     383     516     580     383     926     514    1194     516     414
# Monocytes          6       7      16      20      11      44      27      23      26      10
# NK cells           2       1       0       2       6       6       4       8      10       6
# T cells           34      80     226     103      16      52      45      41      84      30

my.freq.table <- prop.table(x = table(NIA.singlets@meta.data$SingleR_RNAseq, NIA.singlets@meta.data$Sample), margin = 2)
my.freq.table <- my.freq.table[,c(5:7,1:2,8:10,3:4)]
my.freq.table.av <- apply(my.freq.table,1,mean)
my.freq.table.av.sort <- sort(my.freq.table.av, decreasing = T,index.return = T)

my.cols <- c("blue3",
             "dodgerblue",
             "cyan3",
             "chartreuse3",
             "deeppink2",
             "darkgoldenrod1",
             "coral1",
             "aquamarine2",
             "firebrick1",
             "darkviolet",
             "darkorchid1",
             "deepskyblue",
             "darkgreen",
             "yellow",
             "darkgoldenrod",
             "firebrick4",
             "dodgerblue4")

pdf(paste(Sys.Date(),"NIAv3_Peritoneal_10X_SingleR_MouseRNAseq_calls_barplot_sorted_Singlets.pdf",sep = "_"), height = 10, width = 15)
barplot(my.freq.table[my.freq.table.av.sort$ix,], las = 2,
        legend.text = rownames(my.freq.table)[my.freq.table.av.sort$ix],
        col = my.cols[my.freq.table.av.sort$ix],
        ylab = "Cell frequency in 10xGenomics single cell library")
box()
barplot(my.freq.table[my.freq.table.av.sort$ix,], las = 2,
        col = my.cols[my.freq.table.av.sort$ix],
        ylab = "Cell frequency in 10xGenomics single cell library")
box()
dev.off()

my.Bcell.freqs <- list("Bcell_YF" = my.freq.table["B cells",1:3  ],
                       "Bcell_OF" = my.freq.table["B cells",4:5  ],
                       "Bcell_YM" = my.freq.table["B cells",6:8  ],
                       "Bcell_OM" = my.freq.table["B cells",9:10 ])

my.Mph.freqs <- list("Mph_YF" = my.freq.table["Macrophages",1:3  ],
                     "Mph_OF" = my.freq.table["Macrophages",4:5  ],
                     "Mph_YM" = my.freq.table["Macrophages",6:8  ],
                     "Mph_OM" = my.freq.table["Macrophages",9:10 ])


my.Tcell.freqs <- list("Tcell_YF" = my.freq.table["T cells",1:3  ],
                       "Tcell_OF" = my.freq.table["T cells",4:5  ],
                       "Tcell_YM" = my.freq.table["T cells",6:8  ],
                       "Tcell_OM" = my.freq.table["T cells",9:10 ])

pdf(paste(Sys.Date(),"Bcell_Tcell_Mph_freq_boxplot_SingleR_RNAseq.pdf",sep = "_"), height = 4, width = 8)
par(mfrow=c(1,3))
boxplot(my.Bcell.freqs,
        col = c("deeppink","deeppink4","deepskyblue","deepskyblue4"),
        las = 2, ylab = "B cell proportion (Ratio)",
        main = "B cell proportion (scRNAseq)",
        ylim = c(0,1))
beeswarm::beeswarm(my.Bcell.freqs, add = T, pch = 16, cex= 1)
boxplot(my.Mph.freqs,
        col = c("deeppink","deeppink4","deepskyblue","deepskyblue4"),
        las = 2, ylab = "Macrophage proportion (Ratio)",
        main = "Macrophage proportion (scRNAseq)",
        ylim = c(0,1))
beeswarm::beeswarm(my.Mph.freqs, add = T, pch = 16, cex= 1)
boxplot(my.Tcell.freqs,
        col = c("deeppink","deeppink4","deepskyblue","deepskyblue4"),
        las = 2, ylab = "T cell  proportion (Ratio)",
        main = "T cell proportion (scRNAseq)",
        ylim = c(0,1))
beeswarm::beeswarm(my.Tcell.freqs, add = T, pch = 16, cex= 1)
par(mfrow=c(1,1))
dev.off()

# save data for each cluster separately
NIA.singlets.Mph.v2       <- subset(NIA.singlets, subset = SingleR_RNAseq %in% "Macrophages")    # Macrophages, 5621 cells
NIA.singlets.Bcells.v2    <- subset(NIA.singlets, subset = SingleR_RNAseq %in% "B cells")        # B-cells, 4592 cells
NIA.singlets.Tcells.v2    <- subset(NIA.singlets, subset = SingleR_RNAseq %in% "T cells")        # T-cells, 711 cells
NIA.singlets.Monocytes.v2 <- subset(NIA.singlets, subset = SingleR_RNAseq %in% "Monocytes")      # Monocytes, 190 cells

save(NIA.singlets.Mph.v2       , file = paste(Sys.Date(),"Seurat_10x_peritoneal_Singlets_SingleR_MouseRNAseq_Macrophages.RData",sep = "_"))
save(NIA.singlets.Bcells.v2    , file = paste(Sys.Date(),"Seurat_10x_peritoneal_Singlets_SingleR_MouseRNAseq_Bcells.RData",sep = "_"))
save(NIA.singlets.Tcells.v2    , file = paste(Sys.Date(),"Seurat_10x_peritoneal_Singlets_SingleR_MouseRNAseq_Tcells.RData",sep = "_"))
save(NIA.singlets.Monocytes.v2 , file = paste(Sys.Date(),"Seurat_10x_peritoneal_Singlets_SingleR_MouseRNAseq_Monocytes.RData",sep = "_"))
############################################################################################################

############################################################################################################
##### 10. Analyze SingleR annotations, compare the 2 libraries/pokedex for cell types

table(NIA.singlets@meta.data$SingleR_RNAseq, NIA.singlets@meta.data$SingleR_ImmGen)
#               B cells Basophils   DC  ILC Macrophages Mast cells Monocytes NK cells  NKT T cells  Tgd
# B cells         4498         0    9    0          75          1         6        1    0       2    0
# Erythrocytes       2         0    0    0           0          0         0        0    0       0    0
# Macrophages       24         0    1    1        5573         17         0        0    0       4    1
# Monocytes          0         1   49    0          39          0       100        0    0       1    0
# NK cells           0         0    0    7           2          0         0       35    1       0    0
# T cells            6         0    2   28          20          0         1       12    5     552   85

# label cells consistently called across both sets
sum(NIA.singlets@meta.data$SingleR_RNAseq == NIA.singlets@meta.data$SingleR_ImmGen) # 10758 (out of 11161)

# keep only cells with consistent labeling
NIA.singlets.consistent       <- subset(NIA.singlets, subset = SingleR_RNAseq == SingleR_ImmGen )   # 14375
save(NIA.singlets.consistent , file = paste(Sys.Date(),"Seurat_NIAv3_10x_peritoneal_Singlets_SingleR_ConsistentID_SingleR.RData",sep = "_"))

# get proportions
my.freq.table <- prop.table(x = table(NIA.singlets.consistent@meta.data$SingleR_RNAseq, NIA.singlets.consistent@meta.data$Sample), margin = 2)
my.freq.table <- my.freq.table[,c(5:7,1:2,8:10,3:4)]
my.freq.table.av <- apply(my.freq.table,1,mean)
my.freq.table.av.sort <- sort(my.freq.table.av, decreasing = T,index.return = T)

my.Bcell.freqs <- list("Bcell_YF" = my.freq.table["B cells",1:3  ],
                       "Bcell_OF" = my.freq.table["B cells",4:5  ],
                       "Bcell_YM" = my.freq.table["B cells",6:8  ],
                       "Bcell_OM" = my.freq.table["B cells",9:10 ])

my.Mph.freqs <- list("Mph_YF" = my.freq.table["Macrophages",1:3  ],
                     "Mph_OF" = my.freq.table["Macrophages",4:5  ],
                     "Mph_YM" = my.freq.table["Macrophages",6:8  ],
                     "Mph_OM" = my.freq.table["Macrophages",9:10 ])


my.Tcell.freqs <- list("Tcell_YF" = my.freq.table["T cells",1:3  ],
                       "Tcell_OF" = my.freq.table["T cells",4:5  ],
                       "Tcell_YM" = my.freq.table["T cells",6:8  ],
                       "Tcell_OM" = my.freq.table["T cells",9:10 ])

pdf(paste(Sys.Date(),"NIAv3_Bcell_Tcell_Mph_freq_boxplot_SingleR_CONSISTENT.pdf",sep = "_"), height = 4, width = 8)
par(mfrow=c(1,3))
boxplot(my.Bcell.freqs,
        col = c("deeppink","deeppink4","deepskyblue","deepskyblue4"),
        las = 2, ylab = "B cell proportion (Ratio)",
        main = "B cell proportion (scRNAseq)",
        ylim = c(0,1))
beeswarm::beeswarm(my.Bcell.freqs, add = T, pch = 16, cex= 1)
boxplot(my.Mph.freqs,
        col = c("deeppink","deeppink4","deepskyblue","deepskyblue4"),
        las = 2, ylab = "Macrophage proportion (Ratio)",
        main = "Macrophage proportion (scRNAseq)",
        ylim = c(0,1))
beeswarm::beeswarm(my.Mph.freqs, add = T, pch = 16, cex= 1)
boxplot(my.Tcell.freqs,
        col = c("deeppink","deeppink4","deepskyblue","deepskyblue4"),
        las = 2, ylab = "T cell  proportion (Ratio)",
        main = "T cell proportion (scRNAseq)",
        ylim = c(0,1))
beeswarm::beeswarm(my.Tcell.freqs, add = T, pch = 16, cex= 1)
par(mfrow=c(1,1))
dev.off()

NIA.singlets.Mph.v3       <- subset(NIA.singlets.consistent, subset = SingleR_ImmGen %in% "Macrophages")    # Macrophages, 5573 cells
NIA.singlets.Bcells.v3    <- subset(NIA.singlets.consistent, subset = SingleR_ImmGen %in% "B cells")        # B cells, 4498 cells

save(NIA.singlets.Mph.v3       , file = paste(Sys.Date(),"Seurat_10x_peritoneal_Singlets_SingleR_Consistent_Macrophages.RData",sep = "_"))
save(NIA.singlets.Bcells.v3    , file = paste(Sys.Date(),"Seurat_10x_peritoneal_Singlets_SingleR_Consistent_Bcells.RData",sep = "_"))

######################################
# Calculate and plot markers (on consistently called cells)
NIA.singlets.consistent  <- SetIdent(NIA.singlets.consistent , value = "SingleR_ImmGen")

table(NIA.singlets.consistent@meta.data$SingleR_ImmGen)
#     B cells Macrophages   Monocytes    NK cells     T cells 
#        4498        5573         100          35         552 

#### get marker heatmap
perit.markers <- FindAllMarkers(NIA.singlets.consistent, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10.perit.markers <- perit.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

png(paste0(Sys.Date(),"_NIAv3_SingleRconsistent_Peritoneal_aging_marker_heatmap.png"),width = 25, height = 18, units = "cm", res = 300)
DoHeatmap(NIA.singlets.consistent, features = top10.perit.markers$gene)  + theme(axis.text.y = element_text(size = 8))
dev.off()

pdf(paste0(Sys.Date(),"_NIAv3_SingleRconsistent_Peritoneal_aging_cluster_marker_Dotplot.pdf"), width = 20, height = 8)
DotPlot(NIA.singlets.consistent, features = unique(top10.perit.markers$gene), cols = c("deeppink1","deeppink4","deepskyblue1","deepskyblue4"), dot.scale = 8, split.by = "Condition") + RotatedAxis()
dev.off()

pdf(paste0(Sys.Date(),"_NIAv3_SingleRconsistent_Peritoneal_aging_KNOWN_MARKERS_Dotplot.pdf"), width = 10, height = 8)
DotPlot(NIA.singlets.consistent, features = c("Adgre1","Itgam", "Bhlhe40", "Spi1","Cd79a","Cd79b","Cd19","Cd3d","Cd3e","Cd3g","Ccr2","Ly6c2"), cols = c("deeppink1","deeppink4","deepskyblue1","deepskyblue4"), dot.scale = 8, split.by = "Condition") + RotatedAxis()
dev.off()

png(paste0(Sys.Date(),"_NIAv3_SingleRconsistent_Peritoneal_aging_KNOWN_MARKER_heatmap.png"),width = 25, height = 18, units = "cm", res = 300)
DoHeatmap(NIA.singlets.consistent,
          features = c("Cd79a","Cd79b","Cd19","Cd3d","Cd3e","Cd3g","Adgre1","Itgam", "Bhlhe40", "Spi1","Ccr2","Ly6c2"))  + theme(axis.text.y = element_text(size = 8))
dev.off()
############################################################################################################

#######################
sink(file = paste0(Sys.Date(),"_NIAv3_Perit_Seurat_Singlet_processing_session_Info.txt"))
sessionInfo()
sink()


