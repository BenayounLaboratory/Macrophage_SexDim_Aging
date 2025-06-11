setwd('/Volumes/BB_HQ_1/Immune_sex_dimorphism_Aging/Macrophages/Clean_Macrophage_Data/Peritoneal_scRNAseq/v3_set/Seurat/')
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
cellRanger.JAX_PL.1  <- load10X("/Volumes/BB_HQ_1/Immune_sex_dimorphism_Aging/Macrophages/Clean_Macrophage_Data/Peritoneal_scRNAseq/v3_set/Cell_Ranger/JAX_perit_pool1/outs/")
cellRanger.JAX_PL.2  <- load10X("/Volumes/BB_HQ_1/Immune_sex_dimorphism_Aging/Macrophages/Clean_Macrophage_Data/Peritoneal_scRNAseq/v3_set/Cell_Ranger/JAX_perit_pool2/outs/")
cellRanger.JAX_PL.3  <- load10X("/Volumes/BB_HQ_1/Immune_sex_dimorphism_Aging/Macrophages/Clean_Macrophage_Data/Peritoneal_scRNAseq/v3_set/Cell_Ranger/JAX_perit_pool3/outs/")

# estimate RNA soup from bottom 10% droplets
cellRanger.JAX_PL.1 <- autoEstCont(cellRanger.JAX_PL.1) # Estimated global rho of 0.02
cellRanger.JAX_PL.2 <- autoEstCont(cellRanger.JAX_PL.2) # Estimated global rho of 0.02
cellRanger.JAX_PL.3 <- autoEstCont(cellRanger.JAX_PL.3) # Estimated global rho of 0.01

# adjust counts based on RNA soup (and round)
out.JAX_PL.1 <- adjustCounts(cellRanger.JAX_PL.1, roundToInt=TRUE)
out.JAX_PL.2 <- adjustCounts(cellRanger.JAX_PL.2, roundToInt=TRUE)
out.JAX_PL.3 <- adjustCounts(cellRanger.JAX_PL.3, roundToInt=TRUE)

# get seurat objects
seurat.JAX_PL.1 <- CreateSeuratObject( out.JAX_PL.1 )
seurat.JAX_PL.2 <- CreateSeuratObject( out.JAX_PL.2 )
seurat.JAX_PL.3 <- CreateSeuratObject( out.JAX_PL.3 )

# Merge cleaned seurat objects
my.jax.perit <- merge(seurat.JAX_PL.1, 
                        y =  c(seurat.JAX_PL.2,
                               seurat.JAX_PL.3), 
                        add.cell.ids = c("PL_pool1",
                                         "PL_pool2",
                                         "PL_pool3"), 
                        project = "10x_JAX_peritoneal")
my.jax.perit
# An object of class Seurat 
# 32285 features across 18618 samples within 1 assay 
# Active assay: RNA (32285 features, 0 variable features)

# create batch label
Batch <- rep("NA", length(colnames(my.jax.perit@assays$RNA)))
Batch[ grep("PL_pool1"  , colnames(my.jax.perit@assays$RNA)) ]   <- "Pool_1"
Batch[ grep("PL_pool2"  , colnames(my.jax.perit@assays$RNA)) ]   <- "Pool_2"
Batch[ grep("PL_pool3"  , colnames(my.jax.perit@assays$RNA)) ]   <- "Pool_3"
Batch <- data.frame(Batch)
rownames(Batch) <- colnames(my.jax.perit@assays$RNA)

# update Seurat with metadata
my.jax.perit <- AddMetaData(object = my.jax.perit, metadata = as.vector(Batch)    , col.name = "Batch"       )

################################################################################################################################################################
#### 2. QC on mitochondrial reads and depth
# The number of genes and UMIs (nGene and nUMI) are automatically calculated for every object by Seurat.  
# The % of UMI mapping to MT-genes is a common scRNA-seq QC metric.
### QC on mitochondrial reads
my.jax.perit[["percent.mito"]] <- PercentageFeatureSet(my.jax.perit, pattern = "^mt-")

pdf(paste(Sys.Date(),"Peritoneal_cells_aging_violinPlots_QC_gene_UMI_mito.pdf", sep = "_"), height = 5, width = 10)
VlnPlot(object = my.jax.perit, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
dev.off()

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(my.jax.perit, feature1 = "nCount_RNA", feature2 = "percent.mito")
plot2 <- FeatureScatter(my.jax.perit, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

pdf(paste(Sys.Date(),"Peritoneal_cells_aging_QC_scatter.pdf", sep = "_"), height = 5, width = 10)
plot1 + plot2
dev.off()

# We filter out cells that have unique gene counts or less than 300  
my.jax.perit <- subset(my.jax.perit, subset = nFeature_RNA > 300 & percent.mito < 25 & nFeature_RNA < 6000)
my.jax.perit
# An object of class Seurat 
# 32285 features across 17212 samples within 1 assay 
# Active assay: RNA (32285 features, 0 variable features)

plot1 <- FeatureScatter(my.jax.perit, feature1 = "nCount_RNA", feature2 = "percent.mito")
plot2 <- FeatureScatter(my.jax.perit, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

pdf(paste(Sys.Date(),"Peritoneal_cells_aging_QC_scatter_post_filter.pdf", sep = "_"), height = 5, width = 10)
plot1 + plot2
dev.off()

################################################################################################################################################################
#### 3. Normalizing the data
# global-scaling normalization method ???LogNormalize??? that normalizes the gene expression measurements for each cell by the total expression, 
# multiplies this by a scale factor (10,000 by default), and log-transforms the result.
my.jax.perit <- NormalizeData(object = my.jax.perit, normalization.method = "LogNormalize",  scale.factor = 10000)

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
my.jax.perit <- CellCycleScoring(object = my.jax.perit, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# view cell cycle scores and phase assignments
head(x = my.jax.perit@meta.data)
#                                 orig.ident nCount_RNA nFeature_RNA  Batch percent.mito      S.Score    G2M.Score Phase     old.ident
#  PL_pool1_AAACCCAAGCCTCTTC-1 SeuratProject      18797         3779 Pool_1     2.761079 -0.054081575 -0.060391099    G1 SeuratProject
#  PL_pool1_AAACCCAAGTTGGCGA-1 SeuratProject       7845         2107 Pool_1     4.168260 -0.038520959 -0.004598773    G1 SeuratProject
#  PL_pool1_AAACCCACAAGTGGAC-1 SeuratProject       4310         1672 Pool_1     4.106729  0.015678614  0.091023100   G2M SeuratProject
#  PL_pool1_AAACCCACAGCTCCTT-1 SeuratProject       8562         2626 Pool_1     1.494978 -0.064303478 -0.058182994    G1 SeuratProject
#  PL_pool1_AAACCCAGTATGCTTG-1 SeuratProject       5610         1926 Pool_1     3.279857 -0.004336547  0.047163760   G2M SeuratProject
#  PL_pool1_AAACCCAGTTAACAGA-1 SeuratProject       2436         1098 Pool_1     2.216749  0.012379906 -0.098389758     S SeuratProject

################################################################################################################################################################
##### 5.Find and remove doublets 

# Normalize data for further processing
my.jax.perit <- SCTransform(object = my.jax.perit, vars.to.regress = c("nFeature_RNA", "percent.mito"))
my.jax.perit <- RunPCA(my.jax.perit)

# Determine the ‘dimensionality’ of the dataset
pdf(paste0(Sys.Date(), "_Peritoneal_Lavage_JAX_elbowplot.pdf"))
ElbowPlot(my.jax.perit, ndims = 50)
dev.off()

################################################################################
# https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html
# To give us an idea of the number of PCs needed to be included:
# We can calculate where the principal components start to elbow by taking the larger value of:
#    - The point where the principal components only contribute 5% of standard deviation
#    - The principal components cumulatively contribute 90% of the standard deviation.
#    - The point where the percent change in variation between the consecutive PCs is less than 0.1%.

# Determine percent of variation associated with each PC
pct <- my.jax.perit[["pca"]]@stdev / sum(my.jax.perit[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
co1 # 41

# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

# last point where change of % of variation is more than 0.1%.
co2 # 16

# Minimum of the two calculation
pcs <- min(co1, co2)
pcs # 16

# Based on these metrics, first 16 PCs to generate the clusters.
# We can plot the elbow plot again and overlay the information determined using our metrics:

# Create a dataframe with values
plot_df <- data.frame(pct  = pct,
                      cumu = cumu,
                      rank = 1:length(pct))

# Elbow plot to visualize
pdf(paste0(Sys.Date(), "_Peritoneal_Lavage_JAX_elbowplot_threshold_analysis.pdf"), height = 5, width= 6)
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) +
  geom_text() +
  geom_vline(xintercept = 90, color = "grey") +
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()
dev.off()
###############################################################################


my.jax.perit <- RunUMAP(my.jax.perit, dims = 1:pcs)
my.jax.perit <- FindNeighbors(my.jax.perit, dims = 1:pcs)
my.jax.perit <- FindClusters(object = my.jax.perit)


#############################################################################################################
################################  A. Doublet Finder [based on gene expression]  ################################ 
#############################################################################################################

#### need to split by original 10x sample to make sure to identify real doublets when running doublet finding

# Setup the Seurat objects as a list
perit.list <- SplitObject(my.jax.perit, split.by = "Batch")

## Assume doublet rate based on 10x information (+ 20% cells to be safe)
pred.dblt.rate <- predict(pred_dblt_lm, data.frame("cell_number" = 1.2 * unlist(lapply(perit.list, ncol))))/100
pred.dblt.rate
#    Pool_1    Pool_2    Pool_3 
# 0.0501408 0.0562848 0.0588096 

# loop over sample
for (i in 1:length(perit.list)) {
  
  ## pK Identification (no ground-truth)
  sweep.res.list_perit <- paramSweep_v3(perit.list[[i]], PCs = 1:30, sct = TRUE, num.cores	 = 4)
  sweep.stats_perit    <- summarizeSweep(sweep.res.list_perit , GT = FALSE)
  bcmvn_perit          <- find.pK(sweep.stats_perit)
  
  # need some R gymnastics since the Pk is stored as a factor for some reason
  # to get the pK number, need to first convert to character and THEN to numeric
  # numeric first yeild row number
  pk.perit <- as.numeric(as.character(bcmvn_perit[as.numeric(bcmvn_perit$pK[bcmvn_perit$BCmetric == max(bcmvn_perit$BCmetric)]),"pK"]))
  
  ## Homotypic Doublet Proportion Estimate
  homotypic.prop <- modelHomotypic(perit.list[[i]]@meta.data$seurat_clusters)     ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
  nExp_poi       <- round(pred.dblt.rate[i]*nrow(perit.list[[i]]@meta.data))        ## Assuming doublet formation rate based on 10x manual
  
  ## Run DoubletFinder with varying classification stringencies
  perit.list[[i]] <- doubletFinder_v3(perit.list[[i]], PCs = 1:30, pN = 0.25, pK = pk.perit, nExp = nExp_poi,     reuse.pANN = FALSE, sct = T)
  
  # rename column to enable subsetting
  colnames(perit.list[[i]]@meta.data)[grep("DF.classifications_0.25",colnames(perit.list[[i]]@meta.data))] <- "DoubletFinder"
  
}

# plot UMAP plots for each pool
for (i in 1:length(perit.list)) {
  
  # get classification name
  pdf(paste(Sys.Date(),"JAX_Peritoneal",names(perit.list)[i],"Doublet_Finder_UMAP.pdf", sep = "_"), height = 5, width = 5)
  print(DimPlot(perit.list[[i]], reduction = "umap", group.by = "DoubletFinder"))
  dev.off()
}

table(perit.list[[1]]@meta.data$DoubletFinder)
# Doublet Singlet 
#     262    4961 

table(perit.list[[2]]@meta.data$DoubletFinder)
# Doublet Singlet 
# 330    5531 

table(perit.list[[3]]@meta.data$DoubletFinder)
# Doublet Singlet 
# 360    5763 

save(perit.list, file = paste0(Sys.Date(),"_JAX_Peritoneal_cells_aging_Seurat_object_DoubletFinderList.RData"))


#############################################################################################################
##############################  A. HTO analysis [based on TotalSeq antibodies]  ############################# 
#############################################################################################################

#################################################################
############### a. read in cite-seq tools results ###############
#################################################################

################################ read HTO libraries 1
HTO_p1 <- Read10X('/Volumes/BB_HQ_1/Immune_sex_dimorphism_Aging/Macrophages/Clean_Macrophage_Data/Peritoneal_scRNAseq/v3_set/HTO_analysis/perit_p_HTO_1_HTO_parsed/umi_count', gene.column=1)
HTO_p1 <- HTO_p1[1:4,] # remove unmapped

# need to clean up cells with no HTOs
my.zero <- which(colSums(as.matrix(HTO_p1)) < 100)
HTO_p1 <- HTO_p1[,-my.zero] # remove cells with no hashtag
dim(HTO_p1) # 4 6369

################################ read HTO libraries 2
HTO_p2 <- Read10X('/Volumes/BB_HQ_1/Immune_sex_dimorphism_Aging/Macrophages/Clean_Macrophage_Data/Peritoneal_scRNAseq/v3_set/HTO_analysis/perit_p_HTO_2_HTO_parsed/umi_count', gene.column=1)
HTO_p2 <- HTO_p2[1:4,] # remove unmapped

# need to clean up cells with no HTOs
my.zero <- which(colSums(as.matrix(HTO_p2)) < 100)
HTO_p2 <- HTO_p2[,-my.zero] # remove cells with no hashtag
dim(HTO_p2) #  4 7076

################################ read HTO libraries 3
HTO_p3 <- Read10X('/Volumes/BB_HQ_1/Immune_sex_dimorphism_Aging/Macrophages/Clean_Macrophage_Data/Peritoneal_scRNAseq/v3_set/HTO_analysis/perit_p_HTO_3_HTO_parsed/umi_count', gene.column=1)
HTO_p3 <- HTO_p3[1:4,] # remove unmapped

# need to clean up cells with no HTOs
my.zero <- which(colSums(as.matrix(HTO_p3)) < 100)
HTO_p3 <- HTO_p3[,-my.zero] # remove cells with no hashtag
dim(HTO_p3) # 4 7322


################################################################
###############  b. Filter and demultiplex data  ###############
################################################################

# extract barcode information
perit.list[[1]] <- AddMetaData(object = perit.list[[1]], metadata = as.vector(unlist(lapply(strsplit(rownames(perit.list[[1]]@meta.data),"_|-"),'[[',3))), col.name = "barcode")
perit.list[[2]] <- AddMetaData(object = perit.list[[2]], metadata = as.vector(unlist(lapply(strsplit(rownames(perit.list[[2]]@meta.data),"_|-"),'[[',3))), col.name = "barcode")
perit.list[[3]] <- AddMetaData(object = perit.list[[3]], metadata = as.vector(unlist(lapply(strsplit(rownames(perit.list[[3]]@meta.data),"_|-"),'[[',3))), col.name = "barcode")


################################ process Pool #1 ################################
# Select cell barcodes detected by both RNA and HTO
joint.bcs.1 <- intersect(colnames(HTO_p1), perit.list[[1]]@meta.data$barcode) # 5090

# Subset RNA Suerat object by joint cell barcodes
JAX.p1.umis <- perit.list[[1]][, perit.list[[1]]@meta.data$barcode  %in% joint.bcs.1]

# Subset HTO counts by joint cell barcodes (correct order from subsetted seurat above)
JAX.p1.htos <- as.matrix(HTO_p1[, JAX.p1.umis@meta.data$barcode])

# make names that will correspond to Seurat cell names
colnames(JAX.p1.htos) <- rownames(JAX.p1.umis@meta.data)

# Confirm that the HTO have the correct names
rownames(JAX.p1.htos) # "HTO1_F_4m_1-ACCCACCAGTAAGAC"  "HTO2_M_4m_1-GGTCGAGAGCATTCA"  "HTO3_F_20m_1-CTTGCCGCATGTCAT" "HTO4_M_20m_1-AAAGCATTCTTCACG"

# create a Seurat object copy
JAX.p1.hashtag <- JAX.p1.umis

# Normalize RNA data with log normalization
JAX.p1.hashtag <- NormalizeData(JAX.p1.hashtag, normalization.method = "LogNormalize",  scale.factor = 10000)

# Find and scale variable features
JAX.p1.hashtag <- FindVariableFeatures(JAX.p1.hashtag, selection.method = "mean.var.plot")
JAX.p1.hashtag <- ScaleData(JAX.p1.hashtag, features = VariableFeatures(JAX.p1.hashtag))

# Add HTO data as a new assay independent from RNA
JAX.p1.hashtag[["HTO"]] <- CreateAssayObject(counts = JAX.p1.htos)

# Normalize HTO data, here we use centered log-ratio (CLR) transformation
JAX.p1.hashtag <- NormalizeData(JAX.p1.hashtag, assay = "HTO", normalization.method = "CLR")

# Demultiplex cells based on HTO enrichment
JAX.p1.hashtag <- MULTIseqDemux(JAX.p1.hashtag, assay = "HTO") # https://github.com/satijalab/seurat/issues/2549

# Global classification results
table(JAX.p1.hashtag$MULTI_ID)
# Doublet  HTO1-F-4m-1-ACCCACCAGTAAGAC  HTO2-M-4m-1-GGTCGAGAGCATTCA HTO3-F-20m-1-CTTGCCGCATGTCAT 
#     214                         1127                         1096                         1298 
# HTO4-M-20m-1-AAAGCATTCTTCACG                     Negative 
#                         1244                          111 

# Visualize enrichment for selected HTOs with ridge plots, and group cells based on the max HTO signal
Idents(JAX.p1.hashtag) <- "MULTI_ID"

pdf(paste0(Sys.Date(),"_JAX_Pool1_HTO_id_ridge_plot.pdf"), height = 10, width = 15)
RidgePlot(JAX.p1.hashtag, assay = "HTO", features = rownames(JAX.p1.hashtag[["HTO"]]), ncol = 2)
dev.off()

pdf(paste0(Sys.Date(),"_JAX_Pool1_HTO1_2_feature_Scatter_plot.pdf"), height = 5, width = 10)
FeatureScatter(JAX.p1.hashtag, feature1 = "HTO3-F-20m-1-CTTGCCGCATGTCAT", feature2 = "HTO4-M-20m-1-AAAGCATTCTTCACG")
dev.off()

# Compare number of UMIs for singlets, doublets and negative cells
pdf(paste0(Sys.Date(),"_JAX_Pool1_nCount_RNA_violin_based_on_HTO_Demux.pdf"), height = 10, width = 15)
VlnPlot(JAX.p1.hashtag, features = "nCount_RNA", pt.size = 0.1, log = TRUE)
dev.off()

# Generate a two dimensional UMAP embedding for HTOs (grouping cells by singlets and doublets for simplicity)
# First, we will remove HTO-negative cells from the object (since they can't be attributed anyway)
JAX.p1.hashtag.subset <- subset(JAX.p1.hashtag, idents = "Negative", invert = TRUE)

# Calculate a UMAP embedding of the HTO data
DefaultAssay(JAX.p1.hashtag.subset) <- "HTO"
JAX.p1.hashtag.subset <- ScaleData(JAX.p1.hashtag.subset, features = rownames(JAX.p1.hashtag.subset), verbose = FALSE)
JAX.p1.hashtag.subset <- RunPCA(JAX.p1.hashtag.subset, features = rownames(JAX.p1.hashtag.subset), approx = FALSE)
JAX.p1.hashtag.subset <- RunUMAP(JAX.p1.hashtag.subset, dims = 1:4, perplexity = 100, check_duplicates = FALSE)

pdf(paste0(Sys.Date(),"_JAX_Pool1_UMAP_hashtag.pdf"), height = 6, width = 10)
DimPlot(JAX.p1.hashtag.subset)
dev.off()

# create metadata column to plot HTO heatmap
JAX.p1.hashtag@meta.data$HTO.global <- rep("Singlet",dim(JAX.p1.hashtag@meta.data)[1])
JAX.p1.hashtag@meta.data$HTO.global[JAX.p1.hashtag@meta.data$MULTI_ID %in% "Doublet"] <- "Doublet"
JAX.p1.hashtag@meta.data$HTO.global[JAX.p1.hashtag@meta.data$MULTI_ID %in% "Negative"] <- "Negative"

# Create an HTO heatmap, based on Figure 1C in the Cell Hashing paper.
pdf(paste0(Sys.Date(),"_JAX_Pool1_HTO_heatmap.pdf"), height = 5, width = 9)
HTOHeatmap(JAX.p1.hashtag, assay = "HTO",
           classification = "MULTI_ID",
           global.classification = "HTO.global")
dev.off()

save(JAX.p1.hashtag, file = paste0(Sys.Date(),"_JAX_Pool1_post_HTO_Demux.RData"))

################################ process Pool #2 ################################
# Select cell barcodes detected by both RNA and HTO
joint.bcs.2 <- intersect(colnames(HTO_p2), perit.list[[2]]@meta.data$barcode) # 5746

# Subset RNA Suerat object by joint cell barcodes
JAX.p2.umis <- perit.list[[2]][, perit.list[[2]]@meta.data$barcode  %in% joint.bcs.2]

# Subset HTO counts by joint cell barcodes (correct order from subsetted seurat above)
JAX.p2.htos <- as.matrix(HTO_p2[, JAX.p2.umis@meta.data$barcode])

# make names that will correspond to Seurat cell names
colnames(JAX.p2.htos) <- rownames(JAX.p2.umis@meta.data)

# Confirm that the HTO have the correct names
rownames(JAX.p2.htos) # "HTO6_F_20m_2-TATGCTGCCACGGTA"  "HTO7_M_4m_2-GAGTCTGCCAGTATC"   "HTO8_F_4m_2-TATAGAACGCCAGGC"   "HTO10_M_20m_3-CCGATTGTAACAGAC"

# create a Seurat object copy
JAX.p2.hashtag <- JAX.p2.umis

# Normalize RNA data with log normalization
JAX.p2.hashtag <- NormalizeData(JAX.p2.hashtag, normalization.method = "LogNormalize",  scale.factor = 10000)

# Find and scale variable features
JAX.p2.hashtag <- FindVariableFeatures(JAX.p2.hashtag, selection.method = "mean.var.plot")
JAX.p2.hashtag <- ScaleData(JAX.p2.hashtag, features = VariableFeatures(JAX.p2.hashtag))

# Add HTO data as a new assay independent from RNA
JAX.p2.hashtag[["HTO"]] <- CreateAssayObject(counts = JAX.p2.htos)

# Normalize HTO data, here we use centered log-ratio (CLR) transformation
JAX.p2.hashtag <- NormalizeData(JAX.p2.hashtag, assay = "HTO", normalization.method = "CLR")

# Demultiplex cells based on HTO enrichment
JAX.p2.hashtag <- MULTIseqDemux(JAX.p2.hashtag, assay = "HTO") # https://github.com/satijalab/seurat/issues/2549

# Global classification results
table(JAX.p2.hashtag$MULTI_ID)
#    Doublet HTO10-M-20m-3-CCGATTGTAACAGAC  HTO6-F-20m-2-TATGCTGCCACGGTA   HTO7-M-4m-2-GAGTCTGCCAGTATC   HTO8-F-4m-2-TATAGAACGCCAGGC 
#        236                          1301                          1150                          1479                          1486 
#   Negative 
#         94 

# Visualize enrichment for selected HTOs with ridge plots, and group cells based on the max HTO signal
Idents(JAX.p2.hashtag) <- "MULTI_ID"

pdf(paste0(Sys.Date(),"_JAX_Pool2_HTO_id_ridge_plot.pdf"), height = 10, width = 15)
RidgePlot(JAX.p2.hashtag, assay = "HTO", features = rownames(JAX.p2.hashtag[["HTO"]]), ncol = 2)
dev.off()

pdf(paste0(Sys.Date(),"_JAX_Pool2_HTO10_6_feature_Scatter_plot.pdf"), height = 5, width = 10)
FeatureScatter(JAX.p2.hashtag, feature1 = "HTO10-M-20m-3-CCGATTGTAACAGAC", feature2 = "HTO6-F-20m-2-TATGCTGCCACGGTA")
dev.off()

# Compare number of UMIs for singlets, doublets and negative cells
pdf(paste0(Sys.Date(),"_JAX_Pool2_nCount_RNA_violin_based_on_HTO_Demux.pdf"), height = 10, width = 15)
VlnPlot(JAX.p2.hashtag, features = "nCount_RNA", pt.size = 0.1, log = TRUE)
dev.off()

# Generate a two dimensional UMAP embedding for HTOs (grouping cells by singlets and doublets for simplicity)
# First, we will remove HTO-negative cells from the object (since they can't be attributed anyway)
JAX.p2.hashtag.subset <- subset(JAX.p2.hashtag, idents = "Negative", invert = TRUE)

# Calculate a UMAP embedding of the HTO data
DefaultAssay(JAX.p2.hashtag.subset) <- "HTO"
JAX.p2.hashtag.subset <- ScaleData(JAX.p2.hashtag.subset, features = rownames(JAX.p2.hashtag.subset), verbose = FALSE)
JAX.p2.hashtag.subset <- RunPCA(JAX.p2.hashtag.subset, features = rownames(JAX.p2.hashtag.subset), approx = FALSE)
JAX.p2.hashtag.subset <- RunUMAP(JAX.p2.hashtag.subset, dims = 1:4, perplexity = 100, check_duplicates = FALSE)

pdf(paste0(Sys.Date(),"_JAX_Pool2_UMAP_hashtag.pdf"), height = 6, width = 10)
DimPlot(JAX.p2.hashtag.subset)
dev.off()

# create metadata column to plot HTO heatmap
JAX.p2.hashtag@meta.data$HTO.global <- rep("Singlet",dim(JAX.p2.hashtag@meta.data)[1])
JAX.p2.hashtag@meta.data$HTO.global[JAX.p2.hashtag@meta.data$MULTI_ID %in% "Doublet"] <- "Doublet"
JAX.p2.hashtag@meta.data$HTO.global[JAX.p2.hashtag@meta.data$MULTI_ID %in% "Negative"] <- "Negative"

# Create an HTO heatmap, based on Figure 1C in the Cell Hashing paper.
pdf(paste0(Sys.Date(),"_JAX_Pool2_HTO_heatmap.pdf"), height = 5, width = 9)
HTOHeatmap(JAX.p2.hashtag, assay = "HTO",
           classification = "MULTI_ID",
           global.classification = "HTO.global")
dev.off()

save(JAX.p2.hashtag, file = paste0(Sys.Date(),"_JAX_Pool2_post_HTO_Demux.RData"))


################################ process Pool #3 ################################
# Select cell barcodes detected by both RNA and HTO
joint.bcs.3 <- intersect(colnames(HTO_p3), perit.list[[3]]@meta.data$barcode) # 5892

# Subset RNA Suerat object by joint cell barcodes
JAX.p3.umis <- perit.list[[3]][, perit.list[[3]]@meta.data$barcode  %in% joint.bcs.3] 

# Subset HTO counts by joint cell barcodes (correct order from subsetted seurat above)
JAX.p3.htos <- as.matrix(HTO_p3[, JAX.p3.umis@meta.data$barcode])

# make names that will correspond to Seurat cell names
colnames(JAX.p3.htos) <- rownames(JAX.p3.umis@meta.data)

# Confirm that the HTO have the correct names
rownames(JAX.p3.htos) # "HTO5_M_20m_2-CTTTGTCTTTGTGAG" "HTO9_F_4m_3-TGCCTATGAAACAAG"  "HTO10_M_4m_3-CCGATTGTAACAGAC" "HTO7_F_20m_3-GAGTCTGCCAGTATC"

# create a Seurat object copy
JAX.p3.hashtag <- JAX.p3.umis

# Normalize RNA data with log normalization
JAX.p3.hashtag <- NormalizeData(JAX.p3.hashtag, normalization.method = "LogNormalize",  scale.factor = 10000)

# Find and scale variable features
JAX.p3.hashtag <- FindVariableFeatures(JAX.p3.hashtag, selection.method = "mean.var.plot")
JAX.p3.hashtag <- ScaleData(JAX.p3.hashtag, features = VariableFeatures(JAX.p3.hashtag))

# Add HTO data as a new assay independent from RNA
JAX.p3.hashtag[["HTO"]] <- CreateAssayObject(counts = JAX.p3.htos)

# Normalize HTO data, here we use centered log-ratio (CLR) transformation
JAX.p3.hashtag <- NormalizeData(JAX.p3.hashtag, assay = "HTO", normalization.method = "CLR")

# Demultiplex cells based on HTO enrichment
JAX.p3.hashtag <- MULTIseqDemux(JAX.p3.hashtag, assay = "HTO") # https://github.com/satijalab/seurat/issues/2549

# Global classification results
table(JAX.p3.hashtag$MULTI_ID)
#    Doublet HTO10-M-4m-3-CCGATTGTAACAGAC HTO5-M-20m-2-CTTTGTCTTTGTGAG HTO7-F-20m-3-GAGTCTGCCAGTATC  HTO9-F-4m-3-TGCCTATGAAACAAG 
#        345                         1428                         1228                         1379                         1427 
#   Negative 
#         85 

# Visualize enrichment for selected HTOs with ridge plots, and group cells based on the max HTO signal
Idents(JAX.p3.hashtag) <- "MULTI_ID"

pdf(paste0(Sys.Date(),"_JAX_Pool3_HTO_id_ridge_plot.pdf"), height = 10, width = 15)
RidgePlot(JAX.p3.hashtag, assay = "HTO", features = rownames(JAX.p3.hashtag[["HTO"]]), ncol = 2)
dev.off()

pdf(paste0(Sys.Date(),"_JAX_Pool3_HTO10_5_feature_Scatter_plot.pdf"), height = 5, width = 10)
FeatureScatter(JAX.p3.hashtag, feature1 = "HTO10-M-4m-3-CCGATTGTAACAGAC", feature2 = "HTO5-M-20m-2-CTTTGTCTTTGTGAG")
dev.off()

# Compare number of UMIs for singlets, doublets and negative cells
pdf(paste0(Sys.Date(),"_JAX_Pool3_nCount_RNA_violin_based_on_HTO_Demux.pdf"), height = 10, width = 15)
VlnPlot(JAX.p3.hashtag, features = "nCount_RNA", pt.size = 0.1, log = TRUE)
dev.off()

# Generate a two dimensional UMAP embedding for HTOs (grouping cells by singlets and doublets for simplicity)
# First, we will remove HTO-negative cells from the object (since they can't be attributed anyway)
JAX.p3.hashtag.subset <- subset(JAX.p3.hashtag, idents = "Negative", invert = TRUE)

# Calculate a UMAP embedding of the HTO data
DefaultAssay(JAX.p3.hashtag.subset) <- "HTO"
JAX.p3.hashtag.subset <- ScaleData(JAX.p3.hashtag.subset, features = rownames(JAX.p3.hashtag.subset), verbose = FALSE)
JAX.p3.hashtag.subset <- RunPCA(JAX.p3.hashtag.subset, features = rownames(JAX.p3.hashtag.subset), approx = FALSE)
JAX.p3.hashtag.subset <- RunUMAP(JAX.p3.hashtag.subset, dims = 1:4, perplexity = 100, check_duplicates = FALSE)

pdf(paste0(Sys.Date(),"_JAX_Pool3_UMAP_hashtag.pdf"), height = 6, width = 10)
DimPlot(JAX.p3.hashtag.subset)
dev.off()

# create metadata column to plot HTO heatmap
JAX.p3.hashtag@meta.data$HTO.global <- rep("Singlet",dim(JAX.p3.hashtag@meta.data)[1])
JAX.p3.hashtag@meta.data$HTO.global[JAX.p3.hashtag@meta.data$MULTI_ID %in% "Doublet"] <- "Doublet"
JAX.p3.hashtag@meta.data$HTO.global[JAX.p3.hashtag@meta.data$MULTI_ID %in% "Negative"] <- "Negative"

# Create an HTO heatmap, based on Figure 1C in the Cell Hashing paper.
pdf(paste0(Sys.Date(),"_JAX_Pool3_HTO_heatmap.pdf"), height = 5, width = 9)
HTOHeatmap(JAX.p3.hashtag, assay = "HTO",
           classification = "MULTI_ID",
           global.classification = "HTO.global")
dev.off()

save(JAX.p3.hashtag, file = paste0(Sys.Date(),"_JAX_Pool3_post_HTO_Demux.RData"))


#############################################################################
###############  c. Merge split objects and compare doublets  ###############
#############################################################################

# combined the 3 objects
JAX.combined.calls <- merge(JAX.p1.hashtag,
                            y = c(JAX.p2.hashtag,
                                  JAX.p3.hashtag),
                            project = "10x_JAX_Peritoneal_cells")
JAX.combined.calls
# An object of class Seurat 
# 49302 features across 16728 samples within 3 assays 
# Active assay: SCT (17005 features, 0 variable features)
#  2 other assays present: RNA, HTO

# Compare HTO and scds calls
table(JAX.combined.calls@meta.data$DoubletFinder,JAX.combined.calls@meta.data$HTO.global)
#            Doublet Negative Singlet
#    Doublet     217       29     580
#    Singlet     575      268   15055

# Keep cells marked as singlets from both methods
JAX.singlets <- JAX.combined.calls[, bitAnd(JAX.combined.calls@meta.data$DoubletFinder == "Singlet", JAX.combined.calls@meta.data$HTO.global == "Singlet")>0]
JAX.singlets
# An object of class Seurat 
# 49302 features across 14924 samples within 3 assays 
# Active assay: SCT (17005 features, 0 variable features)
#  2 other assays present: RNA, HTO

# table(JAX.singlets@meta.data$DoubletFinder,JAX.singlets@meta.data$HTO.global) ### ALL singlets, sanity check
save(JAX.singlets, file = paste0(Sys.Date(),"_JAX_Combined_Singlets_ONLY_Seurat_object.RData"))


################################################################################################################################################################
##### 6.Scaling the data and removing unwanted sources of variation

# Run SCT normalization, regress potential problems
JAX.singlets <- SCTransform(object = JAX.singlets, vars.to.regress = c("nFeature_RNA", "percent.mito", "Batch","S.Score", "G2M.Score"))

# find variable genes
JAX.singlets <- FindVariableFeatures(JAX.singlets, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(JAX.singlets), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(JAX.singlets)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

pdf(paste(Sys.Date(),"JAX_Peritoneal_cells_aging_Singlets_Variable genes.pdf", sep = "_"), height = 5, width = 10)
plot1 + plot2
dev.off()


################################################################################################################################################################
##### 7. Update metadata information based on HTO and output information on cell cycling

table(JAX.singlets$MULTI_classification)
#  HTO1-F-4m-1-ACCCACCAGTAAGAC HTO10-M-20m-3-CCGATTGTAACAGAC  HTO10-M-4m-3-CCGATTGTAACAGAC   HTO2-M-4m-1-GGTCGAGAGCATTCA  HTO3-F-20m-1-CTTGCCGCATGTCAT 
#                         1081                          1265                          1373                          1052                          1256 
# HTO4-M-20m-1-AAAGCATTCTTCACG  HTO5-M-20m-2-CTTTGTCTTTGTGAG  HTO6-F-20m-2-TATGCTGCCACGGTA  HTO7-F-20m-3-GAGTCTGCCAGTATC   HTO7-M-4m-2-GAGTCTGCCAGTATC 
#                         1201                          1205                          1110                          1352                          1393 
#  HTO8-F-4m-2-TATAGAACGCCAGGC   HTO9-F-4m-3-TGCCTATGAAACAAG 
#                         1422                          1345 

#######################
# create condition label
my.YF <- grep("F-4m" , JAX.singlets$MULTI_classification)
my.OF <- grep("F-20m", JAX.singlets$MULTI_classification)
my.YM <- grep("M-4m" , JAX.singlets$MULTI_classification)
my.OM <- grep("M-20m", JAX.singlets$MULTI_classification)

Condition <- rep("NA", length(colnames(JAX.singlets@assays$RNA)))
Condition[ my.YF  ]   <- "YF"
Condition[ my.OF  ]   <- "OF"
Condition[ my.YM  ]   <- "YM"
Condition[ my.OM  ]   <- "OM"

Condition <- data.frame(Condition)
rownames(Condition) <- colnames(JAX.singlets@assays$RNA)

JAX.singlets <- AddMetaData(object = JAX.singlets, metadata = as.vector(Condition), col.name = "Condition")


#######################
# create sample label
my.YF1 <- grep("F-4m-1" , JAX.singlets$MULTI_classification)
my.YF2 <- grep("F-4m-2" , JAX.singlets$MULTI_classification)
my.YF3 <- grep("F-4m-3" , JAX.singlets$MULTI_classification)
my.OF1 <- grep("F-20m-1", JAX.singlets$MULTI_classification)
my.OF2 <- grep("F-20m-2", JAX.singlets$MULTI_classification)
my.OF3 <- grep("F-20m-3", JAX.singlets$MULTI_classification)
my.YM1 <- grep("M-4m-1" , JAX.singlets$MULTI_classification)
my.YM2 <- grep("M-4m-2" , JAX.singlets$MULTI_classification)
my.YM3 <- grep("M-4m-3" , JAX.singlets$MULTI_classification)
my.OM1 <- grep("M-20m-1", JAX.singlets$MULTI_classification)
my.OM2 <- grep("M-20m-2", JAX.singlets$MULTI_classification)
my.OM3 <- grep("M-20m-3", JAX.singlets$MULTI_classification)

Sample <- rep("NA", length(colnames(JAX.singlets@assays$RNA)))
Sample[ my.YF1  ]   <- "JAX_YF1"
Sample[ my.YF2  ]   <- "JAX_YF2"
Sample[ my.YF3  ]   <- "JAX_YF3"
Sample[ my.OF1  ]   <- "JAX_OF1"
Sample[ my.OF2  ]   <- "JAX_OF2"
Sample[ my.OF3  ]   <- "JAX_OF3"
Sample[ my.YM1  ]   <- "JAX_YM1"
Sample[ my.YM2  ]   <- "JAX_YM2"
Sample[ my.YM3  ]   <- "JAX_YM3"
Sample[ my.OM1  ]   <- "JAX_OM1"
Sample[ my.OM2  ]   <- "JAX_OM2"
Sample[ my.OM3  ]   <- "JAX_OM3"

Sample <- data.frame(Sample)
rownames(Sample) <- colnames(JAX.singlets@assays$RNA)

JAX.singlets <- AddMetaData(object = JAX.singlets, metadata = as.vector(Sample), col.name = "Sample_ID")

#### check
head(JAX.singlets@meta.data)


#### Parse Cell Cyle information

# write predictions to file
write.table(JAX.singlets@meta.data, file = paste(Sys.Date(),"Peritoneal_cells_aging_v3JAX_CellCycle_predictions.txt", sep = "_"), sep = "\t", quote = F)


pdf(paste0(Sys.Date(), "Peritoneal_aging_v3JAX_cell_cycle_pie_charts.pdf"), height = 13, width = 8)
par(mfrow=c(2,2))
par(oma=c(0.1,0.1,0.1,0.1))
pie(table(JAX.singlets@meta.data$Phase[JAX.singlets@meta.data$Condition == "YF"]) , main = "YF")
pie(table(JAX.singlets@meta.data$Phase[JAX.singlets@meta.data$Condition == "YM"]) , main = "YM")
pie(table(JAX.singlets@meta.data$Phase[JAX.singlets@meta.data$Condition == "OF"]) , main = "OF")
pie(table(JAX.singlets@meta.data$Phase[JAX.singlets@meta.data$Condition == "OM"]) , main = "OM")
par(mfrow=c(1,1))
dev.off()


################################################################################################################################################################
##### 8. Cluster the cells

# rerun dimensionality reduction on merged/clean object
JAX.singlets <- RunPCA(JAX.singlets)
JAX.singlets <- RunUMAP(JAX.singlets, dims = 1:16)

# look at samples
pdf(paste(Sys.Date(),"Peritoneal_cells_aging_v3JAX_Singlets_UMAP_by_Group.pdf", sep = "_"), height = 5, width = 5)
DimPlot(JAX.singlets, reduction = "umap", group.by = "Batch")
DimPlot(JAX.singlets, reduction = "umap", group.by = "Condition")
DimPlot(JAX.singlets, reduction = "umap", group.by = "Sample_ID")
dev.off()

# perform clustering again now that only singlets remain
JAX.singlets <- FindNeighbors(JAX.singlets, 
                              dims = 1:16)
JAX.singlets <- FindClusters(object = JAX.singlets,
                             reduction = "pca",
                             k.param = 100,
                             dims = 1:16,
                             resolution = 0.3,
                             n.start = 100)
# Number of communities: 10

pdf(paste(Sys.Date(),"JAX_Peritoneal_cells_aging_UMAP_clusters_0.3res_Singlets.pdf", sep = "_"), height = 5, width = 6.5)
DimPlot(JAX.singlets, label = TRUE)
dev.off()

write.table(JAX.singlets@meta.data, file = paste(Sys.Date(),"JAX_Peritoneal_cells_aging_Seurat_PCA_SNN_clusters_and_Metadata_Singlets.txt", sep = "_"), sep = "\t", quote = F, col.names = F, row.names = T)

pdf(paste(Sys.Date(),"JAX_Peritoneal_cells_aging_Feature_plot_Singlets_UMI_counts.pdf", sep = "_"),  height = 5, width = 6.5)
FeaturePlot(object = JAX.singlets, features = c("nFeature_RNA")  )
dev.off()

################################################################################################################################################################
##### 9. Visualization of potential marker gene expression

### key markers
pdf(paste(Sys.Date(),"JAX_Peritoneal_cells_aging_Feature_plot_Singlets.pdf", sep = "_"), width = 15, height = 15)
FeaturePlot(object = JAX.singlets, features = c("Cd79a","Adgre1","Itgam","Bhlhe40","Cd14","Cd3e","Spi1","Cd209a")  )
FeaturePlot(object = JAX.singlets, features = c("Cd86", "Cd68", "Cd80","Tlr2", "Tlr4")  ) ### M1 "Cx3cr1","Nos2","Il1r", "Il1b"
FeaturePlot(object = JAX.singlets, features = c("Arg1", "Cd163")  ) # m2 "Cd206","Fizz1" Ym1/2
dev.off()

pdf(paste(Sys.Date(),"JAX_Peritoneal_cells_aging_Xist_ridge_plot_Singlets.pdf", sep = "_"))
RidgePlot(object = JAX.singlets, features = "Xist", ncol = 1, group.by = "Condition")
dev.off()

pdf(paste(Sys.Date(),"JAX_Peritoneal_cells_aging_Y_genes_ridge_plot_Singlets.pdf", sep = "_"), height = 5, width = 10)
RidgePlot(object = JAX.singlets, features = c("Ddx3y", "Eif2s3y"), ncol = 2, group.by = "Condition")
dev.off()

save(JAX.singlets, file = paste(Sys.Date(),"JAX_Peritoneal_10X_Singlets.Seurat.object.RData",sep = "_"))


################################################################################################################################################################
##### 10. Identify cell types using SingleR

################ run singleR and save object
my.singler = CreateSinglerObject(counts = JAX.singlets@assays$SCT@counts,
                                 annot = JAX.singlets@meta.data$Sample_ID,
                                 min.genes = 250, technology = "10X",
                                 species = "Mouse", project.name = "10X_Peritoneal")

my.singler$seurat               = JAX.singlets # (optional)
my.singler$meta.data$orig.ident = JAX.singlets@meta.data$orig.ident # the original identities, if not supplied in 'annot'

## if using Seurat v3.0 and over use:
my.singler$meta.data$xy       =    JAX.singlets@reductions$umap@cell.embeddings # the Umap coordinates
my.singler$meta.data$clusters =    JAX.singlets@active.ident                    # the Seurat clusters (if 'clusters' not provided)
my.singler$meta.data$Sample   =    JAX.singlets@meta.data$Sample                # the Seurat clusters (if 'clusters' not provided)

save(my.singler, file = paste(Sys.Date(),"JAX_Peritoneal_10X_SingleR_object_SCT.RData",sep = "_"))

#######################
### Main cell types
pdf(paste(Sys.Date(),"JAX_Peritoneal_10X_SingleR_object.main.cell.types.pdf",sep = "_"), height = 5, width = 7)
SingleR.DrawHeatmap(my.singler$singler[[1]]$SingleR.single.main, top.n = 50, clusters = my.singler$meta.data$Sample)
dev.off()

### All cell types
pdf(paste(Sys.Date(),"JAX_Peritoneal_10X_SingleR_object.all.cell.types.pdf",sep = "_"), height = 5, width = 7)
SingleR.DrawHeatmap(my.singler$singler[[1]]$SingleR.single, top.n = 50, clusters = my.singler$meta.data$Sample)
dev.off()

### Main cell types
pdf(paste(Sys.Date(),"JAX_Peritoneal_10X_SingleR_object.main.cell.types_MouseRNAseq.pdf",sep = "_"), height = 5, width = 7)
SingleR.DrawHeatmap(my.singler$singler[[2]]$SingleR.single.main, top.n = Inf,clusters = my.singler$meta.data$Sample)
dev.off()

### All cell types
pdf(paste(Sys.Date(),"JAX_Peritoneal_10X_SingleR_object.all.cell.types_MouseRNAseq.pdf",sep = "_"), height = 5, width = 7)
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

pdf(paste(Sys.Date(),"JAX_Peritoneal_10X_SingleR_UMAP_annotated_IMMGEN.pdf",sep = "_"), height = 10, width = 12)
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

pdf(paste(Sys.Date(),"JAX_Peritoneal_10X_SingleR_UMAP_annotated_MouseRNAseq.pdf",sep = "_"), height = 10, width = 12)
out$p
dev.off()

write.table(my.singler$singler[[1]]$SingleR.single.main$labels,
            file = paste(Sys.Date(),"JAX_Peritoneal_10X_SingleR_annotated_cells_Immgen_Main.txt",sep = "_"), sep = "\t", quote = F)
write.table(my.singler$singler[[1]]$SingleR.single$labels,
            file = paste(Sys.Date(),"JAX_Peritoneal_10X_SingleR_annotated_cells_Immgen_All.txt",sep = "_"), sep = "\t", quote = F)


write.table(my.singler$singler[[2]]$SingleR.single.main$labels,
            file = paste(Sys.Date(),"JAX_Peritoneal_10X_SingleR_annotated_cells_MouseRNAseq_Main.txt",sep = "_"), sep = "\t", quote = F)
write.table(my.singler$singler[[2]]$SingleR.single$labels,
            file = paste(Sys.Date(),"JAX_Peritoneal_10X_SingleR_annotated_cells_MouseRNAseq_All.txt",sep = "_"), sep = "\t", quote = F)

# Transfer SingleR annotations to Seurat object
sum(rownames(JAX.singlets@meta.data) == rownames(my.singler$singler[[1]]$SingleR.single$labels)) # 15055
dim(JAX.singlets@meta.data) # 15055    25

JAX.singlets@meta.data$SingleR_ImmGen          <- my.singler$singler[[1]]$SingleR.single.main$labels
JAX.singlets@meta.data$SingleR_ImmGen_detailed <- my.singler$singler[[1]]$SingleR.single$labels
JAX.singlets@meta.data$SingleR_RNAseq          <- my.singler$singler[[2]]$SingleR.single.main$labels
JAX.singlets@meta.data$SingleR_RNAseq_detailed <- my.singler$singler[[2]]$SingleR.single$labels

save(JAX.singlets, file = paste(Sys.Date(),"JAX_Peritoneal_10X_Singlets.Seurat_SingleR_annotated.object.RData",sep = "_"))

# How does cluster membership vary by replicate?
table(JAX.singlets@meta.data$SingleR_ImmGen, JAX.singlets@meta.data$Sample)
#                JAX_OF1 JAX_OF2 JAX_OF3 JAX_OM1 JAX_OM2 JAX_OM3 JAX_YF1 JAX_YF2 JAX_YF3 JAX_YM1 JAX_YM2 JAX_YM3
#    B cells         853     660    1052     729     704     701     496     712     701     440     457     608
#    Basophils         0       0       1       0       0       0       0       0       0       0       0       0
#    DC                4       8       7       8      14       5      11      11      11       8      11      15
#    ILC               1       1       0       2       2       1       3       1       1       1       3       1
#    Macrophages     310     385     222     310     349     354     505     586     556     530     843     615
#    Mast cells        1       0       0       0       0       0       0       0       2       0       0       1
#    Monocytes        15       8       6      11       7       6      11      17      22      22      18      24
#    Neutrophils       0       2       3       0       0       0       0       0       0       0       0       0
#    NK cells          5       3       1       3       3       2       2       6       3       7       2      13
#    NKT               1       0       3       1       1       0       3       5       2       2       3       1
#    T cells          60      42      48     125     119     186      42      78      43      39      52      83
#    Tgd               6       1       9      12       6      10       8       6       4       3       4      12

my.freq.table <- prop.table(x = table(JAX.singlets@meta.data$SingleR_ImmGen, JAX.singlets@meta.data$Sample), margin = 2)
my.freq.table <- my.freq.table[,c(7:9,1:3,10:12,4:6)]
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

pdf(paste(Sys.Date(),"JAX_Peritoneal_10X_SingleR_IMMGEN_calls_barplot_sorted_Singlets.pdf",sep = "_"), height = 10, width = 15)
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

my.Bcell.freqs <- data.frame("Bcell_YF" = my.freq.table["B cells",1:3  ],
                             "Bcell_OF" = my.freq.table["B cells",4:6  ],
                             "Bcell_YM" = my.freq.table["B cells",7:9  ],
                             "Bcell_OM" = my.freq.table["B cells",10:12])

my.Mph.freqs <- data.frame("Mph_YF" = my.freq.table["Macrophages",1:3  ],
                           "Mph_OF" = my.freq.table["Macrophages",4:6  ],
                           "Mph_YM" = my.freq.table["Macrophages",7:9  ],
                           "Mph_OM" = my.freq.table["Macrophages",10:12])


my.Tcell.freqs <- data.frame("Tcell_YF" = my.freq.table["T cells",1:3  ],
                             "Tcell_OF" = my.freq.table["T cells",4:6  ],
                             "Tcell_YM" = my.freq.table["T cells",7:9  ],
                             "Tcell_OM" = my.freq.table["T cells",10:12])

pdf(paste(Sys.Date(),"JAX_Peritoneal_10X_Bcell_Tcell_Mph_freq_boxplot.pdf",sep = "_"), height = 4, width = 8)
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
JAX.singlets.Mph       <- subset(JAX.singlets, subset = SingleR_ImmGen %in% "Macrophages")    # Macrophages, 5565 cells
JAX.singlets.Bcells    <- subset(JAX.singlets, subset = SingleR_ImmGen %in% "B cells")        # B-cells, 8113 cells
JAX.singlets.Tcells    <- subset(JAX.singlets, subset = SingleR_ImmGen %in% "T cells")        # T-cells, 917 cells
JAX.singlets.Monocytes <- subset(JAX.singlets, subset = SingleR_ImmGen %in% "Monocytes")      # Monocytes, 167 cells

save(JAX.singlets.Mph,      file = paste(Sys.Date(),"Seurat_10x_peritoneal_Singlets_SingleR_IMMGEN_Macrophages.RData",sep = "_"))
save(JAX.singlets.Bcells,   file = paste(Sys.Date(),"Seurat_10x_peritoneal_Singlets_SingleR_IMMGEN_Bcells.RData",sep = "_"))
save(JAX.singlets.Tcells,   file = paste(Sys.Date(),"Seurat_10x_peritoneal_Singlets_SingleR_IMMGEN_Tcells.RData",sep = "_"))
save(JAX.singlets.Monocytes, file = paste(Sys.Date(),"Seurat_10x_peritoneal_Singlets_SingleR_IMMGEN_Monocytes.RData",sep = "_"))


################## Repeat with RNA-seq annotation
# How does cluster membership vary by replicate?
table(JAX.singlets@meta.data$SingleR_RNAseq, JAX.singlets@meta.data$Sample)
#                 JAX_OF1 JAX_OF2 JAX_OF3 JAX_OM1 JAX_OM2 JAX_OM3 JAX_YF1 JAX_YF2 JAX_YF3 JAX_YM1 JAX_YM2 JAX_YM3
#    B cells          841     646    1047     718     688     690     493     710     692     437     456     604
#    Granulocytes       0       2       3       0       0       0       0       0       0       0       0       0
#    Macrophages      319     396     225     317     355     353     500     578     559     524     838     608
#    Monocytes         21      17      15      22      23      12      28      34      39      37      34      50
#    NK cells           1       3       1       2       2       2       4       8       3       7       3       9
#    T cells           74      46      61     142     137     208      56      92      52      47      62     102


my.freq.table <- prop.table(x = table(JAX.singlets@meta.data$SingleR_RNAseq, JAX.singlets@meta.data$Sample), margin = 2)
my.freq.table <- my.freq.table[,c(7:9,1:3,10:12,4:6)]
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

pdf(paste(Sys.Date(),"JAX_Peritoneal_10X_SingleR_MouseRNAseq_calls_barplot_sorted_Singlets.pdf",sep = "_"), height = 10, width = 15)
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


my.Bcell.freqs <- data.frame("Bcell_YF" = my.freq.table["B cells",1:3  ],
                             "Bcell_OF" = my.freq.table["B cells",4:6  ],
                             "Bcell_YM" = my.freq.table["B cells",7:9  ],
                             "Bcell_OM" = my.freq.table["B cells",10:12])

my.Mph.freqs <- data.frame("Mph_YF" = my.freq.table["Macrophages",1:3  ],
                           "Mph_OF" = my.freq.table["Macrophages",4:6  ],
                           "Mph_YM" = my.freq.table["Macrophages",7:9  ],
                           "Mph_OM" = my.freq.table["Macrophages",10:12])


my.Tcell.freqs <- data.frame("Tcell_YF" = my.freq.table["T cells",1:3  ],
                             "Tcell_OF" = my.freq.table["T cells",4:6  ],
                             "Tcell_YM" = my.freq.table["T cells",7:9  ],
                             "Tcell_OM" = my.freq.table["T cells",10:12])

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
JAX.singlets.Mph.v2       <- subset(JAX.singlets, subset = SingleR_RNAseq %in% "Macrophages")    # Macrophages, 5572 cells
JAX.singlets.Bcells.v2    <- subset(JAX.singlets, subset = SingleR_RNAseq %in% "B cells")        # B-cells, 8022 cells
JAX.singlets.Tcells.v2    <- subset(JAX.singlets, subset = SingleR_RNAseq %in% "T cells")        # T-cells, 1079 cells
JAX.singlets.Monocytes.v2 <- subset(JAX.singlets, subset = SingleR_RNAseq %in% "Monocytes")      # Monocytes, 332 cells

save(JAX.singlets.Mph.v2       , file = paste(Sys.Date(),"Seurat_10x_peritoneal_Singlets_SingleR_MouseRNAseq_Macrophages.RData",sep = "_"))
save(JAX.singlets.Bcells.v2    , file = paste(Sys.Date(),"Seurat_10x_peritoneal_Singlets_SingleR_MouseRNAseq_Bcells.RData",sep = "_"))
save(JAX.singlets.Tcells.v2    , file = paste(Sys.Date(),"Seurat_10x_peritoneal_Singlets_SingleR_MouseRNAseq_Tcells.RData",sep = "_"))
save(JAX.singlets.Monocytes.v2 , file = paste(Sys.Date(),"Seurat_10x_peritoneal_Singlets_SingleR_MouseRNAseq_Monocytes.RData",sep = "_"))
############################################################################################################

############################################################################################################
##### 10. Analyze SingleR annotations, compare the 2 libraries/pokedex for cell types

table(JAX.singlets@meta.data$SingleR_RNAseq, JAX.singlets@meta.data$SingleR_ImmGen)
#                 B cells Basophils   DC  ILC Macrophages Mast cells Monocytes Neutrophils NK cells  NKT T cells  Tgd
#    B cells         7999         0   11    0           9          1         1           0        1    0       0    0
#    Granulocytes       0         0    0    0           0          0         0           5        0    0       0    0
#    Macrophages       93         0    0    0        5475          3         0           0        0    0       1    0
#    Monocytes          0         1  100    0          65          0       166           0        0    0       0    0
#    NK cells           0         0    0    1           2          0         0           0       39    3       0    0
#    T cells           21         0    2   16          14          0         0           0       10   19     916   81


# label cells consistently called across both sets
sum(JAX.singlets@meta.data$SingleR_RNAseq == JAX.singlets@meta.data$SingleR_ImmGen) # 14375 (out of 14924)

# keep only cells with consistent labeling
JAX.singlets.consistent       <- subset(JAX.singlets, subset = SingleR_RNAseq == SingleR_ImmGen )   # 14375
save(JAX.singlets.consistent , file = paste(Sys.Date(),"Seurat_JAX_10x_peritoneal_Singlets_SingleR_ConsistentID_SingleR.RData",sep = "_"))

# get proportions
my.freq.table <- prop.table(x = table(JAX.singlets.consistent@meta.data$SingleR_RNAseq, JAX.singlets.consistent@meta.data$Sample), margin = 2)
my.freq.table <- my.freq.table[,c(7:9,1:3,10:12,4:6)]
my.freq.table.av <- apply(my.freq.table,1,mean)
my.freq.table.av.sort <- sort(my.freq.table.av, decreasing = T,index.return = T)

my.Bcell.freqs <- data.frame("Bcell_YF" = my.freq.table["B cells",1:3  ],
                             "Bcell_OF" = my.freq.table["B cells",4:6  ],
                             "Bcell_YM" = my.freq.table["B cells",7:9  ],
                             "Bcell_OM" = my.freq.table["B cells",10:12])

my.Mph.freqs <- data.frame("Mph_YF" = my.freq.table["Macrophages",1:3  ],
                           "Mph_OF" = my.freq.table["Macrophages",4:6  ],
                           "Mph_YM" = my.freq.table["Macrophages",7:9  ],
                           "Mph_OM" = my.freq.table["Macrophages",10:12])


my.Tcell.freqs <- data.frame("Tcell_YF" = my.freq.table["T cells",1:3  ],
                             "Tcell_OF" = my.freq.table["T cells",4:6  ],
                             "Tcell_YM" = my.freq.table["T cells",7:9  ],
                             "Tcell_OM" = my.freq.table["T cells",10:12])

pdf(paste(Sys.Date(),"Bcell_Tcell_Mph_freq_boxplot_SingleR_CONSISTENT.pdf",sep = "_"), height = 4, width = 8)
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

JAX.singlets.Mph.v3       <- subset(JAX.singlets.consistent, subset = SingleR_ImmGen %in% "Macrophages")    # Macrophages, 5475 cells
JAX.singlets.Bcells.v3    <- subset(JAX.singlets.consistent, subset = SingleR_ImmGen %in% "B cells")        # B cells, 7999 cells

save(JAX.singlets.Mph.v3       , file = paste(Sys.Date(),"Seurat_10x_peritoneal_Singlets_SingleR_Consistent_Macrophages.RData",sep = "_"))
save(JAX.singlets.Bcells.v3    , file = paste(Sys.Date(),"Seurat_10x_peritoneal_Singlets_SingleR_Consistent_Bcells.RData",sep = "_"))

######################################
# Calculate and plot markers (on consistently called cells)
JAX.singlets.consistent  <- SetIdent(JAX.singlets.consistent , value = "SingleR_ImmGen")

table(JAX.singlets.consistent@meta.data$SingleR_ImmGen)
#     B cells Macrophages   Monocytes    NK cells     T cells
#        7999        5475         166          39         916 

#### get marker heatmap
perit.markers <- FindAllMarkers(JAX.singlets.consistent, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10.perit.markers <- perit.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

png(paste0(Sys.Date(),"_SingleRconsistent_Peritoneal_aging_marker_heatmap.png"),width = 25, height = 18, units = "cm", res = 300)
DoHeatmap(JAX.singlets.consistent, features = top10.perit.markers$gene)  + theme(axis.text.y = element_text(size = 8))
dev.off()

pdf(paste0(Sys.Date(),"SingleRconsistent_Peritoneal_aging_cluster_marker_Dotplot.pdf"), width = 20, height = 8)
DotPlot(JAX.singlets.consistent, features = unique(top10.perit.markers$gene), cols = c("deeppink1","deeppink4","deepskyblue1","deepskyblue4"), dot.scale = 8, split.by = "Condition") + RotatedAxis()
dev.off()

pdf(paste0(Sys.Date(),"SingleRconsistent_Peritoneal_aging_KNOWN_MARKERS_Dotplot.pdf"), width = 10, height = 8)
DotPlot(JAX.singlets.consistent, features = c("Adgre1","Itgam", "Bhlhe40", "Spi1","Cd79a","Cd79b","Cd19","Cd3d","Cd3e","Cd3g","Ccr2","Ly6c2"), cols = c("deeppink1","deeppink4","deepskyblue1","deepskyblue4"), dot.scale = 8, split.by = "Condition") + RotatedAxis()
dev.off()

png(paste0(Sys.Date(),"_SingleRconsistent_Peritoneal_aging_KNOWN_MARKER_heatmap.png"),width = 25, height = 18, units = "cm", res = 300)
DoHeatmap(JAX.singlets.consistent,
          features = c("Cd79a","Cd79b","Cd19","Cd3d","Cd3e","Cd3g","Adgre1","Itgam", "Bhlhe40", "Spi1","Ccr2","Ly6c2"))  + theme(axis.text.y = element_text(size = 8))
dev.off()
############################################################################################################

#######################
sink(file = paste0(Sys.Date(),"_JAX_Perit_Seurat_Singlet_processing_session_Info.txt"))
sessionInfo()
sink()


