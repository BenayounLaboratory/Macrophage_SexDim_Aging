setwd('/Volumes/BB_HQ_1/Immune_sex_dimorphism_Aging/Macrophages/Clean_Macrophage_Data/Peritoneal_scRNAseq/v2_set/Seurat_Doublet_Identification_clustering/')
options(stringsAsFactors = F)

# General use packages
library('Seurat')
library(bitops)
library(sctransform)
library(SingleR)
library(pheatmap)
library(dplyr)
library('ggplot2')

# packages for doublet calling
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
cellRanger.perit_YM1  <- load10X("/Volumes/BB_HQ_1/Immune_sex_dimorphism_Aging/Macrophages/Clean_Macrophage_Data/Peritoneal_scRNAseq/v2_set/Cell_Ranger/4m_Male_perit_pool_1/outs/"    )
cellRanger.perit_YF1  <- load10X("/Volumes/BB_HQ_1/Immune_sex_dimorphism_Aging/Macrophages/Clean_Macrophage_Data/Peritoneal_scRNAseq/v2_set/Cell_Ranger/4m_Female_perit_pool_1/outs/"  )
cellRanger.perit_OM1  <- load10X("/Volumes/BB_HQ_1/Immune_sex_dimorphism_Aging/Macrophages/Clean_Macrophage_Data/Peritoneal_scRNAseq/v2_set/Cell_Ranger/20m_Male_perit_pool_1/outs/"   )
cellRanger.perit_OF1  <- load10X("/Volumes/BB_HQ_1/Immune_sex_dimorphism_Aging/Macrophages/Clean_Macrophage_Data/Peritoneal_scRNAseq/v2_set/Cell_Ranger/20m_Female_perit_pool_1/outs/" )
cellRanger.perit_YM2  <- load10X("/Volumes/BB_HQ_1/Immune_sex_dimorphism_Aging/Macrophages/Clean_Macrophage_Data/Peritoneal_scRNAseq/v2_set/Cell_Ranger/4m_Male_perit_pool_2/outs/"    )
cellRanger.perit_YF2  <- load10X("/Volumes/BB_HQ_1/Immune_sex_dimorphism_Aging/Macrophages/Clean_Macrophage_Data/Peritoneal_scRNAseq/v2_set/Cell_Ranger/4m_Female_perit_pool_2/outs/"  )
cellRanger.perit_OM2  <- load10X("/Volumes/BB_HQ_1/Immune_sex_dimorphism_Aging/Macrophages/Clean_Macrophage_Data/Peritoneal_scRNAseq/v2_set/Cell_Ranger/20m_Male_perit_pool_2/outs/"   )
cellRanger.perit_OF2  <- load10X("/Volumes/BB_HQ_1/Immune_sex_dimorphism_Aging/Macrophages/Clean_Macrophage_Data/Peritoneal_scRNAseq/v2_set/Cell_Ranger/20m_Female_perit_pool_2/outs/" )

# estimate RNA soup from bottom 10% droplets
cellRanger.perit_YM1  <- autoEstCont(cellRanger.perit_YM1 ) # Estimated global rho of 0.03
cellRanger.perit_YF1  <- autoEstCont(cellRanger.perit_YF1 ) # Estimated global rho of 0.03
cellRanger.perit_OM1  <- autoEstCont(cellRanger.perit_OM1 ) # Estimated global rho of 0.02
cellRanger.perit_OF1  <- autoEstCont(cellRanger.perit_OF1 ) # Estimated global rho of 0.01
cellRanger.perit_YM2  <- autoEstCont(cellRanger.perit_YM2 ) # Estimated global rho of 0.04
cellRanger.perit_YF2  <- autoEstCont(cellRanger.perit_YF2 ) # Estimated global rho of 0.02
cellRanger.perit_OM2  <- autoEstCont(cellRanger.perit_OM2 ) # Estimated global rho of 0.04
cellRanger.perit_OF2  <- autoEstCont(cellRanger.perit_OF2 ) # Estimated global rho of 0.02

# adjust counts based on RNA soup (and round)
out.perit_YM1 <- adjustCounts(cellRanger.perit_YM1, roundToInt=TRUE)
out.perit_YF1 <- adjustCounts(cellRanger.perit_YF1, roundToInt=TRUE)
out.perit_OM1 <- adjustCounts(cellRanger.perit_OM1, roundToInt=TRUE)
out.perit_OF1 <- adjustCounts(cellRanger.perit_OF1, roundToInt=TRUE)
out.perit_YM2 <- adjustCounts(cellRanger.perit_YM2, roundToInt=TRUE)
out.perit_YF2 <- adjustCounts(cellRanger.perit_YF2, roundToInt=TRUE)
out.perit_OM2 <- adjustCounts(cellRanger.perit_OM2, roundToInt=TRUE)
out.perit_OF2 <- adjustCounts(cellRanger.perit_OF2, roundToInt=TRUE)

# get seurat objects
seurat.perit_YM1 <- CreateSeuratObject( out.perit_YM1 )
seurat.perit_YF1 <- CreateSeuratObject( out.perit_YF1 )
seurat.perit_OM1 <- CreateSeuratObject( out.perit_OM1 )
seurat.perit_OF1 <- CreateSeuratObject( out.perit_OF1 )
seurat.perit_YM2 <- CreateSeuratObject( out.perit_YM2 )
seurat.perit_YF2 <- CreateSeuratObject( out.perit_YF2 )
seurat.perit_OM2 <- CreateSeuratObject( out.perit_OM2 )
seurat.perit_OF2 <- CreateSeuratObject( out.perit_OF2 )


# Merge cleaned seurat objects
my.perit <- merge(seurat.perit_YM1, 
                  y =  c(seurat.perit_YF1,
                         seurat.perit_OM1,
                         seurat.perit_OF1,
                         seurat.perit_YM2,
                         seurat.perit_YF2,
                         seurat.perit_OM2,
                         seurat.perit_OF2), 
                  add.cell.ids = c("perit_YM1",
                                   "perit_YF1",
                                   "perit_OM1",
                                   "perit_OF1",
                                   "perit_YM2",
                                   "perit_YF2",
                                   "perit_OM2",
                                   "perit_OF2"), 
                  project = "10x_NIA_peritoneal")
my.perit
# An object of class Seurat 
# 32285 features across 20750 samples within 1 assay 
# Active assay: RNA (32285 features, 0 variable features)


###############
# create library label
Sample_ID <- rep("NA", length(colnames(my.perit@assays$RNA)))
Sample_ID[ grep("perit_YM1"  , colnames(my.perit@assays$RNA)) ]   <- "perit_YM1"
Sample_ID[ grep("perit_YF1"  , colnames(my.perit@assays$RNA)) ]   <- "perit_YF1"
Sample_ID[ grep("perit_OM1"  , colnames(my.perit@assays$RNA)) ]   <- "perit_OM1"
Sample_ID[ grep("perit_OF1"  , colnames(my.perit@assays$RNA)) ]   <- "perit_OF1"
Sample_ID[ grep("perit_YM2"  , colnames(my.perit@assays$RNA)) ]   <- "perit_YM2"
Sample_ID[ grep("perit_YF2"  , colnames(my.perit@assays$RNA)) ]   <- "perit_YF2"
Sample_ID[ grep("perit_OM2"  , colnames(my.perit@assays$RNA)) ]   <- "perit_OM2"
Sample_ID[ grep("perit_OF2"  , colnames(my.perit@assays$RNA)) ]   <- "perit_OF2"
Sample_ID <- data.frame(Sample_ID)
rownames(Sample_ID) <- colnames(my.perit@assays$RNA)

# create batch label
Batch <- rep("NA", length(colnames(my.perit@assays$RNA)))
Batch[ grep("1_"  , colnames(my.perit@assays$RNA)) ]   <- "Cohort1"
Batch[ grep("2_"  , colnames(my.perit@assays$RNA)) ]   <- "Cohort2"
Batch <- data.frame(Batch)
rownames(Batch) <- colnames(my.perit@assays$RNA)

# create condition label
Condition <- rep("NA", length(colnames(my.perit@assays$RNA)))
Condition[ grep("perit_YM"  , colnames(my.perit@assays$RNA)) ]   <- "YM"
Condition[ grep("perit_YF"  , colnames(my.perit@assays$RNA)) ]   <- "YF"
Condition[ grep("perit_OM"  , colnames(my.perit@assays$RNA)) ]   <- "OM"
Condition[ grep("perit_OF"  , colnames(my.perit@assays$RNA)) ]   <- "OF"
Condition <- data.frame(Condition)
rownames(Condition) <- colnames(my.perit@assays$RNA)


# update Seurat with metadata
my.perit <- AddMetaData(object = my.perit, metadata = as.vector(Sample_ID)    , col.name = "Sample_ID"       )
my.perit <- AddMetaData(object = my.perit, metadata = as.vector(Batch)        , col.name = "Batch"           )
my.perit <- AddMetaData(object = my.perit, metadata = as.vector(Condition)    , col.name = "Condition"       )

###############

head(my.perit@meta.data)
#                                 orig.ident nCount_RNA nFeature_RNA Sample_ID   Batch Condition
# perit_YM1_AAACCTGCAGCCTATA-1 SeuratProject       1900          719 perit_YM1 Cohort1        YM
# perit_YM1_AAACGGGAGGCATTGG-1 SeuratProject       6016         1559 perit_YM1 Cohort1        YM
# perit_YM1_AAACGGGAGTGTCCAT-1 SeuratProject       6039         1880 perit_YM1 Cohort1        YM
# perit_YM1_AAACGGGTCTGAGGGA-1 SeuratProject       4232         1104 perit_YM1 Cohort1        YM
# perit_YM1_AAAGATGAGAATGTGT-1 SeuratProject       5130         1602 perit_YM1 Cohort1        YM
# perit_YM1_AAAGATGAGATCCGAG-1 SeuratProject       9266         1592 perit_YM1 Cohort1        YM

my.perit <- SetIdent(my.perit, value = "Sample_ID")


################################################################################################################################################################
#### 2. QC on mitochondrial reads and depth
# The number of genes and UMIs (nGene and nUMI) are automatically calculated for every object by Seurat.  
# The % of UMI mapping to MT-genes is a common scRNA-seq QC metric.
### QC on mitochondrial reads
my.perit[["percent.mito"]] <- PercentageFeatureSet(my.perit, pattern = "^mt-")

pdf(paste(Sys.Date(),"Peritoneal_cells_aging_violinPlots_QC_gene_UMI_mito.pdf", sep = "_"), height = 5, width = 10)
VlnPlot(object = my.perit, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
dev.off()

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(my.perit, feature1 = "nCount_RNA", feature2 = "percent.mito")
plot2 <- FeatureScatter(my.perit, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

pdf(paste(Sys.Date(),"Peritoneal_cells_aging_QC_scatter.pdf", sep = "_"), height = 5, width = 10)
plot1 + plot2
dev.off()

# We filter out cells that have unique gene counts or less than 300  
my.perit <- subset(my.perit, subset = nFeature_RNA > 300 & percent.mito < 25 & nFeature_RNA < 3500)
my.perit
# An object of class Seurat 
# 32285 features across 19839 samples within 1 assay 
# Active assay: RNA (32285 features, 0 variable features)

plot1 <- FeatureScatter(my.perit, feature1 = "nCount_RNA", feature2 = "percent.mito")
plot2 <- FeatureScatter(my.perit, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

pdf(paste(Sys.Date(),"Peritoneal_cells_aging_QC_scatter_post_filter.pdf", sep = "_"), height = 5, width = 10)
plot1 + plot2
dev.off()

################################################################################################################################################################
#### 3. Normalizing the data
# global-scaling normalization method ???LogNormalize??? that normalizes the gene expression measurements for each cell by the total expression, 
# multiplies this by a scale factor (10,000 by default), and log-transforms the result.
my.perit <- NormalizeData(object = my.perit, normalization.method = "LogNormalize",  scale.factor = 10000)

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
my.perit <- CellCycleScoring(object = my.perit, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# view cell cycle scores and phase assignments
head(x = my.perit@meta.data)

# write predictions to file
write.table(my.perit@meta.data, file = paste(Sys.Date(),"Peritoneal_cells_aging_CellCycle_predictions.txt", sep = "_"), sep = "\t", quote = F)


pdf(paste0(Sys.Date(), "Peritoneal_aging_cell_cycle_pie_charts.pdf"), height = 13, width = 8)
par(mfrow=c(2,2))
par(oma=c(0.1,0.1,0.1,0.1))
pie(table(my.perit@meta.data$Phase[my.perit@meta.data$Condition == "YF"]) , main = "YF")
pie(table(my.perit@meta.data$Phase[my.perit@meta.data$Condition == "YM"]) , main = "YM")
pie(table(my.perit@meta.data$Phase[my.perit@meta.data$Condition == "OF"]) , main = "OF")
pie(table(my.perit@meta.data$Phase[my.perit@meta.data$Condition == "OM"]) , main = "OM")
par(mfrow=c(1,1))
dev.off()

pdf(paste0(Sys.Date(), "Peritoneal_aging_cell_cycle_pie_charts_by_batch.pdf"), height = 13, width = 8)
par(mfrow=c(4,2))
par(oma=c(0.1,0.1,0.1,0.1))
pie(table(my.perit@meta.data$Phase[my.perit@meta.data$Sample_ID == "perit_YM1"]) , main = "YM1")
pie(table(my.perit@meta.data$Phase[my.perit@meta.data$Sample_ID == "perit_YF1"]) , main = "YF1")
pie(table(my.perit@meta.data$Phase[my.perit@meta.data$Sample_ID == "perit_OM1"]) , main = "OM1")
pie(table(my.perit@meta.data$Phase[my.perit@meta.data$Sample_ID == "perit_OF1"]) , main = "OF1")
pie(table(my.perit@meta.data$Phase[my.perit@meta.data$Sample_ID == "perit_YM2"]), main =  "YM2")
pie(table(my.perit@meta.data$Phase[my.perit@meta.data$Sample_ID == "perit_YF2"]), main =  "YF2")
pie(table(my.perit@meta.data$Phase[my.perit@meta.data$Sample_ID == "perit_OM2"]), main =  "OM2")
pie(table(my.perit@meta.data$Phase[my.perit@meta.data$Sample_ID == "perit_OF2"]), main =  "OF2")
par(mfrow=c(1,1))
dev.off()



################################################################################################################################################################
##### 5.Find and remove doublets
# https://github.com/chris-mcginnis-ucsf/DoubletFinder
my.perit <- SCTransform(object = my.perit, vars.to.regress = c("nFeature_RNA", "percent.mito", "Batch"))
my.perit <- RunPCA(my.perit)

# Determine the ‘dimensionality’ of the dataset
pdf(paste0(Sys.Date(), "_Peritoneal_Lavage_NIA_elbowplot.pdf"))
ElbowPlot(my.perit, ndims = 50)
dev.off()

################################################################################
# https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html
# To give us an idea of the number of PCs needed to be included:
# We can calculate where the principal components start to elbow by taking the larger value of:
#    - The point where the principal components only contribute 5% of standard deviation
#    - The principal components cumulatively contribute 90% of the standard deviation.
#    - The point where the percent change in variation between the consecutive PCs is less than 0.1%.

# Determine percent of variation associated with each PC
pct <- my.perit[["pca"]]@stdev / sum(my.perit[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
co1 # 42

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
pdf(paste0(Sys.Date(), "_Peritoneal_Lavage_NIA_elbowplot_threshold_analysis.pdf"), height = 5, width= 6)
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) +
  geom_text() +
  geom_vline(xintercept = 90, color = "grey") +
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()
dev.off()
###############################################################################


my.perit <- RunUMAP(my.perit, dims = 1:pcs)
my.perit <- FindNeighbors(my.perit, dims = 1:pcs)
my.perit <- FindClusters(object = my.perit)


####################################################################################
#### need to split by original 10x sample to make sure to identify real doublets when running doublet finding

# Setup the Seurat objects as a list
perit.list <- SplitObject(my.perit, split.by = "Sample_ID")

## Assume doublet rate based on 10x information (+ 20% cells to be safe)
pred.dblt.rate <- predict(pred_dblt_lm, data.frame("cell_number" = 1.2 * unlist(lapply(perit.list, ncol))))/100
pred.dblt.rate
# perit_YM1 perit_YF1 perit_OM1 perit_OF1 perit_YM2 perit_YF2 perit_OM2 perit_OF2 
# 0.0106464 0.0394272 0.0362976 0.0337440 0.0145440 0.0178464 0.0168288 0.0211584 

# loop over sample
for (i in 1:length(perit.list)) {
  
  ## pK Identification (no ground-truth)
  sweep.res.list <- paramSweep_v3(perit.list[[i]], PCs = 1:30, sct = TRUE, num.cores	 = 2)
  sweep.stats    <- summarizeSweep(sweep.res.list , GT = FALSE)
  bcmvn          <- find.pK(sweep.stats)
  
  # need some R gymnastics since the Pk is stored as a factor for some reason
  # to get the pK number, need to first convert to character and THEN to numeric
  # numeric first yeild row number
  pk.perit <- as.numeric(as.character(bcmvn[as.numeric(bcmvn$pK[bcmvn$BCmetric == max(bcmvn$BCmetric)]),"pK"]))
  
  ## Homotypic Doublet Proportion Estimate
  homotypic.prop <- modelHomotypic(perit.list[[i]]@meta.data$seurat_clusters)     ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
  nExp_poi       <- round((pred.dblt.rate[i])*length(perit.list[[i]]@meta.data$orig.ident))        ## Assuming doublet formation rate based on 10x manual
  
  ## Run DoubletFinder with varying classification stringencies
  perit.list[[i]] <- doubletFinder_v3(perit.list[[i]], PCs = 1:30, pN = 0.25, pK = pk.perit, nExp = nExp_poi,     reuse.pANN = FALSE, sct = T)
  
  # rename column to enable subsetting
  colnames(perit.list[[i]]@meta.data)[grep("DF.classifications_0.25",colnames(perit.list[[i]]@meta.data))] <- "DoubletFinder"
  
}

# plot UMAP plots for each pool
for (i in 1:length(perit.list)) {
  
  # get classification name
  pdf(paste(Sys.Date(),"Library",names(perit.list)[i],"Doublet_Finder_UMAP.pdf", sep = "_"), height = 5, width = 5)
  print(DimPlot(perit.list[[i]], reduction = "umap", group.by = "DoubletFinder"))
  dev.off()
}

table(perit.list[[1]]@meta.data$DoubletFinder)
# Doublet Singlet
#  12    1097 

table(perit.list[[2]]@meta.data$DoubletFinder)
# Doublet Singlet
#    162    3945  

table(perit.list[[3]]@meta.data$DoubletFinder)
# Doublet Singlet
#  137    3644 

table(perit.list[[4]]@meta.data$DoubletFinder)
# Doublet Singlet
#  119    3396 

table(perit.list[[5]]@meta.data$DoubletFinder)
# Doublet Singlet
#  22    1493 

table(perit.list[[6]]@meta.data$DoubletFinder)
# Doublet Singlet
#  33    1826 

table(perit.list[[7]]@meta.data$DoubletFinder)
# Doublet Singlet
#   30    1723

table(perit.list[[8]]@meta.data$DoubletFinder)
# Doublet Singlet
#   47    2157 

# save scds
save(perit.list, file = paste0(Sys.Date(),"_NIA_Peritoneal_Lavage_Seurat_object_DoubletFinder_List.RData"))

## gate back scds annotation to original Seurat object
my.perit@meta.data$DoubletFinder <- NA # initialize

for (i in 1:length(perit.list)) {
  # for each object compare and move doublet annotations over
  my.perit@meta.data[colnames(perit.list[[i]]), ]$DoubletFinder <- perit.list[[i]]@meta.data$DoubletFinder
}



##### Inspect doublets
table(my.perit@meta.data$DoubletFinder, my.perit@meta.data$Sample_ID)
#           perit_OF1 perit_OF2 perit_OM1 perit_OM2 perit_YF1 perit_YF2 perit_YM1 perit_YM2
#   Doublet       119        47       137        30       162        33        12        22
#   Singlet      3396      2157      3644      1723      3945      1826      1097      1493

# save data for singlets
my.perit.singlets    <- subset(my.perit, subset = DoubletFinder %in% "Singlet") 
my.perit.singlets
# An object of class Seurat 
# 48287 features across 19281 samples within 2 assays 
# Active assay: SCT (16002 features, 3000 variable features)
#  1 other assay present: RNA
#  2 dimensional reductions calculated: pca, umap

pdf(paste(Sys.Date(),"Peritoneal_cells_aging_Singlets_UMAP_by_Group.pdf", sep = "_"), height = 5, width = 5)
DimPlot(my.perit.singlets, reduction = "umap", group.by = "Batch", shuffle = T)
DimPlot(my.perit.singlets, reduction = "umap", group.by = "Condition", shuffle = T)
dev.off()


################################################################################################################################################################
##### 6.Scaling the data and removing unwanted sources of variation

# Run SCT normalization, regress potential problems
my.perit.singlets <- SCTransform(object = my.perit.singlets, vars.to.regress = c("nFeature_RNA", "percent.mito", "Batch","S.Score", "G2M.Score"))

# find variable genes
my.perit.singlets <- FindVariableFeatures(my.perit.singlets, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(my.perit.singlets), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(my.perit.singlets)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

pdf(paste(Sys.Date(),"Peritoneal_cells_aging_Singlets_Variable genes.pdf", sep = "_"), height = 5, width = 10)
plot1 + plot2
dev.off()

################################################################################################################################################################
##### 7. Cluster the cells

# perform clustering again now that only singlets remain
my.perit.singlets <- FindNeighbors(my.perit.singlets, dims = 1:16)

my.perit.singlets <- FindClusters(object         = my.perit.singlets,
                                  reduction.type = "pca",
                                  k.param        = 100,
                                  dims.use       = 1:16,
                                  resolution     = 0.3,
                                  print.output   = 1,
                                  save.SNN       = TRUE,
                                  plot.SNN       = TRUE,
                                  n.start        = 100)
# Number of communities: 12

pdf(paste(Sys.Date(),"Peritoneal_cells_aging_UMAP_clusters_0.3res_Singlets.pdf", sep = "_"), height = 5, width = 6.5)
DimPlot(my.perit.singlets, label = TRUE)
dev.off()

write.table(my.perit.singlets@meta.data, file = paste(Sys.Date(),"Seurat_PCA_SNN_clusters_and_Metadata_Singlets.txt", sep = "_"), sep = "\t", quote = F, col.names = F, row.names = T)


################################################################################################################################################################
##### 7. Visualization of potential marker gene expression

pdf(paste(Sys.Date(),"Peritoneal_cells_aging_Feature_plot_Singlets.pdf", sep = "_"), width = 15, height = 15)
FeaturePlot(object = my.perit.singlets, features = c("Cd79a","Adgre1","Itgam","Bhlhe40","Cd14","Cd3e","Spi1","Cd209a")  )
FeaturePlot(object = my.perit.singlets, features = c("Cd86", "Cd68", "Cd80","Tlr2", "Tlr4")  ) ### M1 "Cx3cr1","Nos2","Il1r", "Il1b"
FeaturePlot(object = my.perit.singlets, features = c("Arg1", "Cd163")  ) # m2 "Cd206","Fizz1" Ym1/2
dev.off()

pdf(paste(Sys.Date(),"Peritoneal_cells_aging_Xist_ridge_plot_Singlets.pdf", sep = "_"))
RidgePlot(object = my.perit.singlets, features = "Xist", ncol = 1, group.by = "Condition")
dev.off()

pdf(paste(Sys.Date(),"Peritoneal_cells_aging_Y_genes_ridge_plot_Singlets.pdf", sep = "_"), height = 5, width = 10)
RidgePlot(object = my.perit.singlets, features = c("Ddx3y", "Eif2s3y"), ncol = 2, group.by = "Condition")
dev.off()


save(my.perit.singlets, file = paste(Sys.Date(),"Peritoneal_10X_Singlets.Seurat.object.RData",sep = "_"))

################################################################################################################################################################
##### 8. Identify cell types using SingleR

library(SingleR) # SingleR_1.0.1

my.perit.singlets@meta.data$Sample_ID <- paste0(my.perit.singlets@meta.data$Condition, my.perit.singlets@meta.data$Batch)

################ run singleR and save object
my.singler = CreateSinglerObject(counts    = my.perit.singlets@assays$SCT@counts,
                                 annot     = my.perit.singlets@meta.data$Sample_ID,
                                 min.genes = 250, technology = "10X",
                                 species   = "Mouse", project.name = "10X_Peritoneal")

my.singler$seurat = my.perit.singlets # (optional)
my.singler$meta.data$orig.ident = my.perit.singlets@meta.data$orig.ident # the original identities, if not supplied in 'annot'

## if using Seurat v3.0 and over use:
my.singler$meta.data$xy       =    my.perit.singlets@reductions$umap@cell.embeddings # the Umap coordinates
my.singler$meta.data$clusters =    my.perit.singlets@active.ident # the Seurat clusters (if 'clusters' not provided)
my.singler$meta.data$Sample   =    my.perit.singlets@meta.data$Sample_ID # the Seurat clusters (if 'clusters' not provided)

save(my.singler, file = paste(Sys.Date(),"Peritoneal_10X_SingleR_object_SCT.RData",sep = "_"))

dev.off()

#######################
### Main cell types
pdf(paste(Sys.Date(),"Peritoneal_10X_SingleR_object.main.cell.types.pdf",sep = "_"), height = 5, width = 7)
SingleR.DrawHeatmap(my.singler$singler[[1]]$SingleR.single.main, top.n = 50, clusters = my.singler$meta.data$Sample)
dev.off()

### All cell types
pdf(paste(Sys.Date(),"Peritoneal_10X_SingleR_object.all.cell.types.pdf",sep = "_"), height = 5, width = 7)
SingleR.DrawHeatmap(my.singler$singler[[1]]$SingleR.single, top.n = 50, clusters = my.singler$meta.data$Sample)
dev.off()

### Main cell types
pdf(paste(Sys.Date(),"Peritoneal_10X_SingleR_object.main.cell.types_MouseRNAseq.pdf",sep = "_"), height = 5, width = 7)
SingleR.DrawHeatmap(my.singler$singler[[2]]$SingleR.single.main, top.n = Inf,clusters = my.singler$meta.data$Sample)
dev.off()

### All cell types
pdf(paste(Sys.Date(),"Peritoneal_10X_SingleR_object.all.cell.types_MouseRNAseq.pdf",sep = "_"), height = 5, width = 7)
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

pdf(paste(Sys.Date(),"Peritoneal_10X_SingleR_tSNE_annotated_IMMGEN.pdf",sep = "_"), height = 10, width = 12)
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

pdf(paste(Sys.Date(),"Peritoneal_10X_SingleR_tSNE_annotated_MouseRNAseq.pdf",sep = "_"), height = 10, width = 12)
out$p
dev.off()

write.table(my.singler$singler[[1]]$SingleR.single.main$labels,
            file = paste(Sys.Date(),"Peritoneal_10X_SingleR_annotated_cells_Immgen_Main.txt",sep = "_"), sep = "\t", quote = F)
write.table(my.singler$singler[[1]]$SingleR.single$labels,
            file = paste(Sys.Date(),"Peritoneal_10X_SingleR_annotated_cells_Immgen_All.txt",sep = "_"), sep = "\t", quote = F)


write.table(my.singler$singler[[2]]$SingleR.single.main$labels,
            file = paste(Sys.Date(),"Peritoneal_10X_SingleR_annotated_cells_MouseRNAseq_Main.txt",sep = "_"), sep = "\t", quote = F)
write.table(my.singler$singler[[2]]$SingleR.single$labels,
            file = paste(Sys.Date(),"Peritoneal_10X_SingleR_annotated_cells_MouseRNAseq_All.txt",sep = "_"), sep = "\t", quote = F)

#  Transfer SingleR annotations to Seurat object
sum(rownames(my.perit.singlets@meta.data) %in% rownames(my.singler$singler[[1]]$SingleR.single$labels)) # 19278
dim(my.perit.singlets@meta.data) # 19278    17

my.perit.singlets@meta.data$SingleR_ImmGen          <- my.singler$singler[[1]]$SingleR.single.main$labels
my.perit.singlets@meta.data$SingleR_ImmGen_detailed <- my.singler$singler[[1]]$SingleR.single$labels
my.perit.singlets@meta.data$SingleR_RNAseq          <- my.singler$singler[[2]]$SingleR.single.main$labels
my.perit.singlets@meta.data$SingleR_RNAseq_detailed <- my.singler$singler[[2]]$SingleR.single$labels

save(my.perit.singlets, file = paste(Sys.Date(),"Peritoneal_10X_Singlets.Seurat_SingleR_annotated.object.RData",sep = "_"))

################################################################################################################################################################
##### 9. Analyze SingleR annotations


# How does cluster membership vary by replicate?
table(my.perit.singlets@meta.data$SingleR_ImmGen, my.perit.singlets@meta.data$Sample_ID)
#                  OFCohort1 OFCohort2 OMCohort1 OMCohort2 YFCohort1 YFCohort2 YMCohort1 YMCohort2
#   B cells               2696      1623      2234      1053      2067      1007       577       646
#   Basophils                0         0         0         0         0         1         0         0
#   DC                      18        13        30        19        17         9        13        15
#   Epithelial cells         0         1         0         0         0         0         0         0
#   ILC                      3         2         6         6         6         7         2         6
#   Macrophages            269       272       709       361      1662       596       420       743
#   Monocytes               23        14        36        15        41        45        12        24
#   Neutrophils              2         2         2         5         1         0         0         0
#   NK cells                 8         4        12         1        11        12         1         6
#   NKT                      8         2         8         3         7         4         3         3
#   T cells                355       216       571       246       125       137        62        46
#   Tgd                     14         8        36        14         8         8         7         4

my.freq.table <- prop.table(x = table(my.perit.singlets@meta.data$SingleR_ImmGen, my.perit.singlets@meta.data$Sample_ID), margin = 2)
my.freq.table <- my.freq.table[,c(5:6,1:2,7:8,3:4)]
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

pdf(paste(Sys.Date(),"Peritoneal_10X_SingleR_IMMGEN_calls_barplot_sorted_Singlets.pdf",sep = "_"), height = 10, width = 15)
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

Bcell.YF <- my.freq.table["B cells",1:2]
Bcell.OF <- my.freq.table["B cells",3:4]
Bcell.YM <- my.freq.table["B cells",5:6]
Bcell.OM <- my.freq.table["B cells",7:8]

my.Bcell.freqs <- data.frame("Bcell_YF" = my.freq.table["B cells",1:2],
                             "Bcell_OF" = my.freq.table["B cells",3:4],
                             "Bcell_YM" = my.freq.table["B cells",5:6],
                             "Bcell_OM" = my.freq.table["B cells",7:8])

my.Mph.freqs <- data.frame("Mph_YF" = my.freq.table["Macrophages",1:2],
                           "Mph_OF" = my.freq.table["Macrophages",3:4],
                           "Mph_YM" = my.freq.table["Macrophages",5:6],
                           "Mph_OM" = my.freq.table["Macrophages",7:8])


my.Tcell.freqs <- data.frame("Tcell_YF" = my.freq.table["T cells",1:2],
                             "Tcell_OF" = my.freq.table["T cells",3:4],
                             "Tcell_YM" = my.freq.table["T cells",5:6],
                             "Tcell_OM" = my.freq.table["T cells",7:8])

pdf(paste(Sys.Date(),"Bcell_Tcell_Mph_freq_boxplot.pdf",sep = "_"), height = 4, width = 8)
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
my.perit.singlets.Mph       <- subset(my.perit.singlets, subset = SingleR_ImmGen %in% "Macrophages")    # Macrophages, 5032 cells
my.perit.singlets.Bcells    <- subset(my.perit.singlets, subset = SingleR_ImmGen %in% "B cells")        # B-cells, 11903 cells
my.perit.singlets.Tcells    <- subset(my.perit.singlets, subset = SingleR_ImmGen %in% "T cells")        # T-cells, 1758 cells
my.perit.singlets.Monocytes <- subset(my.perit.singlets, subset = SingleR_ImmGen %in% "Monocytes")      # Monocytes, 210 cells

save(my.perit.singlets.Mph,      file = paste(Sys.Date(),"Seurat_10x_peritoneal_Singlets_SingleR_IMMGEN_Macrophages.RData",sep = "_"))
save(my.perit.singlets.Bcells,   file = paste(Sys.Date(),"Seurat_10x_peritoneal_Singlets_SingleR_IMMGEN_Bcells.RData",sep = "_"))
save(my.perit.singlets.Tcells,   file = paste(Sys.Date(),"Seurat_10x_peritoneal_Singlets_SingleR_IMMGEN_Tcells.RData",sep = "_"))
save(my.perit.singlets.Monocytes, file = paste(Sys.Date(),"Seurat_10x_peritoneal_Singlets_SingleR_IMMGEN_Monocytes.RData",sep = "_"))



###### Repeat with RNA-seq annotation
# How does cluster membership vary by replicate?
table(my.perit.singlets@meta.data$SingleR_RNAseq, my.perit.singlets@meta.data$Sample_ID)
#                   OFCohort1 OFCohort2 OMCohort1 OMCohort2 YFCohort1 YFCohort2 YMCohort1 YMCohort2
#   B cells              2693      1622      2226      1049      2061      1003       578       645
#   Dendritic cells         0         0         0         0         1         0         0         0
#   Erythrocytes            0         0         0         2         1         0         0         0
#   Granulocytes            2         2         1         5         1         1         0         0
#   Macrophages           269       267       702       354      1645       583       395       721
#   Monocytes              41        32        72        39        73        67        48        61
#   NK cells                2         3         7         1        11        16         2         7
#   T cells               389       231       636       273       152       156        74        59


my.freq.table <- prop.table(x = table(my.perit.singlets@meta.data$SingleR_RNAseq, my.perit.singlets@meta.data$Sample_ID), margin = 2)
my.freq.table <- my.freq.table[,c(5:6,1:2,7:8,3:4)]
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

pdf(paste(Sys.Date(),"Peritoneal_10X_SingleR_MouseRNAseq_calls_barplot_sorted_Singlets.pdf",sep = "_"), height = 10, width = 15)
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

Bcell.YF <- my.freq.table["B cells",1:2]
Bcell.OF <- my.freq.table["B cells",3:4]
Bcell.YM <- my.freq.table["B cells",5:6]
Bcell.OM <- my.freq.table["B cells",7:8]

my.Bcell.freqs <- data.frame("Bcell_YF" = my.freq.table["B cells",1:2],
                             "Bcell_OF" = my.freq.table["B cells",3:4],
                             "Bcell_YM" = my.freq.table["B cells",5:6],
                             "Bcell_OM" = my.freq.table["B cells",7:8])

my.Mph.freqs <- data.frame("Mph_YF" = my.freq.table["Macrophages",1:2],
                           "Mph_OF" = my.freq.table["Macrophages",3:4],
                           "Mph_YM" = my.freq.table["Macrophages",5:6],
                           "Mph_OM" = my.freq.table["Macrophages",7:8])


my.Tcell.freqs <- data.frame("Tcell_YF" = my.freq.table["T cells",1:2],
                             "Tcell_OF" = my.freq.table["T cells",3:4],
                             "Tcell_YM" = my.freq.table["T cells",5:6],
                             "Tcell_OM" = my.freq.table["T cells",7:8])

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
my.perit.singlets.Mph.v2       <- subset(my.perit.singlets, subset = SingleR_RNAseq %in% "Macrophages")    # Macrophages, 4936 cells
my.perit.singlets.Bcells.v2    <- subset(my.perit.singlets, subset = SingleR_RNAseq %in% "B cells")        # B-cells, 11877 cells
my.perit.singlets.Tcells.v2    <- subset(my.perit.singlets, subset = SingleR_RNAseq %in% "T cells")        # T-cells, 1970 cells
my.perit.singlets.Monocytes.v2 <- subset(my.perit.singlets, subset = SingleR_RNAseq %in% "Monocytes")      # Monocytes, 433 cells

save(my.perit.singlets.Mph.v2       , file = paste(Sys.Date(),"Seurat_10x_peritoneal_Singlets_SingleR_MouseRNAseq_Macrophages.RData",sep = "_"))
save(my.perit.singlets.Bcells.v2    , file = paste(Sys.Date(),"Seurat_10x_peritoneal_Singlets_SingleR_MouseRNAseq_Bcells.RData",sep = "_"))
save(my.perit.singlets.Tcells.v2    , file = paste(Sys.Date(),"Seurat_10x_peritoneal_Singlets_SingleR_MouseRNAseq_Tcells.RData",sep = "_"))
save(my.perit.singlets.Monocytes.v2 , file = paste(Sys.Date(),"Seurat_10x_peritoneal_Singlets_SingleR_MouseRNAseq_Monocytes.RData",sep = "_"))
############################################################################################################

############################################################################################################
##### 10. Analyze SingleR annotations, compare the 2 libraries

table(my.perit.singlets@meta.data$SingleR_RNAseq, my.perit.singlets@meta.data$SingleR_ImmGen)
#                   B cells Basophils    DC Epithelial cells   ILC Macrophages Monocytes Neutrophils NK cells   NKT T cells   Tgd
#    B cells           11866         0     6                0     1           1         1           0        0     0       1     1
#    Dendritic cells       0         0     1                0     0           0         0           0        0     0       0     0
#    Erythrocytes          3         0     0                0     0           0         0           0        0     0       0     0
#    Granulocytes          0         1     0                0     0           0         0          11        0     0       0     0
#    Macrophages           1         0     0                0     0        4935         0           0        0     0       0     0
#    Monocytes             1         0   126                1     0          96       208           1        0     0       0     0
#    NK cells              0         0     0                0    10           0         0           0       39     0       0     0
#    T cells              32         0     1                0    27           0         1           0       16    38    1757    98

# label cells consistently called across both sets
sum(my.perit.singlets@meta.data$SingleR_RNAseq == my.perit.singlets@meta.data$SingleR_ImmGen) # 18805 (out of 19281)

# keep only cells with consistent labeling
my.perit.singlets.consistent       <- subset(my.perit.singlets, subset = SingleR_RNAseq == SingleR_ImmGen)    # 18805
save(my.perit.singlets.consistent, file = paste(Sys.Date(),"Seurat_10x_peritoneal_Singlets_SingleR_ConsistentID_SingleR.RData",sep = "_"))

# get proportions
my.freq.table <- prop.table(x = table(my.perit.singlets.consistent@meta.data$SingleR_RNAseq, my.perit.singlets.consistent@meta.data$Sample_ID), margin = 2)
my.freq.table <- my.freq.table[,c(5:6,1:2,7:8,3:4)]
my.freq.table.av <- apply(my.freq.table,1,mean)
my.freq.table.av.sort <- sort(my.freq.table.av, decreasing = T,index.return = T)

my.Bcell.freqs <- data.frame("Bcell_YF" = my.freq.table["B cells",1:2],
                             "Bcell_OF" = my.freq.table["B cells",3:4],
                             "Bcell_YM" = my.freq.table["B cells",5:6],
                             "Bcell_OM" = my.freq.table["B cells",7:8])

my.Mph.freqs <- data.frame("Mph_YF" = my.freq.table["Macrophages",1:2],
                           "Mph_OF" = my.freq.table["Macrophages",3:4],
                           "Mph_YM" = my.freq.table["Macrophages",5:6],
                           "Mph_OM" = my.freq.table["Macrophages",7:8])


my.Tcell.freqs <- data.frame("Tcell_YF" = my.freq.table["T cells",1:2],
                             "Tcell_OF" = my.freq.table["T cells",3:4],
                             "Tcell_YM" = my.freq.table["T cells",5:6],
                             "Tcell_OM" = my.freq.table["T cells",7:8])

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

my.perit.singlets.Mph.v3       <- subset(my.perit.singlets.consistent, subset = SingleR_ImmGen %in% "Macrophages")    # Macrophages, 4935 cells
my.perit.singlets.Bcells.v3    <- subset(my.perit.singlets.consistent, subset = SingleR_ImmGen %in% "B cells")        # B cells, 11866 cells


save(my.perit.singlets.Mph.v3       , file = paste(Sys.Date(),"Seurat_10x_peritoneal_Singlets_SingleR_Consistent_Macrophages.RData",sep = "_"))
save(my.perit.singlets.Bcells.v3    , file = paste(Sys.Date(),"Seurat_10x_peritoneal_Singlets_SingleR_Consistent_Bcells.RData",sep = "_"))

######################################
# Calculate and plot markers (on consistently called cells)
my.perit.singlets.consistent  <- SetIdent(my.perit.singlets.consistent , value = "SingleR_ImmGen")

#### get marker heatmap
perit.markers <- FindAllMarkers(my.perit.singlets.consistent, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10.perit.markers <- perit.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

png(paste0(Sys.Date(),"_SingleRconsistent_Peritoneal_aging_marker_heatmap.png"),width = 25, height = 18, units = "cm", res = 300)
DoHeatmap(my.perit.singlets.consistent, features = top10.perit.markers$gene)  + theme(axis.text.y = element_text(size = 8))
dev.off()

pdf(paste0(Sys.Date(),"SingleRconsistent_Peritoneal_aging_cluster_marker_Dtoplot.pdf"), width = 20, height = 8)
DotPlot(my.perit.singlets.consistent, features = unique(top10.perit.markers$gene), cols = c("deeppink1","deeppink4","deepskyblue1","deepskyblue4"), dot.scale = 8, split.by = "Condition") + RotatedAxis()
dev.off()

pdf(paste0(Sys.Date(),"SingleRconsistent_Peritoneal_aging_KNOWN_MARKERS_Dtoplot.pdf"), width = 10, height = 8)
DotPlot(my.perit.singlets.consistent, features = c("Adgre1","Itgam", "Bhlhe40", "Spi1","Cd79a","Cd79b","Cd19","Cd3d","Cd3e","Cd3g","Ccr2","Ly6c2"), cols = c("deeppink1","deeppink4","deepskyblue1","deepskyblue4"), dot.scale = 8, split.by = "Condition") + RotatedAxis()
dev.off()
############################################################################################################

#######################
sink(file = paste(Sys.Date(),"_scRNAseq_peritoneal_analysis_session_Info.txt", sep =""))
sessionInfo()
sink()
