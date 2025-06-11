setwd('/Volumes/BB_HQ_1/Immune_sex_dimorphism_Aging/Macrophages/Clean_Macrophage_Data/Peritoneal_scRNAseq/Combined_Analysis/SCENIC')
options(stringsAsFactors = F)

# load libraries for analysis
library('Seurat')
library(SCENIC)
library(AUCell)
library(RcisTarget)
library(SingleCellExperiment)


################################################################################################################################
#### 1. load and format dataset
# Load Seurat Object
load('../Cell_Proportions_AUGUR/2023-07-24_Seurat_COMBINED_10x_peritoneal_Singlets_SingleR_Consistent_Macrophages.RData')
All.perit.singlets.Mph
# An object of class Seurat 
# 52639 features across 15983 samples within 3 assays 
# Active assay: integrated (3000 features, 3000 variable features)
# 2 other assays present: RNA, SCT
# 2 dimensional reductions calculated: pca, umap

################################################################################################################################
#### 2. run SCENIC
# SCENIC workflow

# This tutorial goes through the steps in the **SCENIC workflow**:
# Building the **gene regulatory network (GRN)**: 
# 1. Identify potential targets for each TF based on co-expression.
# - Filtering the expression matrix and running GENIE3/GRNBoost. 
# - Formatting the targets from GENIE3/GRNBoost into co-expression modules. 

# 2.  Select potential direct-binding targets (regulons) based on DNA-motif analysis (*RcisTarget*: TF motif analysis) 
# Identify **cell states**:

# 3. Analyzing the network activity in each individual cell (*AUCell*)
# - Scoring regulons in the cells (calculate AUC)
# - Optional: Convert the network activity into ON/OFF (binary activity matrix)

# 4. Identify stable cell states based on their gene regulatory network activity (cell clustering) and exploring the results...

# Input

## Expression matrix
exprMat <- as.data.frame(All.perit.singlets.Mph@assays$SCT@data)

## Cell info
# In Step 3-4 (scoring the GRN and clustering), you can plot some information about the cells on the heatmaps or t-SNE. 
# You can choose which variables from the phenodata to plot, and assign them a specific color (otherwise one will be assigned automatically):
cellInfo       <- All.perit.singlets.Mph@meta.data[,c("Condition","Phase","source")]
cellInfo$nGene <- colSums(exprMat>0)
cellInfo       <- data.frame(cellInfo)
head(cellInfo)

saveRDS(cellInfo, file=paste(Sys.Date(),"cellInfo.Rds", sep = "_"))
saveRDS(cellInfo, file="int/cellInfo.Rds")

# Color to assign to the variables (same format as for NMF::aheatmap)
colVars <- list(sampleOrigin=setNames(c(rgb(253,154,251, 200, maxColorValue = 255),
                                        rgb(234, 62,122, 200, maxColorValue = 255),
                                        rgb(117,217,242, 200, maxColorValue = 255),
                                        rgb( 61, 93,154, 200, maxColorValue = 255)), 
                                      c("4m_F","20m_F","4m_M","20m_M")))
saveRDS(cellInfo, file=paste(Sys.Date(),"colVars.Rds", sep = "_"))
saveRDS(colVars, file="int/colVars.Rds")

pdf(paste(Sys.Date(),"color_scheme.pdf", sep = "_"))
plot.new(); legend(0,1, fill=colVars$sampleOrigin, legend=names(colVars$sampleOrigin))
dev.off()

# Initialize SCENIC settings
# to keep consistent settings across the multiple steps: common object where the options for the current run are stored
# using mm10 https://github.com/aertslab/SCENIC/issues/43
mm10_dbs <- list('500bp'= 'mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather', 
                 '10kb' = 'mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.feather')

scenicOptions <- initializeScenic(org          = "mgi"                 , # Mouse 
                                  dbDir        = "RcisTarget_databases", # RcisTarget databases location
                                  dbs          = mm10_dbs              , # mm10 databases 
                                  nCores       = 2                     )

scenicOptions@inputDatasetInfo$datasetTitle <- "SCENIC - PeriMac"
scenicOptions@settings$db_mcVersion <- "v9"


# Modify if needed (store parameters)
scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
scenicOptions@inputDatasetInfo$colVars  <- "int/colVars.Rds"

# Save to use at a later time
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 

# Co-expression network
# The first step on SCENIC workflow is to infer potential transcription factor targets based on the expression data. 

# *Subsampling*: When there is a high proportion of low-quality cells, or if the computation time is an issue, 
# it is also possible to infer the regulatory network using a subset of cells (e.g. selecting random or high-quality cells as input to the co-expression analysis). 
# The activity of the regulatory network, *trained* on this subset of cells, can then be evaluated on all the cells in the dataset with AUCell (Step 3).
# Note that to avoid loss of resolution, the subset of cells should be representative of the whole dataset (e.g. contain sufficient representation of all the cell types).

## Gene filter/selection
# To run GENIE3/GRNBoost we recommend to apply soft gene filter, to remove genes that are expressed either at very low levels or in too few cells. 
# Apply a filtering based on the total number of counts of the gene, and the number of cells in which it is detected. 
# Feel free to modify the filters according to your own needs/dataset. 

# The first filter, the total number of reads per gene is meant to remove genes that are most likely noise. 
# The second filter, the number of cells in which the gene is detected is to avoid that such genes gain a lot of weight if they coincide.
# To proceed with these filtering, we will first calculate some basic statistics on the expression matrix:
# Number of cells in which each gene is expressed, and number of counts (in the dataset unit) per gene:
nCellsPerGene  <- apply(exprMat, 1, function(x) sum(x>0))
nCountsPerGene <- apply(exprMat, 1, sum)

summary(nCellsPerGene)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0      32     711    2128    2859   15983 

summary(nCountsPerGene)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#     0.0    29.5   599.9  2488.2  2446.4 92284.6 

max(exprMat) # [1] 9.039787
sum(exprMat>0) / sum(exprMat==0) # 0.1536257

# **First filter:** Keep only the genes with at least `r 3*.01*ncol(exprMat)` UMI counts across all samples 
# Here the totag. thg. thg. thne woul tota tota totaas expre number the gene would have, i number the gene would have, ire more homogeneous) it was e
minReads           <- 2*0.10*ncol(exprMat)
genesLeft_minReads <- names(nCountsPerGene)[which(nCountsPerGene > minReads)]
length(genesLeft_minReads) # 3393

# **Second filter:** Keep the genes that are detected in at least 20% of the cells. (macrophages are more homogeneous)
minSamples         <- ncol(exprMat)*0.20
nCellsPerGene2     <- nCellsPerGene[genesLeft_minReads]
genesLeft_minCells <- names(nCellsPerGene2)[which(nCellsPerGene2 > minSamples)]
length(genesLeft_minCells) # 3376

# **Genes in databases:**
# In upcoming steps (e.g. motif enrichment analysis), only the genes that are available on RcisTarget databases will be used.
# To save some running time for GENIE3/GRNBoost, we can ignore the genes that are not in the databases.
# Load corresponding databases:
motifRankings   <- importRankings(getDatabases(scenicOptions)[[1]]) # either one, they should have the same genes
genesInDatabase <- colnames(getRanking(motifRankings))

# Exclude missing genes:
genesLeft_minCells_inDatabases <- genesLeft_minCells[which(genesLeft_minCells %in% genesInDatabase)]
length(genesLeft_minCells_inDatabases) # 3282

# We can now **filter the expression matrix** to contain only these `r length(genesLeft_minCells_inDatabases)` genes. 
# This matrix is now ready for the co-expression analysis.
genesKept <- genesLeft_minCells_inDatabases
saveRDS(genesKept, file=getIntName(scenicOptions, "genesKept"))
exprMat_filtered <- as.matrix(exprMat[genesKept, ])

# To avoid confusions in the following steps: 
rm(exprMat)

## Correlation
# GENIE3/GRNBoost can detect both positive and negative associations. In order to distinguish potential activation from repression, 
# we will split the targets into positive- and negative-correlated targets (i.e. Spearman correlation between the TF and the potential target).

# Calculate the correlation. This step can be run either before/after or simultaneously to GENIE3/GRNBoost)*
corrMat <- cor(t(exprMat_filtered), method="spearman")

# Only the rows for TFs will be needed needed
allTFs  <- getDbTfs(scenicOptions)
corrMat <- corrMat[which(rownames(corrMat) %in% allTFs),]
saveRDS(corrMat, file=getIntName(scenicOptions, "corrMat"))

## Run GENIE3
# The input to GENIE3 is typically an expression matrix, and a list of candidate regulators. 
# GENIE3 will typically take several hours (or days) to run. If you are running this workflow on an RStudio session, we recommend that you stop here and run the next code chunk in an independent R console (i.e. with `screen`/`tmux`) or in an server/HPC (if available). 
# The upcoming code chunks will resume the workflow by loading GENIE3 output.

# We are using normalized SCT data, so no need to add log

# Run GENIE3
set.seed(123456789) # for reproducibility
runGenie3(exprMat_filtered      , 
          scenicOptions         )
# Using 272 TFs as potential regulators...
# Running GENIE3 part 1
# Running GENIE3 part 2
# Running GENIE3 part 3
# Running GENIE3 part 4
# Running GENIE3 part 5
# Running GENIE3 part 6
# Running GENIE3 part 7
# Running GENIE3 part 8
# Running GENIE3 part 9
# Running GENIE3 part 10
# Finished running GENIE3.
# Warning message:
#   In runGenie3(exprMat_filtered, scenicOptions) :
#   Only 16% of the 1721 TFs in the database were found in the dataset. Do they use the same gene IDs?


# Build and score the GRN (runSCENIC)
# Once the results from GENIE3/GRNBoost (and the correlation) are ready, the remaining steps of SCENIC can be run. 
# 
# The easiest/fastest way is to use the following *wrapper* functions, each of them corresponding to one of the main steps in SCENIC workflow:
#   
# Build the *gene regulatory network*: 
#   1. Get co-expression modules
#   2. Get regulons (with `r  Biocpkg("RcisTarget")`): TF motif analysis)
# 
# Identify *cell states*:
#   3. Score GRN (regulons) in the cells (with `r  Biocpkg("AUCell")`)
#   4. Cluster cells according to the GRN activity

# Run the remaining steps using the *wrapper* functions: 
scenicOptions <- readRDS("int/scenicOptions.Rds")
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores  <- 2
scenicOptions@settings$seed    <- 123456789

runSCENIC_1_coexNetwork2modules(scenicOptions)
# 15:38	Creating TF modules
# 75%         90% 
#   0.004130292 0.005778739 
# Number of links between TFs and targets: 888197
# [,1]
# nTFs          272
# nTargets     3282
# nGeneSets    1395
# nLinks    1250710

runSCENIC_2_createRegulons(scenicOptions)
# 19:28	Step 2. Identifying regulons
# tfModulesSummary:
#  top5perTarget top10perTarget           w005 top50perTarget          top50           w001 
#             62            101            173            183            271            272 
# 19:28	RcisTarget: Calculating AUC
# Scoring database:  [Source file: mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather]
# Scoring database:  [Source file: mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.feather]
# 20:05	RcisTarget: Adding motif annotation
# Number of motifs in the initial enrichment: 758505
# Number of motifs annotated to the matching TF: 9095
# 20:08	RcisTarget: Prunning targets
# 20:45	Number of motifs that support the regulons: 9095
# [WARNING] Deprecated: --self-contained. use --embed-resources --standalone
# Preview of motif enrichment saved as: output/Step2_MotifEnrichment_preview.html
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2.00   19.25   85.00  473.92  868.25 2205.00 

runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered)
# 20:46	Step 3. Analyzing the network activity in each individual cell
# Number of regulons to evaluate on cells: 127
# Biggest (non-extended) regulons: 
#   Elf2 (1780g)
# Ep300 (1363g)
# Etv1 (1221g)
# Chd2 (1016g)
# Yy1 (964g)
# Elf1 (858g)
# Xbp1 (744g)
# Fli1 (713g)
# Nr3c1 (670g)
# Bclaf1 (420g)
# Quantiles for the number of genes detected by cell: 
#   (Non-detected genes are shuffled at the end of the ranking. Keep it in mind when choosing the threshold for calculating the AUC).
# min      1%      5%     10%     50%    100% 
# 384.00  473.82  691.10  944.00 1676.00 2450.00 
# Using 2 cores.
# Using 2 cores.
# 20:52	Finished running AUCell.
# 20:52	Plotting heatmap...
# 20:52	Plotting t-SNEs...
saveRDS(scenicOptions, file="int/scenicOptions.Rds") # To save status


# Clustering / dimensionality reduction on the regulon activity
# The cells can be grouped/clustered based on the regulon activity, either continuous or binarized .
# If using t-SNE as visualization, it is recommended to try different settings to evaluate the stability of the states/clusters. 
# Feel free to use UMAP, other clustering methods (or trajectory inference methods, if appropriate) instead.
# 
# The function included in SCENIC package runs multiple t-SNEs with different settings; 
# It will create all combinations between the selected “number of PCs” and “perplexity” 
# (expected running time: few minutes to hours, depending on the number of cells):
nPcs <- c(5,10,15,30)
scenicOptions@settings$seed <- 123456789 # same seed for all of them

# Run t-SNE with different settings:
fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=nPcs, perpl=c(5,15,50))
fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=nPcs, perpl=c(5,15,50), onlyHighConf=TRUE, filePrefix="int/tSNE_oHC")

# Plot as pdf (individual files in int/):
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_", list.files("int"), value=T), value=T))

# and to view/compare them…
par(mfrow=c(length(nPcs), 3))
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_AUC", list.files("int"), value=T, perl = T), value=T))
plotTsne_compareSettings(fileNames, scenicOptions, showLegend=FALSE, varName="Condition", cex=.5)

# Using only "high-confidence" regulons (normally similar)
par(mfrow=c(3,3))
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_oHC_AUC", list.files("int"), value=T, perl = T), value=T))
plotTsne_compareSettings(fileNames, scenicOptions, showLegend=FALSE, varName="Condition", cex=.5)
# The chosen t-SNE can then be saved as default to use for plots (can also be “binary”, see below):
par(mfrow=c(1,1))

scenicOptions@settings$defaultTsne$aucType <- "AUC"
scenicOptions@settings$defaultTsne$dims    <- 30
scenicOptions@settings$defaultTsne$perpl   <- 15
saveRDS(scenicOptions, file="int/scenicOptions.Rds")



# Cell states
# AUCell provides the activity of the regulons across the cells. By clustering the cells based on this regulon activity 
# (either the continuous or binary AUC matrix), we can see whether there are groups of cells that tend to have the same 
# regulons active, and reveal the network states that are recurrent across multiple cells. These states would be equivalent
# to the attractor states of the network. Combining these clustering with different visualization methods, we can explore 
# the association of cell states with specific regulons.
# SCENIC provides some wrapper functions to get a quick overview. For example, projecting the AUC and TF expression onto t-SNEs, 
# and visualizing of the AUC as heatmaps, but feel free to explore alternative clustering and visualization tools.
# 
# Projection the AUC and TF expression onto t-SNEs
# Briefly, a t-SNE is a 2D projection of the cells, where cells (dots) are placed close to each other if they have similar input profiles (in our case, regulon activity). The t-SNE usually allows to get a quick and easy overview of the cell states in the dataset. Note however, that t-SNE works well to identify distinct classes, but it is not appropiate for dinamic/continuous processes (e.g. trajectory-like visualizations).
exprMat         <- as.data.frame(All.perit.singlets.Mph@assays$SCT@data)
aucellApp       <- plotTsne_AUCellApp(scenicOptions, exprMat) # default t-SNE
savedSelections <- shiny::runApp(aucellApp)

# AUCell_plotTSNE() to save static plots:
print(tsneFileName(scenicOptions))


# Regulators for known cell types or clusters
# The regulatory analysis from SCENIC can be combined with other analyses (typically clustering), or focus on regulators for specific cell types. There are multiple options to do these analyses (your imagination is the limit!). Here are some quick examples to start:
#   
# Average Regulon Activity by cluster
# To start from clusters/cell types from Seurat: cellInfo <- data.frame(seuratCluster=Idents(seuratObject)))

regulonAUC                     <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC                     <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byGroup        <- sapply(split(rownames(cellInfo), cellInfo$Condition),
                                         function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byGroup_Scaled <- t(scale(t(regulonActivity_byGroup), center = T, scale=T))

save.image(paste0(Sys.Date(),"_SCENIC_Analysis_workspace.RData"))
####################################################################################################################################


####################################################################################################################################
load('2023-07-28_SCENIC_Analysis_workspace.RData')

library("ComplexHeatmap")
library("grid")
library(circlize)

pdf(paste0(Sys.Date(),"HEATMAP_of_top_SCENIC_Regulons_Macrophage_scRNAseq_Aging.pdf"), height = 11, width = 7)
Heatmap(regulonActivity_byGroup_Scaled[,c("YF", "OF","YM", "OM")],
        name="SCENIC Regulon activity",
        cluster_columns = F,
        border = T, rect_gp = gpar(col = "grey", lwd = 0.5),
        column_title = "SCENIC Regulons (Macrophages)")
dev.off()

# topRegulators <- reshape2::melt(regulonActivity_byGroup_Scaled)
# colnames(topRegulators) <- c("Regulon", "Condition", "RelativeActivity")
# topRegulators <- topRegulators[which(topRegulators$RelativeActivity>0),]
# viewTable(topRegulators)

# save.image("/Volumes/BB_Home_HQ/Immune_sex_dimorphism_Aging/Macrophages/Clean_Macrophage_Data/Peritoneal_scRNAseq/SCENIC/SCENIC_Analysis_workspace.RData")


cellInfo       <- All.perit.singlets.Mph@meta.data[,c("Condition","Phase","source","Sample_ID")]
cellInfo$nGene <- colSums(exprMat>0)
cellInfo       <- data.frame(cellInfo)
head(cellInfo)
#                             Condition Phase source Sample_ID nGene
# perit_YM1_AAAGATGAGAATGTGT-1        YM    G1 NIA_v2 YMCohort1  1596
# perit_YM1_AAAGATGCAGATCTGT-1        YM    G1 NIA_v2 YMCohort1  1578
# perit_YM1_AAAGTAGAGTGTACGG-1        YM     S NIA_v2 YMCohort1   971
# perit_YM1_AAAGTAGTCTTCGAGA-1        YM     S NIA_v2 YMCohort1  1182
# perit_YM1_AACACGTTCATCTGTT-1        YM     S NIA_v2 YMCohort1  1538
# perit_YM1_AACGTTGAGTGGGATC-1        YM    G1 NIA_v2 YMCohort1  1652


regulonActivity_bySample        <- sapply(split(rownames(cellInfo), cellInfo$Sample_ID),
                                          function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
colnames(regulonActivity_bySample)

regulonActivity_bySample_NIAv2   <- regulonActivity_bySample[,c("YFCohort1","YFCohort2","OFCohort1","OFCohort2","YMCohort1","YMCohort2","OMCohort1","OMCohort2")]
regulonActivity_bySample_NIAv3   <- regulonActivity_bySample[,c("NIA_YF1","NIA_YF2","NIA_YF3","NIA_OF1","NIA_OF2","NIA_YM1","NIA_YM2","NIA_YM3","NIA_OM1","NIA_OM2")]
regulonActivity_bySample_JAX     <- regulonActivity_bySample[,c("JAX_YF1","JAX_YF2","JAX_YF3", "JAX_OF1","JAX_OF2","JAX_OF3","JAX_YM1","JAX_YM2","JAX_YM3","JAX_OM1","JAX_OM2","JAX_OM3" )]

regulonActivity_bySample_both  <- regulonActivity_bySample[,c("JAX_YF1","JAX_YF2","JAX_YF3","YFCohort1","YFCohort2","NIA_YF1","NIA_YF2","NIA_YF3",
                                                              "JAX_OF1","JAX_OF2","JAX_OF3","OFCohort1","OFCohort2","NIA_OF1","NIA_OF2",
                                                              "JAX_YM1","JAX_YM2","JAX_YM3","YMCohort1","YMCohort2","NIA_YM1","NIA_YM2","NIA_YM3",
                                                              "JAX_OM1","JAX_OM2","JAX_OM3","OMCohort1","OMCohort2","NIA_OM1","NIA_OM2")]

regulonActivity_bySample_NIAv2_Scaled <- t(scale(t(regulonActivity_bySample_NIAv2), center = T, scale=T))
regulonActivity_bySample_NIAv3_Scaled <- t(scale(t(regulonActivity_bySample_NIAv3), center = T, scale=T))
regulonActivity_bySample_JAX_Scaled   <- t(scale(t(regulonActivity_bySample_JAX)  , center = T, scale=T))
regulonActivity_bySample_both_Scaled  <- t(scale(t(regulonActivity_bySample_both) , center = T, scale=T))


pdf(paste0(Sys.Date(),"HEATMAP_of_top_SCENIC_Regulons_Macrophage_scRNAseq_JAX_by_Animal.pdf"), height = 11, width = 8)
Heatmap(regulonActivity_bySample_JAX_Scaled,
        name="Regulon activity",
        cluster_columns = F,
        border = T, rect_gp = gpar(col = "grey", lwd = 0.5),
        column_title = "SCENIC Regulons (Mph, JAX)")
dev.off()


pdf(paste0(Sys.Date(),"HEATMAP_of_top_SCENIC_Regulons_Macrophage_scRNAseq_NIAv2_by_Pool.pdf"), height = 11, width = 8)
Heatmap(regulonActivity_bySample_NIAv2_Scaled,
        name="Regulon activity",
        cluster_columns = F,
        border = T, rect_gp = gpar(col = "grey", lwd = 0.5),
        column_title = "SCENIC Regulons (Mph, NIA)")
dev.off()

pdf(paste0(Sys.Date(),"HEATMAP_of_top_SCENIC_Regulons_Macrophage_scRNAseq_NIAv3_by_Pool.pdf"), height = 11, width = 8)
Heatmap(regulonActivity_bySample_NIAv3_Scaled,
        name="Regulon activity",
        cluster_columns = F,
        border = T, rect_gp = gpar(col = "grey", lwd = 0.5),
        column_title = "SCENIC Regulons (Mph, NIA)")
dev.off()

pdf(paste0(Sys.Date(),"HEATMAP_of_top_SCENIC_Regulons_Macrophage_scRNAseq_BOTH_by_library.pdf"), height = 11, width = 8)
Heatmap(regulonActivity_bySample_both_Scaled,
        name="Regulon activity",
        cluster_columns = F,
        border = T, rect_gp = gpar(col = "grey", lwd = 0.5),
        column_title = "SCENIC Regulons (Mph, both)")
dev.off()



# Group-level for NIA and JAX separately
cellInfo.JAX   <- cellInfo[cellInfo$source == "JAX_v3",]
cellInfo.NIAv2 <- cellInfo[cellInfo$source == "NIA_v2",]
cellInfo.NIAv3 <- cellInfo[cellInfo$source == "NIA_v3",]


regulonActivity_byGroup.JAX        <- sapply(split(rownames(cellInfo.JAX), cellInfo.JAX$Condition),
                                             function(cells) rowMeans(getAUC(regulonAUC)[,cells]))[,c("YF", "OF", "YM", "OM")]
regulonActivity_byGroup_Scaled.JAX  <- t(scale(t(regulonActivity_byGroup.JAX), center = T, scale=T))

pdf(paste0(Sys.Date(),"HEATMAP_of_top_SCENIC_Regulons_Macrophage_scRNAseq_JAX_by_GROUP.pdf"), height = 11, width = 8)
Heatmap(regulonActivity_byGroup_Scaled.JAX,
        name="Regulon activity",
        cluster_columns = F,
        border = T, rect_gp = gpar(col = "grey", lwd = 0.5),
        column_title = "SCENIC Regulons (Mph, JAX)")
dev.off()



regulonActivity_byGroup.NIAv2        <- sapply(split(rownames(cellInfo.NIAv2), cellInfo.NIAv2$Condition),
                                             function(cells) rowMeans(getAUC(regulonAUC)[,cells]))[,c("YF", "OF", "YM", "OM")]
regulonActivity_byGroup_Scaled.NIAv2  <- t(scale(t(regulonActivity_byGroup.NIAv2), center = T, scale=T))

pdf(paste0(Sys.Date(),"HEATMAP_of_top_SCENIC_Regulons_Macrophage_scRNAseq_NIAv2_by_GROUP.pdf"), height = 11, width = 8)
Heatmap(regulonActivity_byGroup_Scaled.NIAv2,
        name="Regulon activity",
        cluster_columns = F,
        border = T, rect_gp = gpar(col = "grey", lwd = 0.5),
        column_title = "SCENIC Regulons (Mph, NIAv2)")
dev.off()



regulonActivity_byGroup.NIAv3        <- sapply(split(rownames(cellInfo.NIAv3), cellInfo.NIAv3$Condition),
                                             function(cells) rowMeans(getAUC(regulonAUC)[,cells]))[,c("YF", "OF", "YM", "OM")]
regulonActivity_byGroup_Scaled.NIAv3  <- t(scale(t(regulonActivity_byGroup.NIAv3), center = T, scale=T))

pdf(paste0(Sys.Date(),"HEATMAP_of_top_SCENIC_Regulons_Macrophage_scRNAseq_NIAv3_by_GROUP.pdf"), height = 11, width = 8)
Heatmap(regulonActivity_byGroup_Scaled.NIAv3,
        name="Regulon activity",
        cluster_columns = F,
        border = T, rect_gp = gpar(col = "grey", lwd = 0.5),
        column_title = "SCENIC Regulons (Mph, NIAv3)")
dev.off()


############################################################################################################
############################################################################################################
## Add AUC info to Seurat Object
auc.regulon <- getAUC(regulonAUC)

sum(rownames(All.perit.singlets.Mph@meta.data) == colnames(auc.regulon))# 15983
All.perit.singlets.Mph@meta.data <- cbind(All.perit.singlets.Mph@meta.data,t(auc.regulon))
All.perit.singlets.Mph@meta.data$Condition2 <- factor(All.perit.singlets.Mph@meta.data$Condition, levels = c("YF", "OF", "YM", "OM"))

YFcells <- rownames(All.perit.singlets.Mph@meta.data)[All.perit.singlets.Mph@meta.data$Condition == "YF"]
OFcells <- rownames(All.perit.singlets.Mph@meta.data)[All.perit.singlets.Mph@meta.data$Condition == "OF"]
YMcells <- rownames(All.perit.singlets.Mph@meta.data)[All.perit.singlets.Mph@meta.data$Condition == "YM"]
OMcells <- rownames(All.perit.singlets.Mph@meta.data)[All.perit.singlets.Mph@meta.data$Condition == "OM"]


# Plot regulons of interest
pdf(paste0(Sys.Date(),"Violin_Irf2_SCENIC_AUC_Mph.pdf"), height = 6, width = 9)
VlnPlot(object = All.perit.singlets.Mph,
        features = c("Irf2 (12g)" ),
        ncol = 2,
        group.by = "Condition2",
        col = c("deeppink","deeppink4","deepskyblue","deepskyblue4"),
        pt.size = 0.25)
dev.off()

irf2.f <- wilcox.test(All.perit.singlets.Mph@meta.data[YFcells,]$`Irf2 (12g)`,All.perit.singlets.Mph@meta.data[OFcells,]$`Irf2 (12g)`)
irf2.m <- wilcox.test(All.perit.singlets.Mph@meta.data[YMcells,]$`Irf2 (12g)`,All.perit.singlets.Mph@meta.data[OMcells,]$`Irf2 (12g)`)
irf2.f$p.value # 1.384503e-35
irf2.m$p.value # 3.013703e-22


# Plot regulons of interest
pdf(paste0(Sys.Date(),"Violin_Tbl1rx1_SCENIC_AUC_Mph.pdf"), height = 6, width = 9)
VlnPlot(object = All.perit.singlets.Mph,
        features = c("Tbl1xr1_extended (306g)" ),
        ncol = 2,
        group.by = "Condition2",
        col = c("deeppink","deeppink4","deepskyblue","deepskyblue4"),
        pt.size = 0.25)
dev.off()

Tbl1xr1.f <- wilcox.test(All.perit.singlets.Mph@meta.data[YFcells,]$`Tbl1xr1_extended (306g)`,All.perit.singlets.Mph@meta.data[OFcells,]$`Tbl1xr1_extended (306g)`)
Tbl1xr1.m <- wilcox.test(All.perit.singlets.Mph@meta.data[YMcells,]$`Tbl1xr1_extended (306g)`,All.perit.singlets.Mph@meta.data[OMcells,]$`Tbl1xr1_extended (306g)`)
Tbl1xr1.f$p.value # 2.638548e-49
Tbl1xr1.m$p.value # 3.882084e-19


# Plot regulons of interest
pdf(paste0(Sys.Date(),"Violin_Mef2c_SCENIC_AUC_Mph.pdf"), height = 6, width = 9)
VlnPlot(object = All.perit.singlets.Mph,
        features = c("Mef2c_extended (54g)" ),
        ncol = 2,
        group.by = "Condition2",
        col = c("deeppink","deeppink4","deepskyblue","deepskyblue4"),
        pt.size = 0.25)
dev.off()

Mef2c.f <- wilcox.test(All.perit.singlets.Mph@meta.data[YFcells,]$`Mef2c_extended (54g)`,All.perit.singlets.Mph@meta.data[OFcells,]$`Mef2c_extended (54g)`)
Mef2c.m <- wilcox.test(All.perit.singlets.Mph@meta.data[YMcells,]$`Mef2c_extended (54g)`,All.perit.singlets.Mph@meta.data[OMcells,]$`Mef2c_extended (54g)`)
Mef2c.f$p.value # 0.0003214881
Mef2c.m$p.value # 1.574523e-05



library(ggplot2)
pdf(paste0(Sys.Date(),"Top_TF_Bulk_with_significant_SCENIC_AUC_Mph_DOTPLOT.pdf"), height = 4, width = 7)
DotPlot(All.perit.singlets.Mph,
        features = rev(c("Irf1 (37g)","Irf2 (12g)","Rela_extended (1841g)","Stat1 (98g)","Stat3 (16g)","Mef2c_extended (54g)","Tbl1xr1_extended (306g)")),
        group.by = "Condition2", cols = c("grey", "red"), col.min = -1,col.max = 1)  + coord_flip()
dev.off()

pdf(paste0(Sys.Date(),"Top_TF_Bulk_F_ONLY_with_significant_SCENIC_AUC_Mph_DOTPLOT.pdf"), height = 3, width = 7)
DotPlot(All.perit.singlets.Mph,
        features = rev(c("Irf2 (12g)","Mef2c_extended (54g)","Tbl1xr1_extended (306g)")),
        group.by = "Condition2", cols = c("grey", "red"), col.min = -1,col.max = 1)   + coord_flip() + ggtitle("SCENIC Regulon Activity") + scale_size_area(limits = c(50,100))
dev.off()


##### TF/TRs
# parse names for comparison
my.parse.1 <- unlist(lapply(strsplit(rownames(regulonActivity_byGroup_Scaled)," "),'[[',1))
my.parse.2 <- unlist(lapply(strsplit(my.parse.1,"_extended"),'[[',1))

# which TFs have aging regulation (use scaled data)
reg_act.aging         <- data.frame(regulonActivity_byGroup_Scaled)
reg_act.aging         <- reg_act.aging[,c("YF", "OF", "YM", "OM")]
reg_act.aging$TF_Name <- my.parse.2

reg_act.aging$F_aging_p <- NA
reg_act.aging$M_aging_p <- NA

# Run stats for female and male changes
for (i in 1:nrow(regulonActivity_byGroup_Scaled)){
  my.regulon <- rownames(regulonActivity_byGroup_Scaled)[i]

  reg_act.aging$F_aging_p[i] <- wilcox.test(All.perit.singlets.Mph@meta.data[YFcells,my.regulon],All.perit.singlets.Mph@meta.data[OFcells,my.regulon])$p.value
  reg_act.aging$M_aging_p[i] <- wilcox.test(All.perit.singlets.Mph@meta.data[YMcells,my.regulon],All.perit.singlets.Mph@meta.data[OMcells,my.regulon])$p.value

}

# Multiple test correction
reg_act.aging$F_aging_FDR <- p.adjust(reg_act.aging$F_aging_p, method = "BH")
reg_act.aging$M_aging_FDR <- p.adjust(reg_act.aging$M_aging_p, method = "BH")

reg_act.aging$F_aging_FDR5 <- ifelse(reg_act.aging$F_aging_FDR < 0.05, "YES", "NO")
reg_act.aging$M_aging_FDR5 <- ifelse(reg_act.aging$M_aging_FDR < 0.05, "YES", "NO")

reg_act.aging$F_aging_FDR1e5 <- ifelse(reg_act.aging$F_aging_FDR < 1e-5, "YES", "NO")
reg_act.aging$M_aging_FDR1e5 <- ifelse(reg_act.aging$M_aging_FDR < 1e-5, "YES", "NO")

write.table(reg_act.aging, file = paste0(Sys.Date(),"_SCENIC_REGULONS_Wilcox_test.txt"), sep = "\t", quote = F)

# Save Seurat
save(All.perit.singlets.Mph, file = paste0(Sys.Date(),"_Seurat_Annotated_by_AUC_SCENIC.RData"))

##############################################################################################################
##############################################################################################################
load('2023-07-28_Seurat_Annotated_by_AUC_SCENIC.RData')
library(ggplot2)
library(ggpubr)

####
# Original combined
pdf(paste0(Sys.Date(),"Top_TF_Bulk_F_ONLY_with_significant_SCENIC_AUC_Mph_DOTPLOT_COMBINED.pdf"), height = 3, width = 7)
DotPlot(All.perit.singlets.Mph,
        features = rev(c("Irf2 (12g)","Mef2c_extended (54g)","Tbl1xr1_extended (306g)")),
        group.by = "Condition2", cols = c("grey", "red"), col.min = -1,col.max = 1)   + coord_flip() + ggtitle("SCENIC Regulon Activity") + scale_size_area(limits = c(50,100))
dev.off()

# split by source
mph.obj.list <- SplitObject(All.perit.singlets.Mph, split.by = "source")
names(mph.obj.list)
# [1] "NIA_v2" "JAX_v3" "NIA_v3"

NIA_v2 <- DotPlot(mph.obj.list$NIA_v2, features = rev(c("Irf2 (12g)","Mef2c_extended (54g)","Tbl1xr1_extended (306g)")), group.by = "Condition", cols = c("grey", "red"), col.min = -1,col.max = 1) + coord_flip() + ggtitle("NIA_v2") + scale_size_area(limits = c(50,100))
NIA_v3 <- DotPlot(mph.obj.list$NIA_v3, features = rev(c("Irf2 (12g)","Mef2c_extended (54g)","Tbl1xr1_extended (306g)")), group.by = "Condition", cols = c("grey", "red"), col.min = -1,col.max = 1) + coord_flip() + ggtitle("NIA_v3") + scale_size_area(limits = c(50,100))
JAX_v3 <- DotPlot(mph.obj.list$JAX_v3, features = rev(c("Irf2 (12g)","Mef2c_extended (54g)","Tbl1xr1_extended (306g)")), group.by = "Condition", cols = c("grey", "red"), col.min = -1,col.max = 1) + coord_flip() + ggtitle("JAX_v3") + scale_size_area(limits = c(50,100))

p.3c   <- ggarrange(NIA_v2, NIA_v3, JAX_v3, ncol = 3, nrow = 1,  common.legend = T, legend = "right")
p.3c

pdf(paste0(Sys.Date(),"_Top_TF_Bulk_F_ONLY_with_significant_SCENIC_AUC_Mph_DOTPLOT_3_cohorts.pdf"), height = 3, width = 12)
plot(p.3c)
dev.off()


# Plot regulons of interest
pdf(paste0(Sys.Date(),"Violin_Irf2_SCENIC_AUC_Mph_COMBINED.pdf"), height = 4.5, width = 8)
VlnPlot(object = All.perit.singlets.Mph,
        features = c("Irf2 (12g)" ),
        ncol = 2,
        group.by = "Condition2",
        col = c("deeppink","deeppink4","deepskyblue","deepskyblue4"),
        pt.size = 0.25)
dev.off()


NIA_v2 <- VlnPlot(mph.obj.list$NIA_v2, features = c("Irf2 (12g)" ), ncol = 2, group.by = "Condition2", col = c("deeppink","deeppink4","deepskyblue","deepskyblue4"), pt.size = 0.25) + ggtitle("NIA_v2")
NIA_v3 <- VlnPlot(mph.obj.list$NIA_v3, features = c("Irf2 (12g)" ), ncol = 2, group.by = "Condition2", col = c("deeppink","deeppink4","deepskyblue","deepskyblue4"), pt.size = 0.25) + ggtitle("NIA_v3")
JAX_v3 <- VlnPlot(mph.obj.list$JAX_v3, features = c("Irf2 (12g)" ), ncol = 2, group.by = "Condition2", col = c("deeppink","deeppink4","deepskyblue","deepskyblue4"), pt.size = 0.25) + ggtitle("JAX_v3")

p.3c   <- ggarrange(NIA_v2, NIA_v3, JAX_v3, ncol = 3, nrow = 1,  common.legend = T, legend = "right")
p.3c


pdf(paste0(Sys.Date(),"_Violin_Irf2_SCENIC_AUC_Mph_COMBINED_3_cohorts.pdf"), height = 4.5, width = 15)
plot(p.3c)
dev.off()

#######################
sink(file = paste(Sys.Date(),"_Macrophage_scRNAseq_SCENIC_analysis_session_Info.txt", sep =""))
sessionInfo()
sink()


