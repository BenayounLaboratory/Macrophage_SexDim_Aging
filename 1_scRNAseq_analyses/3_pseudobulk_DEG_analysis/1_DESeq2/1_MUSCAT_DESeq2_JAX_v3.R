setwd('/Volumes/BB_HQ_1/Immune_sex_dimorphism_Aging/Macrophages/Clean_Macrophage_Data/Peritoneal_scRNAseq/v3_set/Muscat_PB')
options(stringsAsFactors=FALSE)

library(Seurat)
library(beeswarm)
library(ggplot2)

library(muscat)
library(sva)
library(limma)
library(DESeq2)
library(bitops)

library(ggplot2)          # 
library(scales)           # 
library(ComplexHeatmap)   #
library(circlize)         #

theme_set(theme_bw())   

# 2025-08-28
# Run Muscat DE on all QC cell types (JAX v3 set)

###############################################################################################
# 0. preprocess Seurat object for use with muscat

load('../Seurat/2023-02-22_Seurat_JAX_10x_peritoneal_Singlets_SingleR_ConsistentID_SingleR.RData')
JAX.singlets.consistent
# An object of class Seurat 
# 49037 features across 14595 samples within 3 assays 
# Active assay: SCT (16740 features, 2000 variable features)
#  2 other assays present: RNA, HTO
#  2 dimensional reductions calculated: pca, umap


# bring RNA as main assay for processing
DefaultAssay(JAX.singlets.consistent) <- "RNA"

#### Make covariate table
perit.meta <- unique(JAX.singlets.consistent@meta.data[,c("Sample_ID", "Condition")])
rownames(perit.meta) <- perit.meta$Sample_ID
write.table(perit.meta, file = paste0(Sys.Date(),"_JAXv3_sample_metadata_table.txt"), sep = "\t", row.names = F, quote = F)

# convert to SingleCellExperiment
# https://satijalab.org/seurat/archive/v3.1/conversion_vignette.html
JAX.singlets.consistent.sce <- as.SingleCellExperiment(JAX.singlets.consistent)
save(JAX.singlets.consistent.sce, file = paste(Sys.Date(),"SingleCellExperimnents_peritoneal_JAXv3.RData",sep = "_"))

rm(JAX.singlets.consistent) # free up some memory
###############################################################################################

###############################################################################################
# 1. Run muscat for pseudobulking and extraction of samples

###############################################
#######        Data preparation         #######
###############################################
perit.sce.cl <- prepSCE(JAX.singlets.consistent.sce, 
                        cluster_id    = "SingleR_ImmGen",  # population assignments
                        group_id      = "Condition"     ,  # group IDs (ctrl/stim)
                        sample_id     = "Sample_ID"     ,  # sample IDs (ctrl/stim.1234)
                        drop          = TRUE            )   # drop all other colData columns

# store cluster and sample IDs, as well as the number of clusters and samples into the following simple variables:
nk  <- length(kids <- levels(perit.sce.cl$cluster_id))
ns  <- length(sids <- levels(perit.sce.cl$sample_id))
names(kids) <- kids; names(sids) <- sids

# nb. of cells per cluster-sample
t(table(perit.sce.cl$cluster_id, perit.sce.cl$sample_id))
#        B cells Macrophages Monocytes NK cells T cells
# JAX_OF1     839         306        15        1      60
# JAX_OF2     645         383         8        3      42
# JAX_OF3    1045         219         6        1      48
# JAX_OM1     715         302        11        2     124
# JAX_OM2     687         342         7        2     119
# JAX_OM3     690         349         6        2     186
# JAX_YF1     491         496        11        2      42
# JAX_YF2     707         574        17        6      78
# JAX_YF3     690         547        21        2      43
# JAX_YM1     436         522        22        7      39
# JAX_YM2     455         836        18        2      52
# JAX_YM3     599         599        24        9      83

# Aggregation of single-cell to pseudobulk data
pb.perit <- aggregateData(perit.sce.cl, assay = "counts", fun = "sum", by = c("cluster_id", "sample_id"))

# one list item per cell type
assayNames(pb.perit)
# [1] "B cells"     "Macrophages" "Monocytes"   "NK cells"    "T cells"    

# Number of cells in each sample and cell type
cell.per.samp.tab.perit <- t(table(perit.sce.cl$cluster_id, perit.sce.cl$sample_id))

# extract pseudobulk information for samples that pass the cell number cutoff
counts.pb.tmp.perit <- pb.perit@assays@data

# get the genes with no reads in at least half the samples out, they mess up the DESeq2 algorithm
for (i in 1:length(counts.pb.tmp.perit)) {
  my.good <- which(apply(counts.pb.tmp.perit[[i]]>0, 1, sum) >= nrow(cell.per.samp.tab.perit)/2) # see deseq2 vignette, need to remove too low genes
  counts.pb.tmp.perit[[i]] <- counts.pb.tmp.perit[[i]][my.good,]
}
########################################################

####################################################################
#######    Data preparation   ++++   perit QC Cell types     #######
####################################################################

# cell types with at least 5 cells in all samples
perit.celltype.qc <- colnames(cell.per.samp.tab.perit)[colSums(cell.per.samp.tab.perit  >= 5) == nrow(cell.per.samp.tab.perit)]
perit.celltype.qc
# [1] "B cells"     "Macrophages" "Monocytes"   "T cells"    

# extract pseudobulk information for samples that pass the cell number cutoff
counts.pb.perit <- counts.pb.tmp.perit[perit.celltype.qc]

names(counts.pb.perit) <- gsub(" ", "_", names(counts.pb.perit))

#### save counts
save(counts.pb.perit, perit.celltype.qc, file = paste0(Sys.Date(),"_muscat_PB_perit_JAXv3_QC_Clean.RData"))
##############################################################################################

##############################################################################################
# 2. Use SVA to clean up batch effects on expression and DEseq2 for DE analysis

# will run SVA to clean up noise
# run for the cell types with at least 5 cells from every each sample

######################################################
#######   DEG analysis   ++++   peritoneal     #######
######################################################

# import metadata and order it
# perit.meta     <- read.table('2025-08-28_JAXv3_sample_metadata_table.txt', header = T)
perit.meta$Sex   <- ifelse(grepl("F",perit.meta$Condition), "F", "M")
perit.meta$Age   <- ifelse(grepl("Y",perit.meta$Condition), "Y", "O")
perit.meta$Age_m <- ifelse(perit.meta$Age %in% "Y", 4, 20)
perit.meta       <- data.table::setorder(perit.meta, -Condition, Sample_ID)
rownames(perit.meta) <- perit.meta$Sample_ID

# reorder count tables in sensical order
for  (i in 1:length(counts.pb.perit)) {
  counts.pb.perit[[i]] <- counts.pb.perit[[i]][,perit.meta$Sample_ID]
}

# Create list object to receive clean SVA counts
sva.cts.perit        <- vector(mode = "list", length = length(counts.pb.perit))
names(sva.cts.perit) <- names(counts.pb.perit)

# Create list object to receive VST normalized counts
vst.cts.perit        <- vector(mode = "list", length = length(counts.pb.perit))
names(vst.cts.perit) <- names(counts.pb.perit)

# Create list object to receive DESeq2 results
deseq.res.list.perit        <- vector(mode = "list", length = length(counts.pb.perit))
names(deseq.res.list.perit) <- names(counts.pb.perit)

deseq.res.list.perit.F        <- vector(mode = "list", length = length(counts.pb.perit))
names(deseq.res.list.perit.F) <- names(counts.pb.perit)

deseq.res.list.perit.M        <- vector(mode = "list", length = length(counts.pb.perit))
names(deseq.res.list.perit.M) <- names(counts.pb.perit)

deseq.res.list.perit.FM        <- vector(mode = "list", length = length(counts.pb.perit))
names(deseq.res.list.perit.FM) <- names(counts.pb.perit)

# loop over pseudobulk data
for  (i in 1:length(counts.pb.perit)) {

  # get outprefix
  my.outprefix <- paste0(Sys.Date(),"_DEseq2_Pseudobulk_Peritoneal_",names(counts.pb.perit)[[i]])

  ###################################
  #######       Run SVA      #######

  # build design matrix
  sva.dataDesign = data.frame( row.names = perit.meta$Sample_ID ,
                               sex       = perit.meta$Sex       ,
                               age       = perit.meta$Age_m     )

  # Set null and alternative models
  mod1    = model.matrix(~ sex + age , data = sva.dataDesign)
  n.sv.be = num.sv(counts.pb.perit[[i]], mod1, method="be") # b cell is 2

  # apply SVAseq algortihm
  my.svseq = svaseq(as.matrix(counts.pb.perit[[i]]), mod1, n.sv=n.sv.be, constant = 0.1)

  # remove RIN and SV, preserve age and sex
  my.clean <- removeBatchEffect(log2(counts.pb.perit[[i]] + 0.1),
                                batch      = NULL,
                                covariates = cbind(my.svseq$sv),
                                design     = mod1[,1:3])

  # delog and round data for DEseq2 processing
  my.filtered.sva <- round(2^my.clean-0.1)

  # keep only robustly expressed genes
  sva.cts.perit[[i]] <- my.filtered.sva

  # legend
  my.cols  <- rep("",nrow(perit.meta))
  my.cols[perit.meta$Condition %in% "YF"] <- "deeppink"
  my.cols[perit.meta$Condition %in% "OF"] <- "deeppink4"
  my.cols[perit.meta$Condition %in% "YM"] <- "deepskyblue"
  my.cols[perit.meta$Condition %in% "OM"] <- "deepskyblue4"

  my.pch  <- rep(16,nrow(perit.meta))


  # get matrix using age as a modeling covariate
  dds <- DESeqDataSetFromMatrix(countData = sva.cts.perit[[i]],
                                colData   = perit.meta,
                                design    = ~ Age_m + Sex)

  # run DESeq normalizations and export results
  dds.deseq <- DESeq(dds)

  # plot dispersion
  my.disp.out <- paste(my.outprefix,"_dispersion_plot.pdf")

  pdf(my.disp.out)
  plotDispEsts(dds.deseq)
  dev.off()

  # get DESeq2 normalized expression value
  vst.cts.perit[[i]] <- getVarianceStabilizedData(dds.deseq)

  # MDS analysis
  mds.result <- cmdscale(1-cor(vst.cts.perit[[i]],method="spearman"), k = 2, eig = FALSE, add = FALSE, x.ret = FALSE)
  x <- mds.result[, 1]
  y <- mds.result[, 2]

  pdf(paste0(my.outprefix,"_MDS_plot.pdf"))
  plot(x, y,
       xlab = "MDS dimension 1", ylab = "MDS dimension 2",
       main= paste0(names(counts.pb.perit)[[i]]," (MDS)"),
       cex=3, col = my.cols, pch = my.pch,
       cex.lab = 1.25,
       cex.axis = 1.25, las = 1)
  dev.off()

  # extract gene significance by DEseq2 and 
  res.age <- results(dds.deseq, name = "Age_m") # FC

  # exclude genes with NA FDR value
  res.age <- res.age[!is.na(res.age$padj),]

  # store results
  deseq.res.list.perit[[i]]       <- data.frame(res.age)

  ### get sex dimorphic changes at FDR5
  genes.age <- rownames(res.age)[res.age$padj < 0.05]
  my.num.age <- length(genes.age)

  if (my.num.age > 2) {
    # heatmap drawing - only if there is at least 2 gene
    my.heatmap.out <- paste0(my.outprefix,"_AGING_Heatmap_FDR5_GENES.pdf")

    pdf(my.heatmap.out, onefile = F, height = 10, width = 10)
    my.heatmap.title <- paste0(names(counts.pb.perit)[[i]], " aging significant (FDR<5%), ", my.num.age, " genes")
    pheatmap::pheatmap(vst.cts.perit[[i]][genes.age,],
                       cluster_cols = F,
                       cluster_rows = T,
                       colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
                       show_rownames = F, scale="row",
                       main = my.heatmap.title,
                       cellwidth = 15,
                       border    = NA,
                       cellheight = 0.15 )
    dev.off()
  }

  # output result tables of combined analysis to text files
  my.out.ct.mat <- paste0(my.outprefix,"_AGING_VST_log2_counts_matrix.txt")
  write.table(vst.cts.perit[[i]], file = my.out.ct.mat , sep = "\t" , row.names = T, quote = F)

  my.out.stats.age <- paste0(my.outprefix,"_AGING_all_genes_statistics.txt")
  write.table(deseq.res.list.perit[[i]], file = my.out.stats.age , sep = "\t" , row.names = T, quote = F)

  my.out.fdr5.age <- paste0(my.outprefix,"_AGING_FDR5_genes_statistics.txt")
  write.table(deseq.res.list.perit[[i]][genes.age,], file = my.out.fdr5.age, sep = "\t" , row.names = T, quote = F)

  
  ################################
  ########  SPLIT by SEX  ########
  ################################
  
  # get matrix using age as a modeling covariate in FEMALES
  dds.F <- DESeqDataSetFromMatrix(countData = sva.cts.perit[[i]][,perit.meta$Sex %in% 'F'],
                                  colData = perit.meta[perit.meta$Sex %in% 'F',],
                                  design = ~ Age_m)
  
  # run DESeq normalizations and export results
  dds.deseq.F <- DESeq(dds.F)
  res.age.F   <- results(dds.deseq.F, name = "Age_m") # FC
  res.age.F   <- res.age.F[!is.na(res.age.F$padj),]
  
  # store results
  deseq.res.list.perit.F[[i]]       <- data.frame(res.age.F)
  
  # get matrix using age as a modeling covariate in MALES
  dds.M <- DESeqDataSetFromMatrix(countData = sva.cts.perit[[i]][,perit.meta$Sex %in% 'M'],
                                  colData = perit.meta[perit.meta$Sex %in% 'M',],
                                  design = ~ Age_m)
  
  # run DESeq normalizations and export results
  dds.deseq.M <- DESeq(dds.M)
  res.age.M   <- results(dds.deseq.M, name = "Age_m") # FC
  res.age.M   <- res.age.M[!is.na(res.age.M$padj),]
  
  # store results
  deseq.res.list.perit.M[[i]]       <- data.frame(res.age.M)
  
  #%%%%%%%%%%%%%% 
  colnames(res.age.F) <- paste(colnames(res.age.F),"F",sep = "_")
  colnames(res.age.M) <- paste(colnames(res.age.M),"M",sep = "_")
  my.merged.FM <- merge(data.frame(res.age.F), data.frame(res.age.M), by = "row.names")
  
  deseq.res.list.perit.FM[[i]]      <- my.merged.FM
  
  my.spear.cor <- cor.test(my.merged.FM$log2FoldChange_F,my.merged.FM$log2FoldChange_M, method = 'spearman')
  my.rho       <- signif(my.spear.cor$estimate,3)

  #### divergently regulated genes
  my.F_notM <- bitAnd(my.merged.FM$padj_F < 0.05, my.merged.FM$padj_M > 0.1) > 0
  my.M_notF <- bitAnd(my.merged.FM$padj_M < 0.05, my.merged.FM$padj_F > 0.1) > 0
  
  pdf(paste(my.outprefix,"_aging_male_vs_female_FC_scatterplot_divergent_FDR5.pdf", sep = "_"))
  smoothScatter(my.merged.FM$log2FoldChange_F,my.merged.FM$log2FoldChange_M, 
                xlim = c(-0.25,0.25), ylim = c(-0.25,0.25),
                xlab = "log2(FC) in Females with aging",
                ylab = "log2(FC) in Males with aging",
                main = names(counts.pb.perit)[[i]])
  abline(0,1, col = "grey", lty = "dashed")
  abline(h = 0, col = "red", lty = "dashed")
  abline(v = 0, col = "red", lty = "dashed")
  text(-0.25, 0.25, paste("Rho = ",my.rho), pos = 4)
  text(-0.25, 0.22, paste("p = ",signif(my.spear.cor$p.value,2)), pos = 4)
  points(my.merged.FM$log2FoldChange_F[my.F_notM],my.merged.FM$log2FoldChange_M[my.F_notM], cex = 0.5, pch = 16, col = "orchid")
  points(my.merged.FM$log2FoldChange_F[my.M_notF],my.merged.FM$log2FoldChange_M[my.M_notF], cex = 0.5, pch = 16, col = "royalblue")
  legend("bottomleft",
         c(paste("M not F,", sum(my.M_notF)),paste("F not M,", sum(my.F_notM))),
         col = c("royalblue","orchid"), pch = 16, pt.cex = 0.5,bty = 'n')
  dev.off()
  
  write.table(my.merged.FM[my.M_notF,], file = paste(my.outprefix,"AGING_Male.NOT.Female_PB_FDR5_genes_statistics.txt", sep ="_"), sep = "\t" , row.names = F, quote=F)
  write.table(my.merged.FM[my.F_notM,], file = paste(my.outprefix,"AGING_Female.NOT.Male_PB_FDR5_genes_statistics.txt", sep ="_"), sep = "\t" , row.names = F, quote=F)
  
  write.table(my.merged.FM, file = paste(my.outprefix,"AGING_Separated_Sex_Merged_Table_ALL_genes_statistics.txt", sep ="_"), sep = "\t" , row.names = F, quote=F)
  write.table(res.age.F   , file = paste(my.outprefix,"AGING_Female_ALL_genes_statistics.txt", sep ="_"), sep = "\t" , row.names = F, quote=F)
  write.table(res.age.M   , file = paste(my.outprefix,"AGING_Male_ALL_genes_statistics.txt", sep ="_"), sep = "\t" , row.names = F, quote=F)

}

# save R object with all DEseq2 results
my.rdata.age <- paste0(Sys.Date(),"_pseudobulk_JAXv3_DEseq2_objects.RData")
save(deseq.res.list.perit, deseq.res.list.perit.F, deseq.res.list.perit.M, file = my.rdata.age)

my.rdata.age.2 <- paste0(Sys.Date(),"_pseudobulk_JAXv3_MergedFM_DEseq2_objects.RData")
save(deseq.res.list.perit.FM, file = my.rdata.age.2)

my.vst.age <- paste0(Sys.Date(),"_pseudobulk_JAXv3_VST_data_objects.RData")
save(vst.cts.perit, file = my.vst.age)
################################################################################################


#######################
sink(file = paste(Sys.Date(),"_scRNAseq_peritoneal_JAXv3_session_Info.txt", sep =""))
sessionInfo()
sink()

