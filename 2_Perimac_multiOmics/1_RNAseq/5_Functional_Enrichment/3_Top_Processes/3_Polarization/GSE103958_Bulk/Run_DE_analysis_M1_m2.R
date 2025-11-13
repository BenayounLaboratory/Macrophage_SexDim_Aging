setwd('/Volumes/BB_Home_HQ/Immune_sex_dimorphism_Aging/Macrophages/Macrophage_datasets_for_Comparison/M1_M2_macrophages/GSE103958_Bulk/DEseq2')
options(stringsAsFactors = FALSE)

# load libraries for analysis
library(DESeq2)
library(pheatmap)
library('pvclust')
library('bitops')
library(RColorBrewer)
library(fields)

# 2021-11-18
# Analyze to get genes lists of M1- and M2 up genes

my.bmdm.1 <- read.table("../STAR/2021-11-12_GSE103958_BMDM_Polarization_counts.txt", header = T, sep = "\t", skip = 1)
my.bmdm   <- my.bmdm.1[,-c(2:6)]

colnames(my.bmdm) <- c("GeneName" ,
                       "GSM2786903_BMDM_Control_rep1"     ,
                       "GSM2786904_BMDM_Control_rep2"     ,
                       "GSM2786905_BMDM_Control_rep3"     ,
                       "GSM2786906_BMDM_IFNgLPS_1h_rep_1" ,
                       "GSM2786907_BMDM_IFNgLPS_1h_rep_2" ,
                       "GSM2786908_BMDM_IFNgLPS_1h_rep_3" ,
                       "GSM2786909_BMDM_IFNgLPS_4h_rep_1" ,
                       "GSM2786910_BMDM_IFNgLPS_4h_rep_2" ,
                       "GSM2786911_BMDM_IFNgLPS_4h_rep_3" ,
                       "GSM2786912_BMDM_IFNgLPS_12h_rep_1",
                       "GSM2786913_BMDM_IFNgLPS_12h_rep_2",
                       "GSM2786914_BMDM_IFNgLPS_12h_rep_3",
                       "GSM2786915_BMDM_IFNgLPS_24h_rep_1",
                       "GSM2786916_BMDM_IFNgLPS_24h_rep_2",
                       "GSM2786917_BMDM_IFNgLPS_24h_rep_3",
                       "GSM2786918_BMDM_IL4_12h_rep_1"    ,
                       "GSM2786919_BMDM_IL4_12h_rep_2"    ,
                       "GSM2786920_BMDM_IL4_12h_rep_3"    ,
                       "GSM2786921_BMDM_IL13_12h_rep_1"   ,
                       "GSM2786922_BMDM_IL13_12h_rep_2"   ,
                       "GSM2786923_BMDM_IL13_12h_rep_3"   )

rownames(my.bmdm) <- my.bmdm$GeneName

# get the genes with no reads out
my.good            <- which(apply(my.bmdm[,-1]>0, 1, sum) >= 10) # see deseq2 vignette, need to remove too low genes
my.bmdm.filt       <- my.bmdm[my.good,-1] # 13031  genes

# separate based on M1/M2
my.bmdm.lps <- my.bmdm.filt[,c(1:15)]
my.bmdm.IL  <- my.bmdm.filt[,c(1:3,16:21)]
########################################################################################################################


########################################################################################################################
########################################################################################################################
##########  1. M1 genes analysis

# Looking for M1 genes by looking for linear regression by time of LPS
# i.e. increase in expression with time in M2 

my.lps.time <- c(rep(0,3),rep(1,3),rep(4,3),rep(12,3),rep(24,3))

# design matrix
dataDesign.m1 <- data.frame( row.names = colnames( my.bmdm.lps), 
                             lps = my.lps.time)

# get matrix using age as a modeling covariate
dds.m1 <- DESeqDataSetFromMatrix(countData = my.bmdm.lps,
                                 colData   = dataDesign.m1,
                                 design    = ~ lps)

# run DESeq normalizations and export results
dds.deseq.m1 <- DESeq(dds.m1)

# normalized expression value
tissue.cts.m1 <- getVarianceStabilizedData(dds.deseq.m1)

# color-code 
my.colors <- c(rep("burlywood1",3),rep("chocolate1",3),rep("chocolate2",3),rep("chocolate3",3),rep("chocolate4",3))

# do MDS analysis
mds.result <- cmdscale(1-cor(tissue.cts.m1,method="spearman"), k = 2, eig = FALSE, add = FALSE, x.ret = FALSE)
x <- mds.result[, 1]
y <- mds.result[, 2]

my.mds.out <- paste0(Sys.Date(),"_BMDM_LPS_IFNg_MDS_plot.pdf")
pdf(my.mds.out)
plot(x, y,
     xlab = "MDS dimension 1", ylab = "MDS dimension 2",
     main="MDS - BMDM (M1)",
     cex=2, pch = 16, col = my.colors,
     cex.lab = 1.5,
     cex.axis = 1.5)
dev.off()

# barplot(tissue.cts.m1["Xist",], ylab = "Normalized log2(counts) Xist expression", las = 2, col = my.colors)
# barplot(tissue.cts.m1["Ddx3y",], ylab = "Normalized log2(counts) Ddx3y expression", las = 2, col = my.colors)

# get results
res.m1 <- results(dds.deseq.m1, name = "lps")
res.m1 <- res.m1[!is.na(res.m1$padj),]

genes.m1 <- rownames(res.m1)[res.m1$padj < 0.05]
my.num.m1 <- length(genes.m1)  # 3251

my.heatmap.out <- paste0(Sys.Date(),"_BMDM_LPS_IFNg_Heatmap_FDR5.pdf")
pdf(my.heatmap.out, onefile = F)
my.heatmap.title <- paste("LPS/IFNg significant (FDR<5%), ", my.num.m1, " genes",sep="")
pheatmap(tissue.cts.m1[genes.m1,],
         cluster_cols = F,
         cluster_rows = T,
         colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
         show_rownames = F, scale="row",
         main = my.heatmap.title, cellwidth = 15)
dev.off()

save(res.m1, file = paste(Sys.Date(),"BMDM_LPS_IFNg_DEseq2.RData", sep ="_"))




########################################################################################################################
########################################################################################################################
##########  2. M2 genes analysis

# Looking for M2 genes by looking at CTL vs IL4/IL13

my.M2   <- c(rep("M0",3),rep("M2",6))
my.cyto <- c(rep("NONE",3),rep("IL4",3),rep("IL13",3))

# design matrix
dataDesign.m2 <- data.frame( row.names = colnames( my.bmdm.IL), 
                             m2  = my.M2,
                             cyt = my.cyto)

# get matrix using age as a modeling covariate
dds.m2 <- DESeqDataSetFromMatrix(countData = my.bmdm.IL,
                                 colData   = dataDesign.m2,
                                 design    = ~ m2)

# run DESeq normalizations and export results
dds.deseq.m2 <- DESeq(dds.m2)

# normalized expression value
tissue.cts.m2 <- getVarianceStabilizedData(dds.deseq.m2)

# color-code 
my.colors <- c(rep("burlywood1",3),rep("cyan4",3),rep("cornflowerblue",3))

# do MDS analysis
mds.result <- cmdscale(1-cor(tissue.cts.m2,method="spearman"), k = 2, eig = FALSE, add = FALSE, x.ret = FALSE)
x <- mds.result[, 1]
y <- mds.result[, 2]

my.mds.out <- paste0(Sys.Date(),"_BMDM_IL4_IL13_MDS_plot.pdf")
pdf(my.mds.out)
plot(x, y,
     xlab = "MDS dimension 1", ylab = "MDS dimension 2",
     main="MDS - BMDM (m2)",
     cex=2, pch = 16, col = my.colors,
     cex.lab = 1.5,
     cex.axis = 1.5)
dev.off()

# get results
res.m2 <- results(dds.deseq.m2, contrast = c("m2","M2","M0")) # FC in M2 over M0
res.m2 <- res.m2[!is.na(res.m2$padj),]

genes.m2 <- rownames(res.m2)[res.m2$padj < 0.05]
my.num.m2 <- length(genes.m2)  # 1343

my.heatmap.out <- paste0(Sys.Date(),"_BMDM_IL4_IL13_Heatmap_FDR5.pdf")
pdf(my.heatmap.out, onefile = F)
my.heatmap.title <- paste("IL4/IL13 significant (FDR<5%), ", my.num.m2, " genes",sep="")
pheatmap(tissue.cts.m2[genes.m2,],
         cluster_cols = F,
         cluster_rows = T,
         colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
         show_rownames = F, scale="row",
         main = my.heatmap.title, cellwidth = 15)
dev.off()

save(res.m2, file = paste(Sys.Date(),"BMDM_IL4_IL13_DEseq2.RData", sep ="_"))
################################################################################################################################



################################################################################################################################
# get up genes for M1 and M2
my.m1.up <- rownames(res.m1)[bitAnd(res.m1$padj < 0.05, res.m1$log2FoldChange >0) >0] # 1392
my.m2.up <- rownames(res.m2)[bitAnd(res.m2$padj < 0.05, res.m2$log2FoldChange >0) >0] #  760


polar.gset <- list("BMDM_M1_enriched"  = my.m1.up,
                   "BMDM_M2_enriched"  = my.m2.up)


for (i in 1:length(polar.gset)) {
  
  gmtline <- paste(as.vector(c(names(polar.gset)[i],"GSE103958",polar.gset[[i]])),collapse = "\t")
  
  write.table(gmtline, file = "2021-11-18_GSE103958_BMDM_polarization_genes_lists.gmt", quote = F, row.names = F, col.names = F, append = T)
  
}
################################################################################################

#######################
sink(file = paste0(Sys.Date(),"_session_Info.txt"))
sessionInfo()
sink()


