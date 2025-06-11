setwd('/Volumes/BB_HQ_1/Immune_OVXorphism_Aging/Macrophages/OVX/Macrophage/DESeq2')
options(stringsAsFactors = FALSE)

# load libraries for analysis
library(DESeq2)
library(pheatmap)
library('bitops')
library('limma')
library(RColorBrewer)
library(fields)

# 2022-11-23
# Analyze peritoneal macrophage OVX cohort 1 RNAseq

####################################  peritoneal  Macrophages    ####################################
# read in count matrix
my.mph1 <- read.table('../STAR/2022-11-21_PeriMac_OVX_counts.txt', sep = "\t", header = T, skip = 1)
my.mph.1 <- my.mph1[,-c(2:6)]
# apply(my.mph.1[,-1],2,sum)

colnames(my.mph.1) <- c("Geneid",
                        "PeriMac_A1_SHAM"  ,
                        "PeriMac_A3_SHAM"  ,
                        "PeriMac_A5_SHAM"  ,
                        "PeriMac_A7_SHAM"  ,
                        "PeriMac_A9_SHAM"  ,
                        "PeriMac_A2_OVX"   ,
                        "PeriMac_A4_OVX"   ,
                        "PeriMac_A6_OVX"   ,
                        "PeriMac_A8_OVX"   ,
                        "PeriMac_A10_OVX"  ,
                        "PeriMac_A12_OVX"  )

# process RNAseq data and save RData object
my.status = c(rep("SHAM",5),rep("OVX",6))

# see deseq2 vignette, remove genes without consistent expression
my.keep <- apply(my.mph.1[,-1]> 0, 1, sum) >= 5

# Now pull out the null/low expressed genes
my.filtered.matrix           <- my.mph.1[my.keep,-1] # 14479
rownames(my.filtered.matrix) <- my.mph.1[my.keep,1]

# get output file prefixes
my.outprefix <- paste(Sys.Date(),"Peritoneal_Macrophage_OVX_DESeq2",sep="_")

# build design matrix
dataDesign = data.frame(row.names = colnames( my.filtered.matrix ), 
                        status = my.status)

# get matrix using age as a modeling covariate
dds <- DESeqDataSetFromMatrix(countData = my.filtered.matrix,
                              colData   = dataDesign,
                              design    = ~ status)

# run DESeq normalizations and export results
dds.deseq <- DESeq(dds)

# plot dispersion
my.disp.out <- paste(my.outprefix,"_dispersion_plot.pdf")

pdf(my.disp.out)
plotDispEsts(dds.deseq)
dev.off()

# get DESeq2 normalized expression value
tissue.cts <- getVarianceStabilizedData(dds.deseq)

# color-code 
my.colors <- rep("deeppink",dim(my.filtered.matrix)[2])
my.colors[my.status == "OVX"] <- "hotpink4"

# expression range
pdf(paste(my.outprefix,"_Normalized_counts_boxplot.pdf"))
boxplot(tissue.cts,col=my.colors,cex=0.5,ylab="Log2 DESeq2 Normalized counts", las = 2)  
dev.off()

# plot Xist expression
pdf(paste(my.outprefix,"_Normalized_counts_Xist_expression_barplot.pdf"))
barplot(tissue.cts["Xist",], ylab = "Normalized log2(counts) Xist expression", las = 2, col = my.colors)
dev.off()


# MDS analysis
mds.result <- cmdscale(1-cor(tissue.cts,method="spearman"), k = 2, eig = FALSE, add = FALSE, x.ret = FALSE)
x <- mds.result[, 1]
y <- mds.result[, 2]

pdf(paste0(my.outprefix,"_MDS_plot.pdf"))
plot(x, y, 
     xlab = "MDS dimension 1", ylab = "MDS dimension 2",
     main="Multi-dimensional Scaling", 
     cex=3, col= my.colors, pch= 16,
     cex.lab = 1.5,
     cex.axis = 1.5)
dev.off()

# PCA analysis
my.pos.var <- apply(tissue.cts,1,var) > 0
my.pca <- prcomp(t(tissue.cts[my.pos.var,]),scale = TRUE)
x <- my.pca$x[,1]
y <- my.pca$x[,2]

my.summary <- summary(my.pca)

pdf(paste(my.outprefix,"_PCA_plot.pdf",sep=""))
plot(x,y, 
     pch = 16, cex=3, col= my.colors, 
     xlab = paste('PC1 (', round(100*my.summary$importance[,1][2],1),"%)", sep=""),
     ylab = paste('PC2 (', round(100*my.summary$importance[,2][2],1),"%)", sep=""),
     cex.lab = 1.5,
     cex.axis = 1.5) 
dev.off()
#####

## Extract results
res.ovx <- results(dds.deseq, contrast = c("status","OVX","SHAM")) # FC in OVX over SHAM

### get the heatmap of sex dimorphic changes at FDR5; exclude NA
res.ovx <- res.ovx[!is.na(res.ovx$padj),]

genes.ovx <- rownames(res.ovx)[res.ovx$padj < 0.05]
my.num.ovx <- length(genes.ovx)

# heatmap drawing - only if there is at least 2 gene
my.heatmap.out <- paste(my.outprefix,"OVX_regulated_Heatmap_FDR5.pdf", sep = "_")

pdf(my.heatmap.out, onefile = F)
my.heatmap.title <- paste("PeriMac OVX significant (FDR<5%), ", my.num.ovx, " genes",sep="")
pheatmap(tissue.cts[genes.ovx,],
         cluster_cols = F,
         cluster_rows = T,
         colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
         show_rownames = F, scale="row",
         main = my.heatmap.title, cellwidth = 15, border_color = NA)
dev.off()


# output result tables of combined analysis to text files
my.out.ct.mat <- paste0(my.outprefix,"_log2_counts_matrix.txt")
write.table(tissue.cts, file = my.out.ct.mat , sep = "\t" , row.names = T, quote=F)

my.out.stats.ovx <- paste0(my.outprefix,"_OVX_all_genes_statistics.txt")
write.table(res.ovx, file = my.out.stats.ovx , sep = "\t" , row.names = T, quote=F)

my.out.fdr5.ovx <- paste0(my.outprefix,"_OVX_FDR5_genes_statistics.txt")
write.table(res.ovx[genes.ovx,], file = my.out.fdr5.ovx, sep = "\t" , row.names = T, quote=F)

my.rdata.ovx <- paste0(my.outprefix,"_OVX_DEseq2_object.RData")
save(res.ovx, file = my.rdata.ovx)
################################################################################################


#######################
sink(file = paste(Sys.Date(),"Macrophage_OVX_RNAseq_analysis_session_Info.txt", sep =""))
sessionInfo()
sink()

