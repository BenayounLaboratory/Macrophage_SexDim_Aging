# setwd('/Volumes/BB_HQ_1/Immune_sex_dimorphism_Aging/Macrophages/E2_treaments/RNAseq/DEseq2')
setwd('/Volumes/BB_Travel_2/J774_E2/DESeq2')
options(stringsAsFactors = FALSE)

# load libraries for analysis
library(DESeq2)
library(pheatmap)
library('bitops')
library('limma')
library(RColorBrewer)
library(fields)

# 2025-09-01
# Analyze J7 E2 supplementation mRNAseq

# 2025-11-24
# rerun using MM10 genome to match the rest of the study

####################################  peritoneal  Macrophages    ####################################
# read in count matrix
my.mph1 <- read.table('../STAR/2025-11-24_J774A_E2_CLEAN_counts.txt', sep = "\t", header = T, skip = 1)
my.mph.1 <- my.mph1[,-c(2:6)]
# apply(my.mph.1[,-1],2,sum)

colnames(my.mph.1) <- c("Geneid"           ,                       
                        "J774A.1_CTL_1" ,
                        "J774A.1_CTL_2" ,
                        "J774A.1_CTL_3" ,
                        "J774A.1_E2_1" ,
                        "J774A.1_E2_2" ,
                        "J774A.1_E2_3" ,
                        "J774A.1_E2_4" )

# process RNAseq data and save RData object
my.status = c(rep("CTL",3),rep("E2",4))

# see deseq2 vignette, remove genes without consistent expression
my.keep <- apply(my.mph.1[,-1]> 0, 1, sum) >= 5

# Now pull out the null/low expressed genes
my.filtered.matrix           <- my.mph.1[my.keep,-1] # 13125
rownames(my.filtered.matrix) <- my.mph.1[my.keep,1]

# get output file prefixes
my.outprefix <- paste(Sys.Date(),"J774A.1_E2_treatment_DESeq2_MM10",sep="_")

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
my.colors[my.status == "E2"] <- "mediumorchid1"

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
     cex.axis = 1.5,
     xlim = c(-0.02,0.03),
     ylim = c(-0.005,0.005)
     )
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

###########################################################################
## Extract results

### E2
res.e2 <- results(dds.deseq, contrast = c("status","E2","CTL")) # FC in E2 over vehicle old F

### get the heatmap of sex dimorphic changes at FDR5; exclude NA
res.e2 <- res.e2[!is.na(res.e2$padj),]

genes.e2 <- rownames(res.e2)[res.e2$padj < 0.05]
my.num.e2 <- length(genes.e2) # 32

# heatmap drawing - only if there is at least 2 gene
my.heatmap.out <- paste(my.outprefix,"E2_regulated_Heatmap_FDR5.pdf", sep = "_")

pdf(my.heatmap.out, onefile = F)
my.heatmap.title <- paste("E2 significant (FDR<5%), ", my.num.e2, " genes",sep="")
pheatmap(tissue.cts[genes.e2,],
         cluster_cols = F,
         cluster_rows = T,
         colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
         show_rownames = F, scale="row",
         main = my.heatmap.title, cellwidth = 15, border_color = NA)
dev.off()


# output result tables of combined analysis to text files
my.out.stats.e2 <- paste0(my.outprefix,"_E2_all_genes_statistics.txt")
write.table(res.e2, file = my.out.stats.e2 , sep = "\t" , row.names = T, quote=F)

my.out.fdr5.e2 <- paste0(my.outprefix,"_E2_FDR5_genes_statistics.txt")
write.table(res.e2[genes.e2,], file = my.out.fdr5.e2, sep = "\t" , row.names = T, quote=F)

my.rdata.e2 <- paste0(my.outprefix,"_E2_DEseq2_object.RData")
save(res.e2, file = my.rdata.e2)
################################################################################################

################################################################################################
### DEseq2 VST counts
my.out.ct.mat <- paste0(my.outprefix,"_log2_counts_matrix.txt")
write.table(tissue.cts, file = my.out.ct.mat , sep = "\t" , row.names = T, quote=F)
################################################################################################


#######################
sink(file = paste(Sys.Date(),"J774A.1_E2_RNAseq_analysis_session_Info.txt", sep =""))
sessionInfo()
sink()

