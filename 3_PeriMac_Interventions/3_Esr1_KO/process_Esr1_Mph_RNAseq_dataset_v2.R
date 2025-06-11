setwd('/Volumes/BB_HQ_1/Immune_sex_dimorphism_Aging/Macrophages/Clean_Macrophage_Data/Intervention_Data/Esr1_KO/DESeq2/')
options(stringsAsFactors = FALSE)

# load libraries for analysis
library(DESeq2)
library(pheatmap)
library('bitops')
library('limma')
library(RColorBrewer)
library(fields)
# library(sva)

# 2024-09-23
# Analyze peritoneal macrophage Esr1 KO cohort 1 RNAseq
# Note low depth of some libraries, may fail

# 2024-10-14
# Full desired depth on all libraries

####################################  peritoneal  Macrophages    ####################################
# read in count matrix
my.mph1 <- read.table('../STAR/2024-10-1_PeriMac_Esr1KO_deeper_CLEAN_counts.txt', sep = "\t", header = T, skip = 1)
my.mph.1 <- my.mph1[,-c(2:6)]
# apply(my.mph.1[,-1],2,sum)

colnames(my.mph.1) <- gsub ("_STAR_Aligned.sortedByCoord.out.bam", "", colnames(my.mph.1))
colnames(my.mph.1) <- gsub ("_STAR_Aligned.MERGEDLANES.bam", "", colnames(my.mph.1))
colnames(my.mph.1) <- gsub ("RNAseq_", "", colnames(my.mph.1))
colnames(my.mph.1) <- gsub ("PeriMac_YF_", "", colnames(my.mph.1))

# process RNAseq data and save RData object
my.status = c(rep("KO",6),rep("WT",6))

# see deseq2 vignette, remove genes without consistent expression
my.keep <- apply(my.mph.1[,-1]> 0, 1, sum) >= 5

# Now pull out the null/low expressed genes
my.filtered.matrix           <- my.mph.1[my.keep,-1] # 15270
rownames(my.filtered.matrix) <- my.mph.1[my.keep,1]

# get output file prefixes
my.outprefix <- paste(Sys.Date(),"Peritoneal_Macrophage_Esr1_DESeq2",sep="_")

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
my.colors[my.status == "KO"] <- "hotpink4"

# expression range
pdf(paste(my.outprefix,"_Normalized_counts_boxplot.pdf"))
boxplot(tissue.cts,col=my.colors,cex=0.5,ylab="Log2 DESeq2 Normalized counts", las = 2)  
dev.off()

# plot Xist expression
pdf(paste(my.outprefix,"_Normalized_counts_Xist_expression_barplot.pdf"))
barplot(tissue.cts["Xist",], ylab = "Normalized log2(counts) Xist expression", las = 2, col = my.colors)
dev.off()

# plot Esr1 expression
pdf(paste(my.outprefix,"_Normalized_counts_Xist_expression_barplot.pdf"))
barplot(tissue.cts["Esr1",], ylab = "Normalized log2(counts) Xist expression", las = 2, col = my.colors)
dev.off()


# plot Esr1 expression
esr1.exp <- list("WT" = tissue.cts["Esr1",7:12],
                 "KO" = tissue.cts["Esr1",1:6])

pdf(paste(my.outprefix,"_boxplot_Normalized_counts_Esr1_expression.pdf"), width = 3, height = 5)
boxplot(esr1.exp, ylab = "Normalized log2(counts) Esr1 expression", las = 1, ylim = c(7,10), outline = F, col = c("deeppink","hotpink4"))
beeswarm::beeswarm(esr1.exp, pch = 16, add =T)
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
     xlim = c(-0.006,0.008),
     ylim = c(-0.005,0.005),
     las = 1)
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
res.esr <- results(dds.deseq, contrast = c("status","KO","WT")) # FC in ESR1 KO over SHAM

### get the heatmap of sex dimorphic changes at FDR5; exclude NA
res.esr <- res.esr[!is.na(res.esr$padj),]

genes.esr <- rownames(res.esr)[res.esr$padj < 0.05]
my.num.esr <- length(genes.esr) # 1586

# heatmap drawing - only if there is at least 2 gene
my.heatmap.out <- paste(my.outprefix,"ESR1_KO_regulated_Heatmap_FDR5.pdf", sep = "_")

pdf(my.heatmap.out, onefile = F)
my.heatmap.title <- paste("PeriMac ESR1 KO significant (FDR<5%), ", my.num.esr, " genes",sep="")
pheatmap(tissue.cts[genes.esr,c(7:12,1:6)],
         cluster_cols = F,
         cluster_rows = T,
         colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
         show_rownames = F, scale="row",
         main = my.heatmap.title, cellwidth = 15, border_color = NA)
dev.off()


# output result tables of combined analysis to text files
my.out.ct.mat <- paste0(my.outprefix,"_log2_counts_matrix.txt")
write.table(tissue.cts, file = my.out.ct.mat , sep = "\t" , row.names = T, quote=F)

my.out.stats.esr <- paste0(my.outprefix,"_ESR1_KO_all_genes_statistics.txt")
write.table(res.esr, file = my.out.stats.esr , sep = "\t" , row.names = T, quote=F)

my.out.fdr5.esr <- paste0(my.outprefix,"_ESR1_KO_FDR5_genes_statistics.txt")
write.table(res.esr[genes.esr,], file = my.out.fdr5.esr, sep = "\t" , row.names = T, quote=F)

my.rdata.esr <- paste0(my.outprefix,"_ESR1_KO_DEseq2_object.RData")
save(res.esr, file = my.rdata.esr)
################################################################################################


#######################
sink(file = paste(Sys.Date(),"Macrophage_ESR1_KO_RNAseq_analysis_session_Info.txt", sep =""))
sessionInfo()
sink()

