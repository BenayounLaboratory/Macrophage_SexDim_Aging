setwd('/Volumes/BB_HQ_1/Immune_sex_dimorphism_Aging/Macrophages/Clean_Macrophage_Data/Intervention_Data/E2_supplementation/DESeq2')
options(stringsAsFactors = FALSE)

# load libraries for analysis
library(DESeq2)
library(pheatmap)
library('bitops')
library('limma')
library(RColorBrewer)
library(fields)

# 2023-10-24
# Analyze peritoneal macrophage E2 supplementationRNAseq

####################################  peritoneal  Macrophages    ####################################
# read in count matrix
my.mph1 <- read.table('../STAR/2023-10-20_PeriMac_aging_E2_supplementation_experiment_CLEAN_counts.txt', sep = "\t", header = T, skip = 1)
my.mph.1 <- my.mph1[,-c(2:6)]
# apply(my.mph.1[,-1],2,sum)

colnames(my.mph.1) <- c("Geneid"           ,                       
                        "PeriMac_YF_CTL_1" ,
                        "PeriMac_YF_CTL_2" ,
                        "PeriMac_YF_CTL_3" ,
                        "PeriMac_YF_CTL_4" ,
                        "PeriMac_YF_CTL_5" ,
                        "PeriMac_OF_CTL_1" ,
                        "PeriMac_OF_CTL_2" ,
                        "PeriMac_OF_CTL_3" ,
                        "PeriMac_OF_CTL_4" ,
                        "PeriMac_OF_CTL_5" ,
                        "PeriMac_OF_E2_1"  ,
                        "PeriMac_OF_E2_2"  ,
                        "PeriMac_OF_E2_3"  ,
                        "PeriMac_OF_E2_4"  ,
                        "PeriMac_OF_E2_5"  )

# process RNAseq data and save RData object
my.status = c(rep("YF_CTL",5),rep("OF_CTL",5),rep("OF_E2",5))

# see deseq2 vignette, remove genes without consistent expression
my.keep <- apply(my.mph.1[,-1]> 0, 1, sum) >= 10

# Now pull out the null/low expressed genes
my.filtered.matrix           <- my.mph.1[my.keep,-1] # 15085
rownames(my.filtered.matrix) <- my.mph.1[my.keep,1]

# get output file prefixes
my.outprefix <- paste(Sys.Date(),"Peritoneal_Macrophage_Aging_E2sup_DESeq2",sep="_")

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
my.colors[my.status == "OF_CTL"] <- "deeppink4"
my.colors[my.status == "OF_E2"] <- "magenta"

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
     xlim = c(-0.015,0.04),
     ylim = c(-0.015,0.03))
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
### 1. Aging
res.aging <- results(dds.deseq, contrast = c("status","OF_CTL","YF_CTL")) # FC in old F over young F veh

### get the heatmap of sex dimorphic changes at FDR5; exclude NA
res.aging <- res.aging[!is.na(res.aging$padj),]

genes.aging <- rownames(res.aging)[res.aging$padj < 0.05]
my.num.aging <- length(genes.aging) ## 757

# heatmap drawing - only if there is at least 2 gene
my.heatmap.out <- paste(my.outprefix,"Aging_CTL_regulated_Heatmap_FDR5.pdf", sep = "_")

pdf(my.heatmap.out, onefile = F)
my.heatmap.title <- paste("PeriMac CTL Aging significant (FDR<5%), ", my.num.aging, " genes",sep="")
pheatmap(tissue.cts[genes.aging,],
         cluster_cols = F,
         cluster_rows = T,
         colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
         show_rownames = F, scale="row",
         main = my.heatmap.title, cellwidth = 15, border_color = NA)
dev.off()

# output result tables of combined analysis to text files
my.out.stats.aging <- paste0(my.outprefix,"_Aging_CTL_all_genes_statistics.txt")
write.table(res.aging, file = my.out.stats.aging , sep = "\t" , row.names = T, quote=F)

my.out.fdr5.aging <- paste0(my.outprefix,"_Aging_CTL_FDR5_genes_statistics.txt")
write.table(res.aging[genes.aging,], file = my.out.fdr5.aging, sep = "\t" , row.names = T, quote=F)

my.rdata.aging <- paste0(my.outprefix,"_Aging_CTL_DEseq2_object.RData")
save(res.aging, file = my.rdata.aging)


### 2. E2
res.e2 <- results(dds.deseq, contrast = c("status","OF_E2","OF_CTL")) # FC in E2 over vehicle old F

### get the heatmap of sex dimorphic changes at FDR5; exclude NA
res.e2 <- res.e2[!is.na(res.e2$padj),]

genes.e2 <- rownames(res.e2)[res.e2$padj < 0.05]
my.num.e2 <- length(genes.e2) # 144

# heatmap drawing - only if there is at least 2 gene
my.heatmap.out <- paste(my.outprefix,"Old_E2_regulated_Heatmap_FDR5.pdf", sep = "_")

pdf(my.heatmap.out, onefile = F)
my.heatmap.title <- paste("PeriMac OldE2 significant (FDR<5%), ", my.num.e2, " genes",sep="")
pheatmap(tissue.cts[genes.e2,],
         cluster_cols = F,
         cluster_rows = T,
         colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
         show_rownames = F, scale="row",
         main = my.heatmap.title, cellwidth = 15, border_color = NA)
dev.off()


# output result tables of combined analysis to text files
my.out.stats.e2 <- paste0(my.outprefix,"_Old_E2_all_genes_statistics.txt")
write.table(res.e2, file = my.out.stats.e2 , sep = "\t" , row.names = T, quote=F)

my.out.fdr5.e2 <- paste0(my.outprefix,"_Old_E2_FDR5_genes_statistics.txt")
write.table(res.e2[genes.e2,], file = my.out.fdr5.e2, sep = "\t" , row.names = T, quote=F)

my.rdata.e2 <- paste0(my.outprefix,"_Old_E2_DEseq2_object.RData")
save(res.e2, file = my.rdata.e2)
################################################################################################

################################################################################################
### DEseq2 VST counts
my.out.ct.mat <- paste0(my.outprefix,"_log2_counts_matrix.txt")
write.table(tissue.cts, file = my.out.ct.mat , sep = "\t" , row.names = T, quote=F)
################################################################################################


################################################################################################
# Compare aging vs. E2 effect
colnames(res.aging) <- paste(colnames(res.aging),"age",sep = "_")
colnames(res.e2   ) <- paste(colnames(res.e2)   ,"E2",sep = "_")
my.merged <- cbind(res.aging,res.e2)
my.merged <- my.merged[!is.na(my.merged$padj_age),]
my.merged <- my.merged[!is.na(my.merged$padj_E2),]

my.spear.cor <- cor.test(my.merged$log2FoldChange_age,my.merged$log2FoldChange_E2, method = 'spearman')
my.rho        <- signif(my.spear.cor$estimate,3)

#### commonly regulated genes
my.common.5  <- bitAnd(my.merged$padj_age < 0.05, my.merged$padj_E2 < 0.05) > 0


pdf(paste0(my.outprefix,"_FC_scatterplot_FDR5.pdf"))
smoothScatter(my.merged$log2FoldChange_age,my.merged$log2FoldChange_E2, 
              xlim = c(-6,6), ylim = c(-6,6),
              xlab = "log2(FC) in Females with aging (Vehicle)",
              ylab = "log2(FC) in Old Females with E2 vs. Vehicle",
              main = "Perimac aging E2 treatment")
abline(0,1, col = "grey", lty = "dashed")
abline(0,-1, col = "cyan", lty = "dashed")
abline(h = 0, col = "red", lty = "dashed")
abline(v = 0, col = "red", lty = "dashed")
text(-5.7, 5.8, paste("Rho = ",my.rho), pos = 4)
text(-5.7, 5.2, paste("p = ",signif(my.spear.cor$p.value,4)), pos = 4)
points(my.merged$log2FoldChange_age[my.common.5 ],my.merged$log2FoldChange_E2[my.common.5 ], cex = 0.5, pch = 16, col = "gold1")
legend("bottomleft",
       c(paste("FDR5%,", sum(my.common.5))),
       col = "gold1", pch = 16, pt.cex = 0.5,bty = 'n')
dev.off()



#######################
sink(file = paste(Sys.Date(),"Macrophage_E2_Aging_RNAseq_analysis_session_Info.txt", sep =""))
sessionInfo()
sink()

