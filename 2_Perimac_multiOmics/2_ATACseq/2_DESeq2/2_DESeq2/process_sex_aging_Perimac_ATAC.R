setwd('/Volumes/BB_Home_HQ/Immune_sex_dimorphism_Aging/Macrophages/Clean_Macrophage_Data/ATACseq/DEseq2')
options(stringsAsFactors = FALSE)

# load libraries for analysis
library(DESeq2)
library(pheatmap)
library('pvclust')
library('bitops')
library('sva')
library('limma')
library(RColorBrewer)
library(fields)


# 2021-11-05
# process ATAC-seq dataset
# use Diffbind count on meta peaks
# run only high purity (> 95%)
# animal #1 for AC43 animals was < 95% purity QC


####################################   Neutrophils   ####################################
# read in count matrix
my.mph1 <- read.csv('../Diffbind/2021-10-14_PeriMac_ATAC_Sex_Aging_Normalized_count_matrix.txt', sep = "\t", header = T)
my.mph1$PeakName <- paste(my.mph1$seqnames,my.mph1$start,my.mph1$end,sep = ":")
my.mph <- my.mph1[,c(23,6:22)]

colnames(my.mph) <- c("PeakName"    ,
                      "YF_PeriMac08",
                      "YF_PeriMac09",
                      "YF_PeriMac16",
                      "YF_PeriMac17",
                      "OF_PeriMac03",
                      "OF_PeriMac06",
                      "OF_PeriMac11",
                      "OF_PeriMac14",
                      "YM_PeriMac02",
                      "YM_PeriMac07",
                      "YM_PeriMac10",
                      "YM_PeriMac15",
                      "YM_PeriMac18",
                      "OM_PeriMac04",
                      "OM_PeriMac05",
                      "OM_PeriMac12",
                      "OM_PeriMac13")

# get HOMER annotations
my.peak.annot <- read.csv('../Diffbind/HOMER_2021-10-14_PeriMac_ATAC_Sex_Aging_diffbind_peaks.xls', sep = "\t", header = T)
colnames(my.peak.annot)[1] <- "PeakID"
my.peak.annot$PeakName <- paste(my.peak.annot$Chr,my.peak.annot$Start-1,my.peak.annot$End,sep = ":")

# round counts (DESeq needs integers)
my.mph[,2:18] <- round(my.mph[,2:18])
rownames(my.mph) <- my.mph$PeakName

# get the peaks with no reads out
my.good <- rownames(my.mph)[apply(my.mph[,-1] > 0, 1, sum) > 10 ] # see deseq2 vignette, need to remove too low genes
my.filtered.matrix <- my.mph[my.good,-1] # 73981 peaks

# color-code
my.colors <- rep("deeppink",17)
my.colors[grep("OF",colnames(my.filtered.matrix))] <- "deeppink4"
my.colors[grep("YM",colnames(my.filtered.matrix))] <- "deepskyblue"
my.colors[grep("OM",colnames(my.filtered.matrix))] <- "deepskyblue4"

my.Age <- c(rep(4,4),rep(20,4),rep(4,5),rep(20,4))
my.Sex <- c(rep("F",8),rep("M",9))

################################################################################################
################################################################################################

############################################################
############# A. model age and sex together ################
############################################################

#########################
# DESeq2 on cleaned data
my.outprefix <- paste(Sys.Date(),"PeriMac_ATACseq_DESeq2_Analysis",sep="_")

# design matrix
dataDesign = data.frame( row.names = colnames( my.filtered.matrix ),
                         age = my.Age,
                         sex = my.Sex)

# get matrix using age as a modeling covariate
dds <- DESeqDataSetFromMatrix(countData = my.filtered.matrix,
                              colData = dataDesign,
                              design = ~ age + sex)

# run DESeq normalizations and export results
dds.deseq <- DESeq(dds)

# plot dispersion
my.disp.out <- paste(my.outprefix,"dispersion_plot.pdf",sep="_")

pdf(my.disp.out)
plotDispEsts(dds.deseq)
dev.off()

# normalized expression value
tissue.cts <- getVarianceStabilizedData(dds.deseq)

# do MDS analysis
mds.result <- cmdscale(1-cor(tissue.cts,method="spearman"), k = 2, eig = FALSE, add = FALSE, x.ret = FALSE)
x <- mds.result[, 1]
y <- mds.result[, 2]

my.mds.out <- paste(my.outprefix,"MDS_plot.pdf",sep="_")
pdf(my.mds.out)
plot(x, y,
     xlab = "MDS dimension 1", ylab = "MDS dimension 2",
     main="ATACseq height (MDS)",
     cex=3, pch = 16, col = my.colors,
     xlim = c(-0.03,0.06),
     ylim = c(-0.02,0.02),
     cex.lab = 1.5,
     cex.axis = 1.5,
     las = 1)
dev.off()

# PCA analysis
my.pos.var <- apply(tissue.cts,1,var) > 0
my.pca <- prcomp(t(tissue.cts[my.pos.var,]),scale = TRUE)
x <- my.pca$x[,1]
y <- my.pca$x[,2]

my.summary <- summary(my.pca)

my.pca.out <- paste(my.outprefix,"PCA_plot.pdf",sep="_")
pdf(my.pca.out)
plot(x,y,
     cex=3, pch = 16, col = my.colors,
     xlab = paste('PC1 (', round(100*my.summary$importance[,1][2],1),"%)", sep=""),
     ylab = paste('PC2 (', round(100*my.summary$importance[,2][2],1),"%)", sep=""),
     #xlim = c(-120,120),
     #ylim = c(-100,100),
     cex.lab = 1.5,
     cex.axis = 1.5,
     main="ATACseq height (PCA)")
dev.off()


# expression range
pdf(paste(my.outprefix,"_Normalized_ATACseq_boxplot.pdf"))
boxplot(tissue.cts,col=my.colors,cex=0.5,ylab="Log2 DESeq2 Normalized counts", las = 2, outline = F)
dev.off()

# plot X and Y ATACseq peak depths
pdf(paste(my.outprefix,"_Normalized_counts_Xist_ATACseq_barplot.pdf"))
barplot(apply(tissue.cts[rownames(tissue.cts) %in% my.peak.annot$PeakName[my.peak.annot$Gene.Name %in% "Xist"],],2,sum), ylab = "Normalized log2(counts)  Xist loci", las = 2, col = my.colors)
dev.off()

pdf(paste(my.outprefix,"_Normalized_counts_Y_chromosome_ATACseq_barplot.pdf"))
barplot(apply(tissue.cts[rownames(tissue.cts) %in% my.peak.annot$PeakName[my.peak.annot$Chr %in% "chrY"],],2,sum), ylab = "Normalized log2(counts) Y chromosome", las = 2, col = my.colors)
dev.off()

###############################################################################################
## a. model aging with sex as covariate  %%%%%%%%%%%%%%
res.age <- results(dds.deseq, name= "age")

### get the heatmap of aging changes at FDR5; exclude NA
res.age <- res.age[!is.na(res.age$padj),]

genes.aging <- rownames(res.age)[res.age$padj < 0.05]
my.num.aging <- length(genes.aging) # 4522

my.heatmap.out <- paste(my.outprefix,"AGING_Heatmap_FDR5.pdf", sep = "_")
pdf(my.heatmap.out, onefile = F)
my.heatmap.title <- paste("Aging significant (FDR<5%), ", my.num.aging, " ATACseq peaks (height)",sep="")
pheatmap(tissue.cts[genes.aging,],
         cluster_cols = F,
         cluster_rows = T,
         colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
         show_rownames = F, scale="row",
         main = my.heatmap.title, cellwidth = 15)
dev.off()

save(res.age, file = paste(Sys.Date(),"PeriMac_ATACseq_Aging_BOTH.RData", sep ="_"))

###############################################################################################
## b. sex with age as covariate
res.sex <- results(dds.deseq, contrast = c("sex","F","M")) # FC in females over Males

### get the heatmap of sex dimorphic changes at FDR5; exclude NA
res.sex <- res.sex[!is.na(res.sex$padj),]

genes.sex <- rownames(res.sex)[res.sex$padj < 0.05]
my.num.sex <- length(genes.sex) # 4651

my.heatmap.out <- paste(my.outprefix,"SEX_DIM_Heatmap_FDR5.pdf", sep = "_")
pdf(my.heatmap.out, onefile = F)
my.heatmap.title <- paste("Sex significant (FDR<5%), ", my.num.sex, " ATACseq peaks (height)",sep="")
pheatmap(tissue.cts[genes.sex,],
         cluster_cols = F,
         cluster_rows = T,
         colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
         show_rownames = F, scale="row",
         main = my.heatmap.title, cellwidth = 15)
dev.off()

save(res.sex, file = paste(Sys.Date(),"PeriMac_ATACseq_SEX.RData", sep ="_"))


# output result tables of combined analysis to text files
my.out.ct.mat <- paste(my.outprefix,"_log2_counts_matrix_DEseq2.txt", sep = "_")
write.table(tissue.cts, file = my.out.ct.mat , sep = "\t" , row.names = T, quote=F)

my.out.stats.age <- paste(my.outprefix,"AGING_all_genes_statistics.txt", sep = "_")
my.out.stats.sex <- paste(my.outprefix,"SEX_DIM_all_genes_statistics.txt", sep = "_")
write.table(res.age, file = my.out.stats.age , sep = "\t" , row.names = T, quote=F)
write.table(res.sex, file = my.out.stats.sex , sep = "\t" , row.names = T, quote=F)

my.out.fdr5.age <- paste(my.outprefix,"AGING_FDR5_genes_statistics.txt", sep = "_")
my.out.fdr5.sex <- paste(my.outprefix,"SEX_DIM_FDR5_genes_statistics.txt", sep = "_")
write.table(res.age[genes.aging,], file = my.out.fdr5.age, sep = "\t" , row.names = T, quote=F)
write.table(res.sex[genes.sex,], file = my.out.fdr5.sex, sep = "\t" , row.names = T, quote=F)
################################################################################################

################################################################################################
################################################################################################
# get output file prefixes
my.outprefix.age_M <- paste(Sys.Date(),"PeriMac_ATACseq_DESeq2_model_with_AGING_Males",sep="_")
my.outprefix.age_F <- paste(Sys.Date(),"PeriMac_ATACseq_DESeq2_model_with_AGING_Females",sep="_")

###################################################################
############# B. model age in each sex separately  ################
###################################################################

# get matrix using age as a modeling covariate in FEMALES
dds.F <- DESeqDataSetFromMatrix(countData = my.filtered.matrix[,my.Sex %in% 'F'],
                                colData = dataDesign[my.Sex %in% 'F',],
                                design = ~ age)

# run DESeq normalizations and export results
dds.deseq.F <- DESeq(dds.F)
res.age.F <- results(dds.deseq.F, name= "age")

# get matrix using age as a modeling covariate in MALES
dds.M <- DESeqDataSetFromMatrix(countData = my.filtered.matrix[,my.Sex %in% 'M'],
                                colData = dataDesign[my.Sex %in% 'M',],
                                design = ~ age)

# run DESeq normalizations and export results
dds.deseq.M <- DESeq(dds.M)
res.age.M <- results(dds.deseq.M, name= "age")

my.rdata.age.F <- paste0(my.outprefix.age_F,"_DEseq2_object.RData")
my.rdata.age.M <- paste0(my.outprefix.age_M,"_DEseq2_object.RData")
save(res.age.F, file = my.rdata.age.F)
save(res.age.M, file = my.rdata.age.M)

#%%%%%%%%%%%%%% 
colnames(res.age.F) <- paste(colnames(res.age.F),"F",sep = "_")
colnames(res.age.M) <- paste(colnames(res.age.M),"M",sep = "_")
my.merged.FM <- cbind(res.age.F,res.age.M)
my.merged.FM <- my.merged.FM[!is.na(my.merged.FM$padj_F),]
my.merged.FM <- my.merged.FM[!is.na(my.merged.FM$padj_M),]

res.age.F <- res.age.F[!is.na(res.age.F$padj_F),]
res.age.M <- res.age.M[!is.na(res.age.M$padj_M),]

my.spear.cor <- cor.test(my.merged.FM$log2FoldChange_F,my.merged.FM$log2FoldChange_M, method = 'spearman')
my.rho <- signif(my.spear.cor$estimate,3)
my.p <- signif(my.spear.cor$p.value,3)

#### commonly regulated genes
my.F_M.5 <- bitAnd(my.merged.FM$padj_F < 0.05, my.merged.FM$padj_M < 0.05) > 0
my.F_M.10 <- bitAnd(my.merged.FM$padj_F < 0.1, my.merged.FM$padj_M < 0.1) > 0

pdf(paste(Sys.Date(),"PeriMac_ATACseq_aging_male_vs_female_FC_scatterplot_FDR5_10.pdf", sep = "_"))
smoothScatter(my.merged.FM$log2FoldChange_F,my.merged.FM$log2FoldChange_M, 
              xlim = c(-0.4,0.4), ylim = c(-0.4,0.4),
              xlab = "log2(FC per m) in Females with aging",
              ylab = "log2(FC per m) in Males with aging",
              main = "PeriMac ATACseq height")
abline(0,1, col = "grey", lty = "dashed")
abline(h = 0, col = "red", lty = "dashed")
abline(v = 0, col = "red", lty = "dashed")
text(-0.3, 0.28, paste("Rho = ",my.rho))
text(-0.3, 0.23, paste("p = ",my.p))
points(my.merged.FM$log2FoldChange_F[my.F_M.10],my.merged.FM$log2FoldChange_M[my.F_M.10], cex = 0.5, pch = 16, col = "gold1")
points(my.merged.FM$log2FoldChange_F[my.F_M.5],my.merged.FM$log2FoldChange_M[my.F_M.5], cex = 0.5, pch = 16, col = "orange3")
legend("topleft",
       c(paste("FDR5%,", sum(my.F_M.5)),paste("FDR10%,", sum(my.F_M.10))),
       col = c("orange3","gold1"), pch = 16, pt.cex = 0.5,bty = 'n')
dev.off()

write.table(my.merged.FM[my.F_M.10,], file = paste(Sys.Date(),"PeriMac_ATACseq_AGING_SexAgreement_FDR10_genes_statistics.txt", sep ="_"), sep = "\t" , row.names = T, quote=F)


#### divergently regulated genes
my.F_notM <- bitAnd(my.merged.FM$padj_F < 0.05, my.merged.FM$padj_M > 0.1) > 0
my.M_notF <- bitAnd(my.merged.FM$padj_M < 0.05, my.merged.FM$padj_F > 0.1) > 0

pdf(paste(Sys.Date(),"PeriMac_ATACseq_aging_male_vs_female_FC_scatterplot_divergent_FDR5.pdf", sep = "_"))
smoothScatter(my.merged.FM$log2FoldChange_F,my.merged.FM$log2FoldChange_M, 
              xlim = c(-0.3,0.3), ylim = c(-0.3,0.3),
              xlab = "log2(FC per m) in Females with aging",
              ylab = "log2(FC per m) in Males with aging",
              main = "PeriMac ATACseq height")
abline(0,1, col = "grey", lty = "dashed")
abline(h = 0, col = "red", lty = "dashed")
abline(v = 0, col = "red", lty = "dashed")
text(-0.4, 0.4, paste("Rho = ",my.rho))
points(my.merged.FM$log2FoldChange_F[my.F_notM],my.merged.FM$log2FoldChange_M[my.F_notM], cex = 0.5, pch = 16, col = "orchid")
points(my.merged.FM$log2FoldChange_F[my.M_notF],my.merged.FM$log2FoldChange_M[my.M_notF], cex = 0.5, pch = 16, col = "royalblue")
legend("topleft",
       c(paste("M not F,", sum(my.M_notF)),paste("F not M,", sum(my.F_notM))),
       col = c("royalblue","orchid"), pch = 16, pt.cex = 0.5,bty = 'n')
dev.off()

write.table(my.merged.FM[my.M_notF,], file = paste(Sys.Date(),"PeriMac_ATACseq_AGING_Male.NOT.Female_FDR5_genes_statistics.txt", sep ="_"), sep = "\t" , row.names = T, quote=F)
write.table(my.merged.FM[my.F_notM,], file = paste(Sys.Date(),"PeriMac_ATACseq_AGING_Female.NOT.Male_FDR5_genes_statistics.txt", sep ="_"), sep = "\t" , row.names = T, quote=F)

write.table(my.merged.FM, file = paste(Sys.Date(), "PeriMac_ATACseq_AGING_Separated_Sex_Merged_Table_ALL_genes_statistics.txt", sep ="_"), sep = "\t" , row.names = T, quote=F)
write.table(res.age.F, file = paste(Sys.Date()   , "PeriMac_ATACseq_AGING_Female_ALL_genes_statistics.txt", sep ="_"), sep = "\t" , row.names = T, quote=F)
write.table(res.age.M, file = paste(Sys.Date()   , "PeriMac_ATACseq_AGING_Male_ALL_genes_statistics.txt", sep ="_"), sep = "\t" , row.names = T, quote=F)


my.M_notF.peak <- rownames(tissue.cts)[rownames(tissue.cts) %in% rownames(my.merged.FM)[my.M_notF]]
my.F_notM.peak <- rownames(tissue.cts)[rownames(tissue.cts) %in% rownames(my.merged.FM)[my.F_notM]]

my.heatmap.out <- paste(my.outprefix,"Female_Bias_Heatmap_FDR5.pdf", sep = "_")
pdf(my.heatmap.out, onefile = F)
my.heatmap.title <- paste("Female bias ", length(my.F_notM.peak), " ATACseq peaks (height)",sep="")
pheatmap(tissue.cts[my.F_notM.peak,],
         cluster_cols = F,
         cluster_rows = T,
         colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
         show_rownames = F, scale="row",
         main = my.heatmap.title, cellwidth = 15)
dev.off()

################################################################################################
################################################################################################



################################################################################################
################################################################################################
# annotate Peaks, output bed files
library(bitops)

##### Log2 normalized counts
tissue.cts <- cbind(rownames(tissue.cts),tissue.cts)
colnames(tissue.cts)[1] <- "PeakName"
tissue.cts.annot <- data.frame(merge(my.peak.annot[,c("PeakName","Chr", "Start", "End", "Annotation", "Distance.to.TSS", "Gene.Name","Gene.Description" )],tissue.cts))

my.out.ct.mat <- paste(my.outprefix,"log2_counts_matrix_DEseq2_PeakAnnot.txt", sep = "_")
write.table(tissue.cts, file = my.out.ct.mat , sep = "\t" , row.names = T, quote=F)

##### Output SEX bed
res.sex <- cbind(rownames(res.sex),res.sex)
colnames(res.sex)[1] <- "PeakName"
res.sex.annot <- data.frame(merge(my.peak.annot[,c("PeakName","Chr", "Start", "End", "Annotation", "Distance.to.TSS", "Gene.Name","Gene.Description" )],res.sex))

my.out.stats.sex <- paste(my.outprefix,"SEX_DIM_all_genes_statistics_PeakAnnot.txt",sep = "_")
write.table(res.sex.annot, file = my.out.stats.sex , sep = "\t" , row.names = F, quote=F)

my.out.fdr5.sex <- paste(my.outprefix,"SEX_DIM_FDR5_genes_statistics_PeakAnnot.txt",sep = "_")
write.table(res.sex.annot[res.sex.annot$padj < 0.05,], file = my.out.fdr5.sex, sep = "\t" , row.names = F, quote=F)

my.out.bckgd.sex.BED <- paste0(my.outprefix,"_SEX_DIM_Background.bed")
write.table(res.sex.annot[,c("Chr", "Start", "End","PeakName")], file = my.out.bckgd.sex.BED , sep = "\t" , row.names = F, col.names = F, quote=F)

my.out.fdr5.sex.BED <- paste0(my.outprefix,"_SEX_DIM_FDR5.bed")
write.table(res.sex.annot[res.sex.annot$padj < 0.05,c("Chr", "Start", "End","PeakName")], file = my.out.fdr5.sex.BED , sep = "\t" , row.names = F, col.names = F, quote=F)

my.out.fdr5.sex.BED.up <- paste0(my.outprefix,"_SEX_DIM_FDR5_UP.bed")
write.table(res.sex.annot[bitAnd(res.sex.annot$padj < 0.05,res.sex.annot$log2FoldChange >0)>0,c("Chr", "Start", "End","PeakName")], file = my.out.fdr5.sex.BED.up , sep = "\t" , row.names = F, col.names = F, quote=F)

my.out.fdr5.sex.BED.dwn <- paste0(my.outprefix,"_SEX_DIM_FDR5_DWN.bed")
write.table(res.sex.annot[bitAnd(res.sex.annot$padj < 0.05,res.sex.annot$log2FoldChange <0)>0,c("Chr", "Start", "End","PeakName")], file = my.out.fdr5.sex.BED.dwn , sep = "\t" , row.names = F, col.names = F, quote=F)

##### Output AGE bed
res.age <- cbind(rownames(res.age),res.age)
colnames(res.age)[1] <- "PeakName"

res.age.annot <- merge(my.peak.annot[,c("PeakName","Chr", "Start", "End", "Annotation", "Distance.to.TSS", "Gene.Name","Gene.Description" )],res.age)

my.out.stats.age <- paste(my.outprefix,"AGING_all_genes_statistics_PeakAnnot.txt",sep = "_")
write.table(res.age.annot, file = my.out.stats.age , sep = "\t" , row.names = F, quote=F)

my.out.fdr5.age <- paste(my.outprefix,"AGING_FDR5_genes_statistics_PeakAnnot.txt",sep = "_")
write.table(res.age.annot[res.age.annot$padj < 0.05,], file = my.out.fdr5.age, sep = "\t" , row.names = F, quote=F)

my.out.bckgd.age.BED <- paste0(my.outprefix,"_AGING_Background.bed")
write.table(res.age.annot[,c("Chr", "Start", "End","PeakName")], file = my.out.bckgd.age.BED , sep = "\t" , row.names = F, col.names = F, quote=F)

my.out.fdr5.age.BED <- paste0(my.outprefix,"_AGING_FDR5.bed")
write.table(res.age.annot[res.age.annot$padj < 0.05,c("Chr", "Start", "End","PeakName")], file = my.out.fdr5.age.BED , sep = "\t" , row.names = F, col.names = F, quote=F)

my.out.fdr5.age.BED.up <- paste0(my.outprefix,"_AGING_FDR5_UP_with_Age.bed")
write.table(res.age.annot[bitAnd(res.age.annot$padj < 0.05,res.age.annot$log2FoldChange >0)>0,c("Chr", "Start", "End","PeakName")], file = my.out.fdr5.age.BED.up , sep = "\t" , row.names = F, col.names = F, quote=F)

my.out.fdr5.age.BED.dwn <- paste0(my.outprefix,"_AGING_FDR5_DWN_with_Age.bed")
write.table(res.age.annot[bitAnd(res.age.annot$padj < 0.05,res.age.annot$log2FoldChange <0)>0,c("Chr", "Start", "End","PeakName")], file = my.out.fdr5.age.BED.dwn , sep = "\t" , row.names = F, col.names = F, quote=F)


##### Output sex specific bed
res.age.F     <- cbind(rownames(res.age.F),res.age.F)
res.age.M     <- cbind(rownames(res.age.M),res.age.M)
my.merged.FM  <- cbind(rownames(my.merged.FM),my.merged.FM)
colnames(res.age.F     )[1] <- "PeakName"
colnames(res.age.M     )[1] <- "PeakName"
colnames(my.merged.FM  )[1] <- "PeakName"

res.age.F.annot    <- merge(my.peak.annot[,c("PeakName","Chr", "Start", "End", "Annotation", "Distance.to.TSS", "Gene.Name","Gene.Description" )],res.age.F)
res.age.M.annot    <- merge(my.peak.annot[,c("PeakName","Chr", "Start", "End", "Annotation", "Distance.to.TSS", "Gene.Name","Gene.Description" )],res.age.M)
my.merged.FM.annot <- merge(my.peak.annot[,c("PeakName","Chr", "Start", "End", "Annotation", "Distance.to.TSS", "Gene.Name","Gene.Description" )],my.merged.FM)

####
my.out.stats.age.F <- paste(my.outprefix,"AGING_FEMALE_all_genes_statistics_PeakAnnot.txt",sep = "_")
write.table(res.age.F.annot, file = my.out.stats.age.F , sep = "\t" , row.names = F, quote=F)

my.out.fdr5.age.F <- paste(my.outprefix,"AGING_FEMALE_FDR5_genes_statistics_PeakAnnot.txt",sep = "_")
write.table(res.age.F.annot[res.age.F.annot$padj < 0.05,], file = my.out.fdr5.age.F, sep = "\t" , row.names = F, quote=F)

my.out.bckgd.age.F.BED <- paste0(my.outprefix,"_AGING_FEMALE_Background.bed")
write.table(res.age.F.annot[,c("Chr", "Start", "End","PeakName")], file = my.out.bckgd.age.F.BED , sep = "\t" , row.names = F, col.names = F, quote=F)

my.out.fdr5.age.F.BED.up <- paste0(my.outprefix,"_AGING_FEMALE_FDR5_UP_with_Age.bed")
write.table(res.age.F.annot[bitAnd(res.age.F.annot$padj_F < 0.05,res.age.F.annot$log2FoldChange_F >0)>0,c("Chr", "Start", "End","PeakName")], file = my.out.fdr5.age.F.BED.up , sep = "\t" , row.names = F, col.names = F, quote=F)

my.out.fdr5.age.F.BED.dwn <- paste0(my.outprefix,"_AGING_FEMALE_FDR5_DWN_with_Age.bed")
write.table(res.age.F.annot[bitAnd(res.age.F.annot$padj_F < 0.05,res.age.F.annot$log2FoldChange_F <0)>0,c("Chr", "Start", "End","PeakName")], file = my.out.fdr5.age.F.BED.dwn , sep = "\t" , row.names = F, col.names = F, quote=F)

####
my.out.stats.age.M <- paste(my.outprefix,"AGING_MALE_all_genes_statistics_PeakAnnot.txt",sep = "_")
write.table(res.age.M.annot, file = my.out.stats.age.M , sep = "\t" , row.names = F, quote=F)

my.out.fdr5.age.M <- paste(my.outprefix,"AGING_MALE_FDR5_genes_statistics_PeakAnnot.txt",sep = "_")
write.table(res.age.M.annot[res.age.M.annot$padj_M < 0.05,], file = my.out.fdr5.age.M, sep = "\t" , row.names = F, quote=F)

my.out.bckgd.age.M.BED <- paste0(my.outprefix,"_AGING_MALE_Background.bed")
write.table(res.age.M.annot[,c("Chr", "Start", "End","PeakName")], file = my.out.bckgd.age.M.BED , sep = "\t" , row.names = F, col.names = F, quote=F)

my.out.fdr5.age.M.BED.up <- paste0(my.outprefix,"_AGING_MALE_FDR5_UP_with_Age.bed")
write.table(res.age.M.annot[bitAnd(res.age.M.annot$padj_M < 0.05,res.age.M.annot$log2FoldChange_M >0)>0,c("Chr", "Start", "End","PeakName")], file = my.out.fdr5.age.M.BED.up , sep = "\t" , row.names = F, col.names = F, quote=F)

my.out.fdr5.age.M.BED.dwn <- paste0(my.outprefix,"_AGING_MALE_FDR5_DWN_with_Age.bed")
write.table(res.age.M.annot[bitAnd(res.age.M.annot$padj_M < 0.05,res.age.M.annot$log2FoldChange_M <0)>0,c("Chr", "Start", "End","PeakName")], file = my.out.fdr5.age.M.BED.dwn , sep = "\t" , row.names = F, col.names = F, quote=F)

####
my.out.stats.sexdiv <- paste(my.outprefix,"AGING_Separated_Sex_Merged_Table_all_genes_statistics_PeakAnnot.txt",sep = "_")
write.table(my.merged.FM.annot, file = my.out.stats.sexdiv , sep = "\t" , row.names = F, quote=F)

my.out.bckgd.sexdiv.BED <- paste0(my.outprefix,"_AGING_Separated_Sex_Merged_Table_Background.bed")
write.table(my.merged.FM.annot[,c("Chr", "Start", "End","PeakName")], file = my.out.bckgd.sexdiv.BED , sep = "\t" , row.names = F, col.names = F, quote=F)

#### divergently regulated genes
my.F_notM.tab <- bitAnd(my.merged.FM.annot$padj_F < 0.05, my.merged.FM.annot$padj_M > 0.1) > 0
my.M_notF.tab <- bitAnd(my.merged.FM.annot$padj_M < 0.05, my.merged.FM.annot$padj_F > 0.1) > 0


my.out.fdr5.sexdiv.BED.up.F <- paste0(my.outprefix,"_AGING_FDR5_UP_Female_ONLY_with_Age.bed")
write.table(my.merged.FM.annot[bitAnd(my.F_notM.tab,my.merged.FM.annot$log2FoldChange_F >0)>0,c("Chr", "Start", "End","PeakName")], file = my.out.fdr5.sexdiv.BED.up.F , sep = "\t" , row.names = F, col.names = F, quote=F)

my.out.fdr5.sexdiv.BED.dwn.F <- paste0(my.outprefix,"_AGING_FDR5_DWN_Female_ONLY_with_Age.bed")
write.table(my.merged.FM.annot[bitAnd(my.F_notM.tab,my.merged.FM.annot$log2FoldChange_F <0)>0,c("Chr", "Start", "End","PeakName")], file = my.out.fdr5.sexdiv.BED.dwn.F , sep = "\t" , row.names = F, col.names = F, quote=F)


my.out.fdr5.sexdiv.BED.up.M <- paste0(my.outprefix,"_AGING_FDR5_UP_Male_ONLY_with_Age.bed")
write.table(my.merged.FM.annot[bitAnd(my.M_notF.tab,my.merged.FM.annot$log2FoldChange_M >0)>0,c("Chr", "Start", "End","PeakName")], file = my.out.fdr5.sexdiv.BED.up.M , sep = "\t" , row.names = F, col.names = F, quote=F)

my.out.fdr5.sexdiv.BED.dwn.M <- paste0(my.outprefix,"_AGING_FDR5_DWN_Male_ONLY_with_Age.bed")
write.table(my.merged.FM.annot[bitAnd(my.M_notF.tab,my.merged.FM.annot$log2FoldChange_M <0)>0,c("Chr", "Start", "End","PeakName")], file = my.out.fdr5.sexdiv.BED.dwn.M , sep = "\t" , row.names = F, col.names = F, quote=F)
################################################################################################

#######################
sink(file = paste(my.outprefix,"_session_Info.txt", sep =""))
sessionInfo()
sink()
