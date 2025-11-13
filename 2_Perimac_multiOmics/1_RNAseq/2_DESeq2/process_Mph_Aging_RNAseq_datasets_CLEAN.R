setwd('/Volumes/BB_HQ_1//Immune_sex_dimorphism_Aging/Macrophages/Clean_Macrophage_Data/RNAseq/DEseq2_analysis/DESeq2/')
options(stringsAsFactors = FALSE)

# load libraries for analysis
library(DESeq2)
library(pheatmap)
library('bitops')
library('limma')
library(RColorBrewer)
library(fields)


####################################  peritoneal  Macrophages AC41   ####################################
# read in count matrix
my.mph1 <- read.table('../../STAR/2022-03-30_PeriMac_NIA_aging_CLEAN_counts.txt', sep = "\t", header = T, skip = 1)
my.mph.1 <- my.mph1[,-c(2:6)]
# apply(my.mph.1[,-1],2,sum)

colnames(my.mph.1) <- c("Geneid",
                        "PeriMac_YF1",
                        "PeriMac_YF2",
                        "PeriMac_YF3",
                        "PeriMac_YF4",
                        "PeriMac_YF5",
                        "PeriMac_OF1",
                        "PeriMac_OF2",
                        "PeriMac_OF3",
                        "PeriMac_OF4",
                        "PeriMac_OF5",
                        "PeriMac_YM1",
                        "PeriMac_YM2",
                        "PeriMac_YM3",
                        "PeriMac_YM4",
                        "PeriMac_YM5",
                        "PeriMac_OM1",
                        "PeriMac_OM2",
                        "PeriMac_OM3",
                        "PeriMac_OM4",
                        "PeriMac_OM5")

# process RNAseq data and save RData object
# create covariates for analysis
my.age    <- c(rep(4,5),rep(20,5),rep(4,5),rep(20,5))
my.sex    <- c(rep("F",10),rep("M",10))

# see deseq2 vignette, remove genes without consistent expression
my.keep <- apply(my.mph.1[,-1]> 0, 1, sum) >= dim(my.mph.1)[2]/2

# Now pull out the null/low expressed genes
my.filtered.matrix           <- my.mph.1[my.keep,-1] # 14783
rownames(my.filtered.matrix) <- my.mph.1[my.keep,1]
################################################################################################

############################################################
############# A. model age and sex together ################
############################################################
# get output file prefixes
my.outprefix.age_sex <- paste(Sys.Date(),"Peritoneal_Macrophages_DESeq2_model_with_AGE_SEX",sep="_")

# build design matrix
dataDesign.tg = data.frame(row.names = colnames( my.filtered.matrix ), 
                           age = my.age, # age in months
                           sex = my.sex)

# get matrix using age as a modeling covariate
dds <- DESeqDataSetFromMatrix(countData = my.filtered.matrix,
                              colData = dataDesign.tg,
                              design = ~ age + sex)

# run DESeq normalizations and export results
dds.deseq <- DESeq(dds)

# plot dispersion
my.disp.out <- paste(my.outprefix.age_sex,"_dispersion_plot.pdf")

pdf(my.disp.out)
plotDispEsts(dds.deseq)
dev.off()

# get DESeq2 normalized expression value
tissue.cts <- getVarianceStabilizedData(dds.deseq)

# color-code 
my.colors <- rep("deeppink",dim(my.filtered.matrix)[2])
my.colors[bitAnd(my.age == 20, my.sex %in% "F") > 0] <- "deeppink4"
my.colors[bitAnd(my.age == 4 , my.sex %in% "M") > 0] <- "deepskyblue"
my.colors[bitAnd(my.age == 20, my.sex %in% "M") > 0] <- "deepskyblue4"

# expression range
pdf(paste(my.outprefix.age_sex,"_Normalized_counts_boxplot.pdf"))
boxplot(tissue.cts,col=my.colors,cex=0.5,ylab="Log2 DESeq2 Normalized counts", las = 2)  
dev.off()

# plot Xist expression
pdf(paste(my.outprefix.age_sex,"_Normalized_counts_Xist_expression_barplot.pdf"))
barplot(tissue.cts["Xist",], ylab = "Normalized log2(counts) Xist expression", las = 2, col = my.colors)
dev.off()

pdf(paste(my.outprefix.age_sex,"_Normalized_counts_Uty_expression_barplot.pdf"))
barplot(tissue.cts["Uty",], ylab = "Normalized log2(counts) Uty expression", las = 2, col = my.colors)
dev.off()

# MDS analysis
mds.result <- cmdscale(1-cor(tissue.cts,method="spearman"), k = 2, eig = FALSE, add = FALSE, x.ret = FALSE)
x <- mds.result[, 1]
y <- mds.result[, 2]

pdf(paste0(my.outprefix.age_sex,"_MDS_plot.pdf"))
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

pdf(paste(my.outprefix.age_sex,"_PCA_plot.pdf",sep=""))
plot(x,y, 
     pch = 16, cex=3, col= my.colors, 
     xlab = paste('PC1 (', round(100*my.summary$importance[,1][2],1),"%)", sep=""),
     ylab = paste('PC2 (', round(100*my.summary$importance[,2][2],1),"%)", sep=""),
     cex.lab = 1.5,
     cex.axis = 1.5) 
dev.off()
###########################################################################

##%%%%%%%%%%%%%% 1. aging with sex as covariate  %%%%%%%%%%%%%%
res.age <- results(dds.deseq, name= "age")

### get the heatmap of aging changes at FDR5; exclude NA
res.age <- res.age[!is.na(res.age$padj),]

genes.aging <- rownames(res.age)[res.age$padj < 0.05]
my.num.aging <- length(genes.aging)

if (my.num.aging > 1) {
  # heatmap drawing - only if there is at least 2 genes
  my.heatmap.out <- paste(my.outprefix.age_sex,"AGING_Heatmap_FDR5.pdf", sep = "_")
  
  pdf(my.heatmap.out, onefile = F)
  my.heatmap.title <- paste("Peritoneal_Macrophages"," aging significant (FDR<5%), ", my.num.aging, " genes",sep="")
  pheatmap(tissue.cts[genes.aging,],
           cluster_cols = F,
           cluster_rows = T,
           colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
           show_rownames = F, scale="row",
           main = my.heatmap.title, cellwidth = 15)
  dev.off()
}


#%%%%%%%%%%%%%% 2. sex with age as covariate %%%%%%%%%%%%%%
res.sex <- results(dds.deseq, contrast = c("sex","F","M")) # FC in females over Males

### get the heatmap of sex dimorphic changes at FDR5; exclude NA
res.sex <- res.sex[!is.na(res.sex$padj),]

genes.sex <- rownames(res.sex)[res.sex$padj < 0.05]
my.num.sex <- length(genes.sex)

if (my.num.sex > 1) {
  # heatmap drawing - only if there is at least 2 gene
  my.heatmap.out <- paste(my.outprefix.age_sex,"SEX_DIM_Heatmap_FDR5.pdf", sep = "_")
  
  pdf(my.heatmap.out, onefile = F)
  my.heatmap.title <- paste("Peritoneal_Macrophages"," sex significant (FDR<5%), ", my.num.sex, " genes",sep="")
  pheatmap(tissue.cts[genes.sex,],
           cluster_cols = F,
           cluster_rows = T,
           colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
           show_rownames = F, scale="row",
           main = my.heatmap.title, cellwidth = 15)
  dev.off()
}


# output result tables of combined analysis to text files
my.out.ct.mat <- paste0(my.outprefix.age_sex,"_log2_counts_matrix.txt")
write.table(tissue.cts, file = my.out.ct.mat , sep = "\t" , row.names = T, quote=F)

my.out.stats.age <- paste0(my.outprefix.age_sex,"_AGING_all_genes_statistics.txt")
my.out.stats.sex <- paste0(my.outprefix.age_sex,"_SEX_DIM_all_genes_statistics.txt")
write.table(res.age, file = my.out.stats.age , sep = "\t" , row.names = T, quote=F)
write.table(res.sex, file = my.out.stats.sex , sep = "\t" , row.names = T, quote=F)

my.out.fdr5.age <- paste0(my.outprefix.age_sex,"_AGING_FDR5_genes_statistics.txt")
my.out.fdr5.sex <- paste0(my.outprefix.age_sex,"_SEX_DIM_FDR5_genes_statistics.txt")
write.table(res.age[genes.aging,], file = my.out.fdr5.age, sep = "\t" , row.names = T, quote=F)
write.table(res.sex[genes.sex,], file = my.out.fdr5.sex, sep = "\t" , row.names = T, quote=F)

my.rdata.age <- paste0(my.outprefix.age_sex,"_AGING_DEseq2_object.RData")
my.rdata.sex <- paste0(my.outprefix.age_sex,"_SEX_DIM_DEseq2_object.RData")
save(res.age, file = my.rdata.age)
save(res.sex, file = my.rdata.sex)
################################################################################################


################################################################################################
# get output file prefixes
my.outprefix.age_M <- paste(Sys.Date(),"Peritoneal_Macrophages","DESeq2_model_with_AGING_Males",sep="_")
my.outprefix.age_F <- paste(Sys.Date(),"Peritoneal_Macrophages","DESeq2_model_with_AGING_Females",sep="_")

###################################################################
############# B. model age in each sex separately  ################
###################################################################

# get matrix using age as a modeling covariate in FEMALES
dds.F <- DESeqDataSetFromMatrix(countData = my.filtered.matrix[,my.sex %in% 'F'],
                                colData = dataDesign.tg[my.sex %in% 'F',],
                                design = ~ age)

# run DESeq normalizations and export results
dds.deseq.F <- DESeq(dds.F)
res.age.F <- results(dds.deseq.F, name= "age")

# get matrix using age as a modeling covariate in MALES
dds.M <- DESeqDataSetFromMatrix(countData = my.filtered.matrix[,my.sex %in% 'M'],
                                colData = dataDesign.tg[my.sex %in% 'M',],
                                design = ~ age)

# run DESeq normalizations and export results
dds.deseq.M <- DESeq(dds.M)
res.age.M <- results(dds.deseq.M, name= "age")

my.rdata.age.F <- paste0(my.outprefix.age_sex,"_AGING_Female_DEseq2_object.RData")
my.rdata.age.M <- paste0(my.outprefix.age_sex,"_AGING_Male_DEseq2_object.RData")
save(res.age.F, file = my.rdata.age.F)
save(res.age.M, file = my.rdata.age.M)


#%%%%%%%%%%%%%% 
colnames(res.age.F) <- paste(colnames(res.age.F),"F",sep = "_")
colnames(res.age.M) <- paste(colnames(res.age.M),"M",sep = "_")
my.merged.FM <- cbind(res.age.F,res.age.M)
my.merged.FM <- my.merged.FM[!is.na(my.merged.FM$padj_F),]
my.merged.FM <- my.merged.FM[!is.na(my.merged.FM$padj_M),]

my.spear.cor <- cor.test(my.merged.FM$log2FoldChange_F,my.merged.FM$log2FoldChange_M, method = 'spearman')
my.rho       <- signif(my.spear.cor$estimate,3)

#### commonly regulated genes
my.F_M.5 <- bitAnd(my.merged.FM$padj_F < 0.05, my.merged.FM$padj_M < 0.05) > 0
my.F_M.10 <- bitAnd(my.merged.FM$padj_F < 0.1, my.merged.FM$padj_M < 0.1) > 0

pdf(paste(Sys.Date(),"Peritoneal_Macrophages_aging_male_vs_female_FC_scatterplot_FDR5.pdf", sep = "_"))
smoothScatter(my.merged.FM$log2FoldChange_F,my.merged.FM$log2FoldChange_M, 
              xlim = c(-0.25,0.25), ylim = c(-0.25,0.25),
              xlab = "log2(FC per m) in Females with aging",
              ylab = "log2(FC per m) in Males with aging",
              main = "Peritoneal_Macrophages")
abline(0,1, col = "grey", lty = "dashed")
abline(h = 0, col = "red", lty = "dashed")
abline(v = 0, col = "red", lty = "dashed")
text(-0.26, 0.23, paste("Rho = ",my.rho), pos = 4)
text(-0.26, 0.21, paste("p = ",signif(my.spear.cor$p.value,2)), pos = 4)
points(my.merged.FM$log2FoldChange_F[my.F_M.5],my.merged.FM$log2FoldChange_M[my.F_M.5], cex = 0.5, pch = 16, col = "gold1")
legend("bottomleft",
       c(paste("FDR5%,", sum(my.F_M.5))),
       col = "gold1", pch = 16, pt.cex = 0.5,bty = 'n')
dev.off()

write.table(my.merged.FM[my.F_M.10,], file = paste(Sys.Date(),"Peritoneal_Macrophages","AGING_SexAgreement_FDR10_genes_statistics.txt", sep ="_"), sep = "\t" , row.names = T, quote=F)


#### divergently regulated genes
my.F_notM <- bitAnd(my.merged.FM$padj_F < 0.05, my.merged.FM$padj_M > 0.1) > 0
my.M_notF <- bitAnd(my.merged.FM$padj_M < 0.05, my.merged.FM$padj_F > 0.1) > 0

pdf(paste(Sys.Date(),"Peritoneal_Macrophages_aging_male_vs_female_FC_scatterplot_divergent_FDR5.pdf", sep = "_"))
smoothScatter(my.merged.FM$log2FoldChange_F,my.merged.FM$log2FoldChange_M, 
              xlim = c(-0.25,0.25), ylim = c(-0.25,0.25),
              xlab = "log2(FC per m) in Females with aging",
              ylab = "log2(FC per m) in Males with aging",
              main = "Peritoneal_Macrophages")
abline(0,1, col = "grey", lty = "dashed")
abline(h = 0, col = "red", lty = "dashed")
abline(v = 0, col = "red", lty = "dashed")
text(-0.26, 0.23, paste("Rho = ",my.rho), pos = 4)
text(-0.26, 0.21, paste("p = ",signif(my.spear.cor$p.value,2)), pos = 4)
points(my.merged.FM$log2FoldChange_F[my.F_notM],my.merged.FM$log2FoldChange_M[my.F_notM], cex = 0.5, pch = 16, col = "orchid")
points(my.merged.FM$log2FoldChange_F[my.M_notF],my.merged.FM$log2FoldChange_M[my.M_notF], cex = 0.5, pch = 16, col = "royalblue")
legend("bottomleft",
       c(paste("M not F,", sum(my.M_notF)),paste("F not M,", sum(my.F_notM))),
       col = c("royalblue","orchid"), pch = 16, pt.cex = 0.5,bty = 'n')
dev.off()

write.table(my.merged.FM[my.M_notF,], file = paste(Sys.Date(),"Peritoneal_Macrophages_AGING_Male.NOT.Female_FDR5_genes_statistics.txt", sep ="_"), sep = "\t" , row.names = T, quote=F)
write.table(my.merged.FM[my.F_notM,], file = paste(Sys.Date(),"Peritoneal_Macrophages_AGING_Female.NOT.Male_FDR5_genes_statistics.txt", sep ="_"), sep = "\t" , row.names = T, quote=F)

write.table(my.merged.FM, file = paste(Sys.Date(),"Peritoneal_Macrophages_AGING_Separated_Sex_Merged_Table_ALL_genes_statistics.txt", sep ="_"), sep = "\t" , row.names = T, quote=F)
write.table(res.age.F, file = paste(Sys.Date(),"Peritoneal_Macrophages_AGING_Female_ALL_genes_statistics.txt", sep ="_"), sep = "\t" , row.names = T, quote=F)
write.table(res.age.M, file = paste(Sys.Date(),"Peritoneal_Macrophages_AGING_Male_ALL_genes_statistics.txt", sep ="_"), sep = "\t" , row.names = T, quote=F)

################################################################################################
################################################################################################

my.mph1.RNAseq.process <- list("Sex" = res.sex,
                               "Aging" = res.age,
                               "Female_Aging" = res.age.F,
                               "Male_Aging" = res.age.M,
                               "FemaleMaleAging" = my.merged.FM)


save(my.mph1.RNAseq.process, file=paste(Sys.Date(),"Peritoneal_Macrophages_RNA_seq_results.RData",sep = "_"))
#####################################################################################


#######################
sink(file = paste(Sys.Date(),"Macrophage_aging_sex_RNAseq_analysis_session_Info.txt", sep =""))
sessionInfo()
sink()

