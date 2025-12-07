setwd('/Volumes/BB_HQ_1/Immune_sex_dimorphism_Aging/Macrophages/Clean_Macrophage_Data/Peritoneal_scRNAseq/Combined_Analysis/Muscat_PB/')
options(stringsAsFactors=FALSE)

library(ggplot2)          # 
library(scales)           # 
library(ComplexHeatmap)   #
library(circlize)         #
library(bitops)         #

library(metaRNASeq)   #
library(Vennerable)


theme_set(theme_bw())   

# 2025-08-27
# combine cohorts/datasets to analyze with muscat for main cell types

###############################################################################################
# 0. Load DESeq2 results

load("/Volumes/BB_HQ_1/Immune_sex_dimorphism_Aging/Macrophages/Clean_Macrophage_Data/Peritoneal_scRNAseq/v2_set/Muscat_PB/2025-08-28_pseudobulk_NIAv2_MergedFM_DEseq2_objects.RData")
NIAv2.FM <- deseq.res.list.perit.FM; rm(deseq.res.list.perit.FM);

load("/Volumes/BB_HQ_1/Immune_sex_dimorphism_Aging/Macrophages/Clean_Macrophage_Data/Peritoneal_scRNAseq/v3_NIA/Muscat_PB/2025-08-28_pseudobulk_NIAv3_MergedFM_DEseq2_objects.RData")
NIAv3.FM <- deseq.res.list.perit.FM; rm(deseq.res.list.perit.FM);

load("/Volumes/BB_HQ_1/Immune_sex_dimorphism_Aging/Macrophages/Clean_Macrophage_Data/Peritoneal_scRNAseq/v3_set/Muscat_PB/2025-08-28_pseudobulk_JAXv3_MergedFM_DEseq2_objects.RData")
JAXv3.FM <- deseq.res.list.perit.FM; rm(deseq.res.list.perit.FM);
###############################################################################################

###############################################################################################
# 1. Get meta analysis

### Can only do cell types found in all 3 datasets

# get common genes
bcell.com.genes <- intersect(intersect(NIAv2.FM$B_cells$Row.names,NIAv3.FM$B_cells$Row.names),JAXv3.FM$B_cells$Row.names)
mph.com.genes   <- intersect(intersect(NIAv2.FM$Macrophages$Row.names,NIAv3.FM$Macrophages$Row.names),JAXv3.FM$Macrophages$Row.names)
tcell.com.genes <- intersect(intersect(NIAv2.FM$T_cells$Row.names,NIAv3.FM$T_cells$Row.names),JAXv3.FM$T_cells$Row.names)

# put rownames
rownames(NIAv2.FM$B_cells    ) <- NIAv2.FM$B_cells    $Row.names
rownames(NIAv3.FM$B_cells    ) <- NIAv3.FM$B_cells    $Row.names
rownames(JAXv3.FM$B_cells    ) <- JAXv3.FM$B_cells    $Row.names
rownames(NIAv2.FM$Macrophages) <- NIAv2.FM$Macrophages$Row.names
rownames(NIAv3.FM$Macrophages) <- NIAv3.FM$Macrophages$Row.names
rownames(JAXv3.FM$Macrophages) <- JAXv3.FM$Macrophages$Row.names
rownames(NIAv2.FM$T_cells    ) <- NIAv2.FM$T_cells    $Row.names
rownames(NIAv3.FM$T_cells    ) <- NIAv3.FM$T_cells    $Row.names
rownames(JAXv3.FM$T_cells    ) <- JAXv3.FM$T_cells    $Row.names

# restrict to common genes across datasets
NIAv2.FM$B_cells     <- NIAv2.FM$B_cells    [bcell.com.genes,]
NIAv3.FM$B_cells     <- NIAv3.FM$B_cells    [bcell.com.genes,]
JAXv3.FM$B_cells     <- JAXv3.FM$B_cells    [bcell.com.genes,]
NIAv2.FM$Macrophages <- NIAv2.FM$Macrophages[mph.com.genes,]
NIAv3.FM$Macrophages <- NIAv3.FM$Macrophages[mph.com.genes,]
JAXv3.FM$Macrophages <- JAXv3.FM$Macrophages[mph.com.genes,]
NIAv2.FM$T_cells     <- NIAv2.FM$T_cells    [tcell.com.genes,]
NIAv3.FM$T_cells     <- NIAv3.FM$T_cells    [tcell.com.genes,]
JAXv3.FM$T_cells     <- JAXv3.FM$T_cells    [tcell.com.genes,]

my.studies      <- c("NIAv2", "NIAv3", "JAXv3")



############################################################
##########  I. B cells

################ a. Females
# We recommand to store both p-value and Fold Change results in lists in order to perform
# meta-analysis and keep track of the potential conflicts (see section 5)
F.bcell.rawpval <- list("pval.NIAv2" = NIAv2.FM$B_cells$pvalue_F        , "pval.NIAv3" = NIAv3.FM$B_cells$pvalue_F        , "pval.JAXv3" = JAXv3.FM$B_cells$pvalue_F)
F.bcell.FC      <- list("FC.NIAv2"   = NIAv2.FM$B_cells$log2FoldChange_F, "FC.NIAv3"   = NIAv3.FM$B_cells$log2FoldChange_F, "FC.JAXv3"   = JAXv3.FM$B_cells$log2FoldChange_F)

# Differentially expressed genes in each individual study can also be marked in a matrix DE:
F.bcell.padj <- list("padj.NIAv2" = NIAv2.FM$B_cells$padj_F        , "padj.NIAv3" = NIAv3.FM$B_cells$padj_F        , "padj.JAXv3" = JAXv3.FM$B_cells$padj_F)

F.bcell.DE           <- data.frame(mapply(F.bcell.padj, FUN=function(x) ifelse(x <= 0.05, 1, 0)))
colnames(F.bcell.DE) <- paste("DE", my.studies, sep=".")
rownames(F.bcell.DE)     <- rownames(NIAv2.FM$B_cells)


pdf(paste0(Sys.Date(), "_rawp_histogram_Bcells_F.pdf"), height = 5, width = 10)
par(mfrow = c(1,3))
hist(F.bcell.rawpval[[1]], breaks=100, col="grey", main="B cell - F- NIAv2", xlab="Raw p-value")
hist(F.bcell.rawpval[[2]], breaks=100, col="grey", main="B cell - F- NIAv3", xlab="Raw p-value")
hist(F.bcell.rawpval[[3]], breaks=100, col="grey", main="B cell - F- JAXv3", xlab="Raw p-value")
par(mfrow = c(1,1))
dev.off()

F.bcell.invnormcomb <- invnorm(F.bcell.rawpval, nrep=c(4,10,12), BHth = 0.05)

pdf(paste0(Sys.Date(), "_META_pvalue_histogram_Bcell_Faging.pdf"), height = 5, width = 5)
hist(F.bcell.invnormcomb$rawpval, breaks=100, col="grey", main="B cell - F - Inverse normal method", xlab = "Raw p-values (meta-analysis)")
dev.off()

# We build a matrix signsFC gathering all signs of fold changes from individual studies.
F.bcell.signsFC      <- mapply(F.bcell.FC, FUN=function(x) sign(x))
F.bcell.sumsigns     <- apply(F.bcell.signsFC,1,sum)
F.bcell.commonsgnFC  <- ifelse(abs(F.bcell.sumsigns)==dim(F.bcell.signsFC)[2], sign(F.bcell.sumsigns),0)

F.bcell.res <- cbind( F.bcell.DE,
                      "meta_sign"   = F.bcell.commonsgnFC,
                      "InvNorm_FDR" = F.bcell.invnormcomb$adjpval)


################ b. Males
# We recommand to store both p-value and Fold Change results in lists in order to perform
# meta-analysis and keep track of the potential conflicts (see section 5)
M.bcell.rawpval <- list("pval.NIAv2" = NIAv2.FM$B_cells$pvalue_M        , "pval.NIAv3" = NIAv3.FM$B_cells$pvalue_M        , "pval.JAXv3" = JAXv3.FM$B_cells$pvalue_M)
M.bcell.FC      <- list("FC.NIAv2"   = NIAv2.FM$B_cells$log2FoldChange_M, "FC.NIAv3"   = NIAv3.FM$B_cells$log2FoldChange_M, "FC.JAXv3"   = JAXv3.FM$B_cells$log2FoldChange_M)

# Differentially expressed genes in each individual study can also be marked in a matrix DE:
M.bcell.padj <- list("padj.NIAv2" = NIAv2.FM$B_cells$padj_M        , "padj.NIAv3" = NIAv3.FM$B_cells$padj_M        , "padj.JAXv3" = JAXv3.FM$B_cells$padj_M)

M.bcell.DE           <- data.frame(mapply(M.bcell.padj, FUN=function(x) ifelse(x <= 0.05, 1, 0)))
colnames(M.bcell.DE) <- paste("DE", my.studies, sep=".")
rownames(M.bcell.DE)     <- rownames(NIAv2.FM$B_cells)


pdf(paste0(Sys.Date(), "_rawp_histogram_Bcells_M.pdf"), height = 5, width = 10)
par(mfrow = c(1,3))
hist(M.bcell.rawpval[[1]], breaks=100, col="grey", main="B cell - M - NIAv2", xlab="Raw p-value")
hist(M.bcell.rawpval[[2]], breaks=100, col="grey", main="B cell - M - NIAv3", xlab="Raw p-value")
hist(M.bcell.rawpval[[3]], breaks=100, col="grey", main="B cell - M - JAXv3", xlab="Raw p-value")
par(mfrow = c(1,1))
dev.off()

M.bcell.invnormcomb <- invnorm(M.bcell.rawpval, nrep=c(4,10,12), BHth = 0.05)

pdf(paste0(Sys.Date(), "_META_pvalue_histogram_Bcell_Maging.pdf"), height = 5, width = 5)
hist(M.bcell.invnormcomb$rawpval, breaks=100, col="grey", main="B cell - M - Inverse normal method", xlab = "Raw p-values (meta-analysis)")
dev.off()

# We build a matrix signsFC gathering all signs of fold changes from individual studies.
M.bcell.signsFC      <- mapply(M.bcell.FC, FUN=function(x) sign(x))
M.bcell.sumsigns     <- apply(M.bcell.signsFC,1,sum)
M.bcell.commonsgnFC  <- ifelse(abs(M.bcell.sumsigns)==dim(M.bcell.signsFC)[2], sign(M.bcell.sumsigns),0)

M.bcell.res <- cbind( M.bcell.DE,
                      "meta_sign"   = M.bcell.commonsgnFC,
                      "InvNorm_FDR" = M.bcell.invnormcomb$adjpval)


############################################################
##########  II. Macrophages

################ a. Females
# We recommand to store both p-value and Fold Change results in lists in order to perform
# meta-analysis and keep track of the potential conflicts (see section 5)
F.mph.rawpval <- list("pval.NIAv2" = NIAv2.FM$Macrophages$pvalue_F        , "pval.NIAv3" = NIAv3.FM$Macrophages$pvalue_F        , "pval.JAXv3" = JAXv3.FM$Macrophages$pvalue_F)
F.mph.FC      <- list("FC.NIAv2"   = NIAv2.FM$Macrophages$log2FoldChange_F, "FC.NIAv3"   = NIAv3.FM$Macrophages$log2FoldChange_F, "FC.JAXv3"   = JAXv3.FM$Macrophages$log2FoldChange_F)

# Differentially expressed genes in each individual study can also be marked in a matrix DE:
F.mph.padj <- list("padj.NIAv2" = NIAv2.FM$Macrophages$padj_F        , "padj.NIAv3" = NIAv3.FM$Macrophages$padj_F        , "padj.JAXv3" = JAXv3.FM$Macrophages$padj_F)

F.mph.DE           <- data.frame(mapply(F.mph.padj, FUN=function(x) ifelse(x <= 0.05, 1, 0)))
colnames(F.mph.DE) <- paste("DE", my.studies, sep=".")
rownames(F.mph.DE)     <- rownames(NIAv2.FM$Macrophages)


pdf(paste0(Sys.Date(), "_rawp_histogram_Macrophages_F.pdf"), height = 5, width = 10)
par(mfrow = c(1,3))
hist(F.mph.rawpval[[1]], breaks=100, col="grey", main="Macrophages - F - NIAv2", xlab="Raw p-value")
hist(F.mph.rawpval[[2]], breaks=100, col="grey", main="Macrophages - F - NIAv3", xlab="Raw p-value")
hist(F.mph.rawpval[[3]], breaks=100, col="grey", main="Macrophages - F - JAXv3", xlab="Raw p-value")
par(mfrow = c(1,1))
dev.off()

F.mph.invnormcomb <- invnorm(F.mph.rawpval, nrep=c(4,10,12), BHth = 0.05)

pdf(paste0(Sys.Date(), "_META_pvalue_histogram_mph_Faging.pdf"), height = 5, width = 5)
hist(F.mph.invnormcomb$rawpval, breaks=100, col="grey", main="Macrophages - F - Inverse normal method", xlab = "Raw p-values (meta-analysis)")
dev.off()

# We build a matrix signsFC gathering all signs of fold changes from individual studies.
F.mph.signsFC      <- mapply(F.mph.FC, FUN=function(x) sign(x))
F.mph.sumsigns     <- apply(F.mph.signsFC,1,sum)
F.mph.commonsgnFC  <- ifelse(abs(F.mph.sumsigns)==dim(F.mph.signsFC)[2], sign(F.mph.sumsigns),0)

F.mph.res <- cbind( F.mph.DE,
                    "meta_sign"   = F.mph.commonsgnFC,
                    "InvNorm_FDR" = F.mph.invnormcomb$adjpval)


################ b. Males
# We recommand to store both p-value and Fold Change results in lists in order to perform
# meta-analysis and keep track of the potential conflicts (see section 5)
M.mph.rawpval <- list("pval.NIAv2" = NIAv2.FM$Macrophages$pvalue_M        , "pval.NIAv3" = NIAv3.FM$Macrophages$pvalue_M        , "pval.JAXv3" = JAXv3.FM$Macrophages$pvalue_M)
M.mph.FC      <- list("FC.NIAv2"   = NIAv2.FM$Macrophages$log2FoldChange_M, "FC.NIAv3"   = NIAv3.FM$Macrophages$log2FoldChange_M, "FC.JAXv3"   = JAXv3.FM$Macrophages$log2FoldChange_M)

# Differentially expressed genes in each individual study can also be marked in a matrix DE:
M.mph.padj <- list("padj.NIAv2" = NIAv2.FM$Macrophages$padj_M        , "padj.NIAv3" = NIAv3.FM$Macrophages$padj_M        , "padj.JAXv3" = JAXv3.FM$Macrophages$padj_M)

M.mph.DE           <- data.frame(mapply(M.mph.padj, FUN=function(x) ifelse(x <= 0.05, 1, 0)))
colnames(M.mph.DE) <- paste("DE", my.studies, sep=".")
rownames(M.mph.DE)     <- rownames(NIAv2.FM$Macrophages)


pdf(paste0(Sys.Date(), "_rawp_histogram_Macrophages_M.pdf"), height = 5, width = 10)
par(mfrow = c(1,3))
hist(M.mph.rawpval[[1]], breaks=100, col="grey", main="Macrophages - M - NIAv2", xlab="Raw p-value")
hist(M.mph.rawpval[[2]], breaks=100, col="grey", main="Macrophages - M - NIAv3", xlab="Raw p-value")
hist(M.mph.rawpval[[3]], breaks=100, col="grey", main="Macrophages - M - JAXv3", xlab="Raw p-value")
par(mfrow = c(1,1))
dev.off()

M.mph.invnormcomb <- invnorm(M.mph.rawpval, nrep=c(4,10,12), BHth = 0.05)

pdf(paste0(Sys.Date(), "_META_pvalue_histogram_mph_Maging.pdf"), height = 5, width = 5)
hist(M.mph.invnormcomb$rawpval, breaks=100, col="grey", main="Macrophages - M - Inverse normal method", xlab = "Raw p-values (meta-analysis)")
dev.off()

# We build a matrix signsFC gathering all signs of fold changes from individual studies.
M.mph.signsFC      <- mapply(M.mph.FC, FUN=function(x) sign(x))
M.mph.sumsigns     <- apply(M.mph.signsFC,1,sum)
M.mph.commonsgnFC  <- ifelse(abs(M.mph.sumsigns)==dim(M.mph.signsFC)[2], sign(M.mph.sumsigns),0)

M.mph.res <- cbind( M.mph.DE,
                    "meta_sign"   = M.mph.commonsgnFC,
                    "InvNorm_FDR" = M.mph.invnormcomb$adjpval)


############################################################
##########  III. T cells

################ a. Females
# We recommand to store both p-value and Fold Change results in lists in order to perform
# meta-analysis and keep track of the potential conflicts (see section 5)
F.tcell.rawpval <- list("pval.NIAv2" = NIAv2.FM$T_cells$pvalue_F        , "pval.NIAv3" = NIAv3.FM$T_cells$pvalue_F        , "pval.JAXv3" = JAXv3.FM$T_cells$pvalue_F)
F.tcell.FC      <- list("FC.NIAv2"   = NIAv2.FM$T_cells$log2FoldChange_F, "FC.NIAv3"   = NIAv3.FM$T_cells$log2FoldChange_F, "FC.JAXv3"   = JAXv3.FM$T_cells$log2FoldChange_F)

# Differentially expressed genes in each individual study can also be marked in a matrix DE:
F.tcell.padj <- list("padj.NIAv2" = NIAv2.FM$T_cells$padj_F        , "padj.NIAv3" = NIAv3.FM$T_cells$padj_F        , "padj.JAXv3" = JAXv3.FM$T_cells$padj_F)

F.tcell.DE           <- data.frame(mapply(F.tcell.padj, FUN=function(x) ifelse(x <= 0.05, 1, 0)))
colnames(F.tcell.DE) <- paste("DE", my.studies, sep=".")
rownames(F.tcell.DE)     <- rownames(NIAv2.FM$T_cells)


pdf(paste0(Sys.Date(), "_rawp_histogram_tcells_F.pdf"), height = 5, width = 10)
par(mfrow = c(1,3))
hist(F.tcell.rawpval[[1]], breaks=100, col="grey", main="T cell - F- NIAv2", xlab="Raw p-value")
hist(F.tcell.rawpval[[2]], breaks=100, col="grey", main="T cell - F- NIAv3", xlab="Raw p-value")
hist(F.tcell.rawpval[[3]], breaks=100, col="grey", main="T cell - F- JAXv3", xlab="Raw p-value")
par(mfrow = c(1,1))
dev.off()

F.tcell.invnormcomb <- invnorm(F.tcell.rawpval, nrep=c(4,10,12), BHth = 0.05)

pdf(paste0(Sys.Date(), "_META_pvalue_histogram_tcell_Faging.pdf"), height = 5, width = 5)
hist(F.tcell.invnormcomb$rawpval, breaks=100, col="grey", main="T cell - F - Inverse normal method", xlab = "Raw p-values (meta-analysis)")
dev.off()

# We build a matrix signsFC gathering all signs of fold changes from individual studies.
F.tcell.signsFC      <- mapply(F.tcell.FC, FUN=function(x) sign(x))
F.tcell.sumsigns     <- apply(F.tcell.signsFC,1,sum)
F.tcell.commonsgnFC  <- ifelse(abs(F.tcell.sumsigns)==dim(F.tcell.signsFC)[2], sign(F.tcell.sumsigns),0)

F.tcell.res <- cbind( F.tcell.DE,
                      "meta_sign"   = F.tcell.commonsgnFC,
                      "InvNorm_FDR" = F.tcell.invnormcomb$adjpval)


################ b. Males
# We recommand to store both p-value and Fold Change results in lists in order to perform
# meta-analysis and keep track of the potential conflicts (see section 5)
M.tcell.rawpval <- list("pval.NIAv2" = NIAv2.FM$T_cells$pvalue_M        , "pval.NIAv3" = NIAv3.FM$T_cells$pvalue_M        , "pval.JAXv3" = JAXv3.FM$T_cells$pvalue_M)
M.tcell.FC      <- list("FC.NIAv2"   = NIAv2.FM$T_cells$log2FoldChange_M, "FC.NIAv3"   = NIAv3.FM$T_cells$log2FoldChange_M, "FC.JAXv3"   = JAXv3.FM$T_cells$log2FoldChange_M)

# Differentially expressed genes in each individual study can also be marked in a matrix DE:
M.tcell.padj <- list("padj.NIAv2" = NIAv2.FM$T_cells$padj_M        , "padj.NIAv3" = NIAv3.FM$T_cells$padj_M        , "padj.JAXv3" = JAXv3.FM$T_cells$padj_M)

M.tcell.DE           <- data.frame(mapply(M.tcell.padj, FUN=function(x) ifelse(x <= 0.05, 1, 0)))
colnames(M.tcell.DE) <- paste("DE", my.studies, sep=".")
rownames(M.tcell.DE)     <- rownames(NIAv2.FM$T_cells)


pdf(paste0(Sys.Date(), "_rawp_histogram_tcells_M.pdf"), height = 5, width = 10)
par(mfrow = c(1,3))
hist(M.tcell.rawpval[[1]], breaks=100, col="grey", main="T cell - M - NIAv2", xlab="Raw p-value")
hist(M.tcell.rawpval[[2]], breaks=100, col="grey", main="T cell - M - NIAv3", xlab="Raw p-value")
hist(M.tcell.rawpval[[3]], breaks=100, col="grey", main="T cell - M - JAXv3", xlab="Raw p-value")
par(mfrow = c(1,1))
dev.off()

M.tcell.invnormcomb <- invnorm(M.tcell.rawpval, nrep=c(4,10,12), BHth = 0.05)

pdf(paste0(Sys.Date(), "_META_pvalue_histogram_tcell_Maging.pdf"), height = 5, width = 5)
hist(M.tcell.invnormcomb$rawpval, breaks=100, col="grey", main="T cell - M - Inverse normal method", xlab = "Raw p-values (meta-analysis)")
dev.off()

# We build a matrix signsFC gathering all signs of fold changes from individual studies.
M.tcell.signsFC      <- mapply(M.tcell.FC, FUN=function(x) sign(x))
M.tcell.sumsigns     <- apply(M.tcell.signsFC,1,sum)
M.tcell.commonsgnFC  <- ifelse(abs(M.tcell.sumsigns)==dim(M.tcell.signsFC)[2], sign(M.tcell.sumsigns),0)

M.tcell.res <- cbind( M.tcell.DE,
                      "meta_sign"   = M.tcell.commonsgnFC,
                      "InvNorm_FDR" = M.tcell.invnormcomb$adjpval)
##########################################################################################################################################

##########################################################################################################################################
# 2. filter metaanalysis results

F.bcell.res.flt  <- F.bcell.res [ bitAnd(F.bcell.res$meta_sign !=0, F.bcell.res$InvNorm_FDR < 0.05)>0,   ] # 391
M.bcell.res.flt  <- M.bcell.res [ bitAnd(M.bcell.res$meta_sign !=0, M.bcell.res$InvNorm_FDR < 0.05)>0,   ] # 285
F.mph.res.flt    <- F.mph.res   [ bitAnd(F.mph.res  $meta_sign !=0, F.mph.res  $InvNorm_FDR < 0.05)>0,   ] # 141
M.mph.res.flt    <- M.mph.res   [ bitAnd(M.mph.res  $meta_sign !=0, M.mph.res  $InvNorm_FDR < 0.05)>0,   ] # 246
F.tcell.res.flt  <- F.tcell.res [ bitAnd(F.tcell.res$meta_sign !=0, F.tcell.res$InvNorm_FDR < 0.05)>0,   ] # 12
M.tcell.res.flt  <- M.tcell.res [ bitAnd(M.tcell.res$meta_sign !=0, M.tcell.res$InvNorm_FDR < 0.05)>0,   ] # 16

### export to excel
options(java.parameters = "-Xmx16g" )
require(openxlsx)

write.xlsx(list("Bcell_F"      = F.bcell.res.flt, 
                "Bcell_M"      = M.bcell.res.flt,
                "Macrophage_F" = F.mph.res.flt  ,
                "Macrophage_M" = M.mph.res.flt  ,
                "Tcell_F"      = F.tcell.res.flt,
                "Tcell_M"      = M.tcell.res.flt  ), 
           rowNames = TRUE, file = paste0(Sys.Date(),"_Peritoneal_Cells_PB_Aging_by_Sex_DESeq2_Results_METARNASEQ.xlsx"))


### get overall
Bcell.all <- list("B cell F" = rownames(F.bcell.res.flt),
                  "B cell M" = rownames(M.bcell.res.flt))
Mph.all   <- list("Macrophage F" = rownames(F.mph.res.flt),
                  "Macrophage M" = rownames(M.mph.res.flt))
tcell.all <- list("T cell F" = rownames(F.tcell.res.flt),
                  "T cell M" = rownames(M.tcell.res.flt))

####### a. B cell up
Bcell.up <- list("B cell F up" = rownames(F.bcell.res.flt)[F.bcell.res.flt$meta_sign > 0],
                 "B cell M up" = rownames(M.bcell.res.flt)[M.bcell.res.flt$meta_sign > 0])
lapply(Bcell.up,length)
# $`B cell F up`
# [1] 263
# $`B cell M up`
# [1] 160

my.Venn <- Venn(Bcell.up)

pdf(paste0(Sys.Date(),"_META_Bcell_UP.pdf"))
plot(my.Venn, doWeights = TRUE, type = "circles",  show = list(Faces = FALSE))
dev.off()

####### b. B cell down
Bcell.dwn <- list("B cell F down" = rownames(F.bcell.res.flt)[F.bcell.res.flt$meta_sign < 0],
                  "B cell M down" = rownames(M.bcell.res.flt)[M.bcell.res.flt$meta_sign < 0])
lapply(Bcell.dwn,length)
# $`B cell F down`
# [1] 128
# 
# $`B cell M down`
# [1] 125

my.Venn <- Venn(Bcell.dwn)

pdf(paste0(Sys.Date(),"_META_Bcell_Down.pdf"))
plot(my.Venn, doWeights = TRUE, type = "circles",  show = list(Faces = FALSE))
dev.off()




####### c. Macrophage up
Mph.up <- list("Macrophage F up" = rownames(F.mph.res.flt)[F.mph.res.flt$meta_sign > 0],
               "Macrophage M up" = rownames(M.mph.res.flt)[M.mph.res.flt$meta_sign > 0])
lapply(Mph.up,length)
# $`Macrophage F up`
# [1] 107
# $`Macrophage M up`
# [1] 187

my.Venn <- Venn(Mph.up)

pdf(paste0(Sys.Date(),"_META_macrophage_UP.pdf"))
plot(my.Venn, doWeights = TRUE, type = "circles",  show = list(Faces = FALSE))
dev.off()

####### d. Macrophage dwon
Mph.dwn <- list("Macrophage F down" = rownames(F.mph.res.flt)[F.mph.res.flt$meta_sign < 0],
                "Macrophage M down" = rownames(M.mph.res.flt)[M.mph.res.flt$meta_sign < 0])
lapply(Mph.dwn,length)
# $`Macrophage F down`
# [1] 34
# $`Macrophage M down`
# [1] 59

my.Venn <- Venn(Mph.dwn)

pdf(paste0(Sys.Date(),"_META_macrophage_DWN.pdf"))
plot(my.Venn, doWeights = TRUE, type = "circles",  show = list(Faces = FALSE))
dev.off()

####### e. T cell up
tcell.up <- list("T cell F up" = rownames(F.tcell.res.flt)[F.tcell.res.flt$meta_sign > 0],
                 "T cell M up" = rownames(M.tcell.res.flt)[M.tcell.res.flt$meta_sign > 0])
lapply(tcell.up,length)
# $`T cell F up`
# [1] 7
# $`T cell M up`
# [1] 14

my.Venn <- Venn(tcell.up)

pdf(paste0(Sys.Date(),"_META_Tcell_UP.pdf"))
plot(my.Venn, doWeights = TRUE, type = "circles",  show = list(Faces = FALSE))
dev.off()

####### f. T cell dwon
tcell.dwn <- list("T cell F down" = rownames(F.tcell.res.flt)[F.tcell.res.flt$meta_sign < 0],
                  "T cell M down" = rownames(M.tcell.res.flt)[M.tcell.res.flt$meta_sign < 0])
lapply(tcell.dwn,length)
# $`T cell F down`
# [1] 5
# $`T cell M down`
# [1] 2

my.Venn <- Venn(tcell.dwn)

pdf(paste0(Sys.Date(),"_META_Tcell_DWN.pdf"))
plot(my.Venn, doWeights = TRUE, type = "circles",  show = list(Faces = FALSE))
dev.off()



#################################################
# Function for computing Jaccard Similarity
jaccard_similarity <- function(my.list) {
  intersection = length(intersect(my.list[[1]], my.list[[2]]))
  union = length(my.list[[1]]) + length( my.list[[2]]) - intersection
  return (intersection/union)
}
#################################################

# jaccard_similarity(Bcell.up)  # 0.2896341
# jaccard_similarity(Bcell.dwn) # 0.3177083
# jaccard_similarity(Mph.up)    # 0.2838428
# jaccard_similarity(Mph.dwn)   # 0.2077922
# jaccard_similarity(tcell.up)  # 0.5
# jaccard_similarity(tcell.dwn) # 0.4

my.res.table.j           <- data.frame(matrix(0,3,3))
colnames(my.res.table.j) <- c("Upregulated","Downregulated", "Overall")
rownames(my.res.table.j) <- c("B_cells","Macrophages","T_cells")

my.res.table.j["B_cells"    ,"Upregulated"  ] <- jaccard_similarity(Bcell.up) 
my.res.table.j["B_cells"    ,"Downregulated"] <- jaccard_similarity(Bcell.dwn)
my.res.table.j["B_cells"    ,"Overall"]       <- jaccard_similarity(Bcell.all)
my.res.table.j["Macrophages","Upregulated"  ] <- jaccard_similarity(Mph.up)   
my.res.table.j["Macrophages","Downregulated"] <- jaccard_similarity(Mph.dwn)  
my.res.table.j["Macrophages","Overall"]       <- jaccard_similarity(Mph.all)  
my.res.table.j["T_cells"    ,"Upregulated"  ] <- jaccard_similarity(tcell.up) 
my.res.table.j["T_cells"    ,"Downregulated"] <- jaccard_similarity(tcell.dwn)
my.res.table.j["T_cells"    ,"Overall"]       <- jaccard_similarity(tcell.all)
 


pdf(paste0(Sys.Date(), "_META_Jaccard_Female_Male_Aging_snRNAseq_v2.pdf"), height = 5, width = 5)
barplot(t(as.matrix(my.res.table.j[3:1,])), 
        beside = T, horiz = T, col = c("firebrick2", "dodgerblue4", "grey30"),
        xlab = "Jaccard Index for DEGs", las = 1,
        main = "Shared Aging DEGs between sexes")
box()
abline(v = 0, col = "grey", lty = "dashed")
abline(v = 0.5, col = "grey", lty = "dashed")
dev.off()

#######################
sink(file = paste(Sys.Date(),"_scRNAseq_peritoneal_3Cohorts_COMBINED_session_Info.txt", sep =""))
sessionInfo()
sink()

