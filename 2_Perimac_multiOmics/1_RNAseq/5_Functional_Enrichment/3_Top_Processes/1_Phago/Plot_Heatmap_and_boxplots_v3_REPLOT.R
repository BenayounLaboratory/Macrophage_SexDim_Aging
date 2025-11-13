setwd('/Volumes/BB_HQ_1/Immune_sex_dimorphism_Aging/Macrophages/Clean_Macrophage_Data/RNAseq/Functional_enrichment/Top_process_for_validation/Phagocytosis_genes/noSVA')
library('grDevices')
library('bitops')
library('pheatmap')


#### load expression data
tissue.cts <- read.table('../../../../DEseq2_analysis/DEseq2/2022-03-30_Peritoneal_Macrophages_DESeq2_model_with_AGE_SEX_log2_counts_matrix.txt', sep = "\t", header = T)

# load gene lists
my.bassik.phago  <- read.table('/Volumes/BB_HQ_1/Immune_sex_dimorphism_Aging/Macrophages/Macrophage_datasets_for_Comparison/Mike_Bassik_Phagocytosis/2021-09-09_Mouse_homologs_Ens104.txt', sep = "\t", header = T)
my.wyss.phago    <- read.table('/Volumes/BB_HQ_1/Immune_sex_dimorphism_Aging/Macrophages/Macrophage_datasets_for_Comparison/Wyss_Corray_Screen/2022-01-07_Biomart_ens105_Name_conversion_FDR5_hits_BV2_screen_286.txt', sep = "\t", header = T)

my.var.null <- apply(tissue.cts,1,var)

my.bassik.phago.v2  <- intersect(unique(my.bassik.phago$Mouse.gene.name) , rownames(tissue.cts)[my.var.null != 0])
my.wyss.phago.v2    <- intersect(unique(my.wyss.phago$Gene.name) , rownames(tissue.cts)[my.var.null != 0])
my.CIRSPR.phago     <- unique(c(my.bassik.phago.v2,my.wyss.phago.v2))


################################################################################
fem.spe <- read.table('../../../../DEseq2_analysis/DEseq2/2022-03-30_Peritoneal_Macrophages_AGING_Female.NOT.Male_FDR5_genes_statistics.txt')
mal.spe <- read.table('../../../../DEseq2_analysis/DEseq2/2022-03-30_Peritoneal_Macrophages_AGING_Male.NOT.Female_FDR5_genes_statistics.txt')

fem.spe.up  <- rownames(fem.spe)[fem.spe$log2FoldChange_F > 0]
fem.spe.dwn <- rownames(fem.spe)[fem.spe$log2FoldChange_F < 0]

mal.spe.up  <- rownames(mal.spe)[mal.spe$log2FoldChange_F > 0]
mal.spe.dwn <- rownames(mal.spe)[mal.spe$log2FoldChange_F < 0]

library(Vennerable)

###########################
##### Female analysis
my.criteria <- list("Female Aging Up"    = fem.spe.up,
                    "Female Aging Dwn"   = fem.spe.dwn,
                    "Bassik  Phago"      = my.bassik.phago.v2,
                    "Wyss Phago"         = my.wyss.phago.v2)
my.Venn <- Venn(my.criteria)

pdf(paste0(Sys.Date(),"_Female_Specific_Aging.pdf"))
plot(my.Venn, doWeights=F,type="ellipses")
dev.off()

#
my.criteria <- list("Female Aging Up"    = fem.spe.up,
                    "Female Aging Dwn"   = fem.spe.dwn,
                    "CRISPR  Phago"      = my.CIRSPR.phago)
my.Venn <- Venn(my.criteria)

pdf(paste0(Sys.Date(),"_Female_Specific_Aging_MERGED_PHAGO.pdf"))
plot(my.Venn, doWeights=F)
dev.off()

intersect(fem.spe.dwn,my.CIRSPR.phago)
# "Tm9sf3"  "Tbl1xr1" "Man1a2"  "B4galt1" "Dpy19l3" "Rps6ka3"

#####
my.criteria <- list("Female Aging Dwn"   = fem.spe.dwn,
                    "CRISPR  Phago"      = my.CIRSPR.phago)
my.Venn <- Venn(my.criteria)

pdf(paste0(Sys.Date(),"_Female_Specific_Aging_MERGED_PHAGO_DOWN_ONLY.pdf"))
plot(my.Venn, doWeights=T)
dev.off()

intersect(fem.spe.dwn,my.CIRSPR.phago)
# "Tm9sf3"  "Tbl1xr1" "Man1a2"  "B4galt1" "Dpy19l3" "Rps6ka3"







###########################
##### Male analysis
my.criteria <- list("Male Aging Up"      = mal.spe.up  ,
                    "Male Aging Dwn"     = mal.spe.dwn ,
                    "Bassik  Phago"      = my.bassik.phago.v2,
                    "Wyss Phago"         = my.wyss.phago.v2)
my.Venn <- Venn(my.criteria)

pdf(paste0(Sys.Date(),"_Male_Specific_Aging.pdf"))
plot(my.Venn, doWeights=F,type="ellipses")
dev.off()

#
my.criteria <- list("Male Aging Up"    = mal.spe.up  ,
                    "Male Aging Dwn"   = mal.spe.dwn ,
                    "CRISPR  Phago"    = my.CIRSPR.phago)
my.Venn <- Venn(my.criteria)

pdf(paste0(Sys.Date(),"_Male_Specific_Aging_MERGED_PHAGO.pdf"))
plot(my.Venn, doWeights=F)
dev.off()

intersect(mal.spe.dwn,my.CIRSPR.phago)
# # [1] "Ncor1"   "Scpep1"  "Snx29"   "Csnk2a1" "Ttll3"   "Calr"   




###########################
##### Female analysis (FDR10)
fem.mal.age.comp <- read.table('../../../../DEseq2_analysis/DEseq2/2022-03-30_Peritoneal_Macrophages_AGING_Separated_Sex_Merged_Table_ALL_genes_statistics.txt')
my.F_notM.fdr10 <- bitAnd(fem.mal.age.comp$padj_F < 0.1, fem.mal.age.comp$padj_M > 0.1) > 0

fem.spe.up.fdr10  <- rownames(fem.mal.age.comp)[bitAnd(fem.mal.age.comp$log2FoldChange_F > 0, my.F_notM.fdr10) > 0]
fem.spe.dwn.fdr10 <- rownames(fem.mal.age.comp)[bitAnd(fem.mal.age.comp$log2FoldChange_F < 0, my.F_notM.fdr10) > 0]

my.criteria <- list("Female Aging Up FDR10"    = fem.spe.up.fdr10,
                    "Female Aging Dwn FDR10"   = fem.spe.dwn.fdr10,
                    "CRISPR  Phago"             = my.CIRSPR.phago)
my.Venn <- Venn(my.criteria)

pdf(paste0(Sys.Date(),"_Female_Specific_Aging_MERGED_PHAGO_FDR10.pdf"))
plot(my.Venn, doWeights=F)
dev.off()

intersect(fem.spe.dwn.fdr10,my.CIRSPR.phago)
# [1] "Gna13"   "Hmgcs1"  "Tm9sf3"  "Tbl1xr1" "Man1a2"  "B4galt1" "Dpy19l3" "Magt1"   "Rps6ka3"
#not FDR5:
## "Gna13"   "Hmgcs1" "Magt1"


#####
my.criteria <- list("Female Aging Dwn FDR10"   = fem.spe.dwn.fdr10,
                    "CRISPR  Phago"           = my.CIRSPR.phago)
my.Venn <- Venn(my.criteria)

pdf(paste0(Sys.Date(),"_Female_Specific_Aging_MERGED_PHAGO_DOWN_ONLY_FDR10.pdf"))
plot(my.Venn, doWeights=T)
dev.off()

intersect(fem.spe.dwn.fdr10,my.CIRSPR.phago)
# [1] "Gna13"   "Hmgcs1"  "Tm9sf3"  "Tbl1xr1" "Man1a2"  "B4galt1" "Dpy19l3" "Magt1"   "Rps6ka3"
#not FDR5:
## "Gna13"   "Hmgcs1" "Magt1"







################################################################
# get boxplot of genes of interest
fem.mal.age.comp <- read.table('../../../DEseq2_analysis/DEseq2/2022-03-30_Peritoneal_Macrophages_AGING_Separated_Sex_Merged_Table_ALL_genes_statistics.txt')

#### DEBUG
# my.gene    <- "Rps6ka3"
# my.cts     <- tissue.cts
# my.sig.res <- fem.mal.age.comp
# my.ylim    <- c(8,14)
# y.pval     <- 13

get_sex_age_bplots <- function(my.gene, my.cts, my.sig.res, my.ylim = c(0,20), y.pval = 19) {
        # extract data frame of expression
        my.gene.exp <- data.frame("YF" = as.numeric(my.cts[my.gene,1:5]  ) ,
                                  "OF" = as.numeric(my.cts[my.gene,6:10]  ) ,
                                  "YM" = as.numeric(my.cts[my.gene,11:15] ) ,
                                  "OM" = as.numeric(my.cts[my.gene,16:20]) )
        
        pdfname <- paste0(Sys.Date(),"_",my.gene,"_Gene_Expression_boxplot_DEseq2.pdf")
        pdf(pdfname, height = 5, width = 3)
        boxplot(my.gene.exp,
                col = c("deeppink", "deeppink4","deepskyblue", "deepskyblue4"),
                ylab = "DESeq2 VST-normalized log2 counts",
                ylim = my.ylim, outline = F, main = my.gene)
        beeswarm::beeswarm(my.gene.exp, add = T, pch = 16, cex = 0.75)
        text(1.5, y.pval, signif(my.sig.res[my.gene,]$padj_F,2))
        text(3.5, y.pval, signif(my.sig.res[my.gene,]$padj_M,2))
        dev.off()
        
}

get_sex_age_bplots("Rps6ka3", tissue.cts, fem.mal.age.comp,  c(10,14), 13)
get_sex_age_bplots("Tbl1xr1", tissue.cts, fem.mal.age.comp,  c(10,14), 13)
get_sex_age_bplots("Tm9sf3" , tissue.cts, fem.mal.age.comp,  c(10,14), 13)
get_sex_age_bplots("Man1a2" , tissue.cts, fem.mal.age.comp,  c( 8,12), 11.5)
get_sex_age_bplots("Dpy19l3" , tissue.cts, fem.mal.age.comp,  c(6, 10), 9)
get_sex_age_bplots("B4galt1" , tissue.cts, fem.mal.age.comp,  c( 8,12), 11.5)


get_sex_age_bplots("Gna13"  , tissue.cts, fem.mal.age.comp,  c(10,14), 13)
get_sex_age_bplots("Hmgcs1"  , tissue.cts, fem.mal.age.comp, c(6, 12), 11)
get_sex_age_bplots("Magt1"  , tissue.cts, fem.mal.age.comp ,  c(8, 12), 11.5)

## "Gna13"   "Hmgcs1" "Magt1"
########################################################################################
get_sex_age_bplots("Esr1"  , tissue.cts, fem.mal.age.comp,  c(6,10), 10)
get_sex_age_bplots("Xist"  , tissue.cts, fem.mal.age.comp,  c(5,15), 15)

# male biased
get_sex_age_bplots("Rela"  , tissue.cts, fem.mal.age.comp,  c(8,12), 11.5)
get_sex_age_bplots("Hcar2"   , tissue.cts, fem.mal.age.comp,  c(6, 12), 12)


get_sex_age_bplots("Irf2"   , tissue.cts, fem.mal.age.comp,  c(10, 16), 15)


get_sex_age_bplots("Timd4"   , tissue.cts, fem.mal.age.comp,  c(10, 16), 15)
get_sex_age_bplots("Adgre1"   , tissue.cts, fem.mal.age.comp,  c(11, 17), 16)
get_sex_age_bplots("Cd74"   , tissue.cts, fem.mal.age.comp,  c(9, 15), 14)


get_sex_age_bplots("Clec7a"   , tissue.cts, fem.mal.age.comp,  c(8, 12), 11) # dectin 1
get_sex_age_bplots("Tlr2"     , tissue.cts, fem.mal.age.comp,  c(7, 12), 11) # 

tissue.cts["Tlr2",]


#######################
sink(file = paste(Sys.Date(),"_session_Info.txt", sep =""))
sessionInfo()
sink()




