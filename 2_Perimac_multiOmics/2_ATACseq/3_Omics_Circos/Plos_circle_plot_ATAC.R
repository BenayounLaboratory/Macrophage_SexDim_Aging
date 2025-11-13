setwd('/Volumes/BB_Home_HQ/Immune_sex_dimorphism_Aging/Macrophages/Clean_Macrophage_Data/ATACseq//Omics_Circos')
options(stringsAsFactors = FALSE);

library("OmicCircos")
library('bitops')

# 2021-11-05
# process ATAC-seq dataset
# use Diffbind count on meta peaks
# run only high purity (> 95%)
# animal #1 for AC43 animals was < 95% purity QC


###########
## construct mm10 scaffold, letting an opening so the legend can be written in
my.mm10 <- read.table('mouse.mm10.genome', sep="\t")
my.mm10 <- data.frame('chr' = my.mm10$V1,
                      'start'= rep(0,length(my.mm10$V1)),
                      'end'= my.mm10$V2,
                      'bla1'= rep(NA,length(my.mm10$V1)),
                      'bla2'= rep(NA,length(my.mm10$V1))
)
my.mm10 <- my.mm10[1:21,]

mm10.db <- segAnglePo(my.mm10, seg=my.mm10$chr, angle.end = 350)

#########################################################################
my.ATAC.peaks.sex   <- read.csv('../DEseq2/2021-11-05_PeriMac_ATACseq_DESeq2_Analysis_SEX_DIM_all_genes_statistics_PeakAnnot.txt',sep="\t",header=T)

my.annotated.ATAC.sex <- data.frame('chr'            = my.ATAC.peaks.sex$Chr,
                                     'po'             = 0.5*(my.ATAC.peaks.sex$Start + my.ATAC.peaks.sex$End),
                                     'Gene.Name'      = my.ATAC.peaks.sex$PeakName,
                                     'log2FoldChange' = my.ATAC.peaks.sex$log2FoldChange,
                                     'padj'           = my.ATAC.peaks.sex$padj
)

my.sig.mph <- my.annotated.ATAC.sex[my.annotated.ATAC.sex$padj < 0.05,]
my.sig.mph$ABS_log2FoldChange <- abs(my.sig.mph$log2FoldChange)
my.sig.mph$FoldChange <- 2^my.sig.mph$log2FoldChange

# mapping	
# data frame or matrix containing mapping information and values. 
# Column 1: the segment or chromosome ID; 
# column 2: the position; 
# column 3: the position label (optional) or the value and additional columns are the values. 
# such as gene expression and copy number. Missing values are allowed and will be ignored.

#### divergently regulated genes
pdf(paste(Sys.Date(),"Circos_PeriMac_DE_Sex_FDR5_ATAC_height.pdf", sep = "_"), width = 7, height = 7)
par(mar=c(0, 0, 0, 0))
plot(c(1,800), c(1,800), type="n", axes=FALSE, xlab="", ylab="", main="")
circos(R=400, cir=mm10.db, W=5, type="chr", print.chr.lab=TRUE, scale=FALSE, col = "black")
circos(R=270, cir=mm10.db, W=125,
       mapping=my.sig.mph[my.sig.mph$log2FoldChange <0,],
       col.v=6,
       type="b",
       B=F, lwd=0.5,
       col=c("deepskyblue"));
circos(R=150, cir=mm10.db, W=125,
       mapping=my.sig.mph[my.sig.mph$log2FoldChange >0,],
       col.v=6,
       type="b",
       B=F, lwd=0.5,
       col=c("deeppink"));
dev.off()

#########################################################################
my.ATAC.peaks.sexAge <- read.csv('../DEseq2//2021-11-05_PeriMac_ATACseq_DESeq2_Analysis_AGING_Separated_Sex_Merged_Table_all_genes_statistics_PeakAnnot.txt',sep="\t",header=T)

#### divergently regulated genes
my.F_notM.tab <- bitAnd(my.ATAC.peaks.sexAge$padj_F < 0.05, my.ATAC.peaks.sexAge$padj_M > 0.1) > 0
my.M_notF.tab <- bitAnd(my.ATAC.peaks.sexAge$padj_M < 0.05, my.ATAC.peaks.sexAge$padj_F > 0.1) > 0

my.ATAC.peaks.sexAge$SexDimAge <- "No"
my.ATAC.peaks.sexAge$SexDimAge[my.F_notM.tab] <- "F_only"
my.ATAC.peaks.sexAge$SexDimAge[my.M_notF.tab] <- "M_only"

my.ATAC.peaks.sexAge.v2 <- my.ATAC.peaks.sexAge[my.ATAC.peaks.sexAge$SexDimAge != "No",] # 3869


my.annotated.ATAC.sexAge <- data.frame('chr'                  = my.ATAC.peaks.sexAge.v2$Chr,
                                        'po'                   = 0.5*(my.ATAC.peaks.sexAge.v2$Start + my.ATAC.peaks.sexAge.v2$End),
                                        'Gene.Name'            = my.ATAC.peaks.sexAge.v2$PeakName,
                                        'log2FoldChange_F'     = my.ATAC.peaks.sexAge.v2$log2FoldChange_F,
                                        'padj_F'               = my.ATAC.peaks.sexAge.v2$padj_F,
                                        'ABS_log2FoldChange_F' = abs(my.ATAC.peaks.sexAge.v2$log2FoldChange_F),
                                        
                                        'log2FoldChange_M'     = my.ATAC.peaks.sexAge.v2$log2FoldChange_M,
                                        'padj_M'               = my.ATAC.peaks.sexAge.v2$padj_M,
                                        'ABS_log2FoldChange_M' = abs(my.ATAC.peaks.sexAge.v2$log2FoldChange_M),
                                        'sexDim'               = my.ATAC.peaks.sexAge.v2$SexDimAge
                                        )

# mapping	
# data frame or matrix containing mapping information and values. 
# Column 1: the segment or chromosome ID; 
# column 2: the position; 
# column 3: the position label (optional) or the value and additional columns are the values. 
# such as gene expression and copy number. Missing values are allowed and will be ignored.

#### divergently regulated genes
pdf(paste(Sys.Date(),"Circos_PeriMax_Aging_DE_SexSpecific_ATAC_height_NoDir.pdf", sep = "_"), width = 7, height = 7)
par(mar=c(0, 0, 0, 0))
plot(c(1,800), c(1,800), type="n", axes=FALSE, xlab="", ylab="", main="")
circos(R=400, cir=mm10.db, W=5, type="chr", print.chr.lab=TRUE, scale=FALSE, col = "black")
circos(R=300, cir=mm10.db, W=100,
       mapping=my.annotated.ATAC.sexAge[my.annotated.ATAC.sexAge$sexDim == 'F_only',],
       col.v=6,
       type="b",
       B=F, lwd=0.5,
       col=c("deeppink"));
circos(R=180, cir=mm10.db, W=100,
       mapping=my.annotated.ATAC.sexAge[my.annotated.ATAC.sexAge$sexDim == 'M_only',],
       col.v=9,
       type="b",
       B=F, lwd=0.5,
       col=c("deepskyblue"));
dev.off()


pdf(paste(Sys.Date(),"Circos_PeriMax_Aging_DE_SexSpecific_ATAC_height_WithDir.pdf", sep = "_"), width = 7, height = 7)
par(mar=c(0, 0, 0, 0))
plot(c(1,800), c(1,800), type="n", axes=FALSE, xlab="", ylab="", main="")
circos(R=400, cir=mm10.db, W=5, type="chr", print.chr.lab=TRUE, scale=FALSE, col = "black")
circos(R=320, cir=mm10.db, W=80,
       mapping=my.annotated.ATAC.sexAge[bitAnd(my.annotated.ATAC.sexAge$sexDim == 'F_only', my.annotated.ATAC.sexAge$log2FoldChange_F >0)>0,],
       col.v=6,
       type="b",
       B=F, lwd=0.5,
       col=c("deeppink4"));
circos(R=240, cir=mm10.db, W=80,
       mapping=my.annotated.ATAC.sexAge[bitAnd(my.annotated.ATAC.sexAge$sexDim == 'F_only', my.annotated.ATAC.sexAge$log2FoldChange_F <0)>0,],
       col.v=6,
       type="b",
       B=F, lwd=0.5,
       col=c("deeppink"));


circos(R=140, cir=mm10.db, W=80,
       mapping=my.annotated.ATAC.sexAge[bitAnd(my.annotated.ATAC.sexAge$sexDim == 'M_only', my.annotated.ATAC.sexAge$log2FoldChange_M >0)>0,],
       col.v=6,
       type="b",
       B=F, lwd=0.5,
       col=c("deepskyblue4"));
circos(R=50, cir=mm10.db, W=80,
       mapping=my.annotated.ATAC.sexAge[bitAnd(my.annotated.ATAC.sexAge$sexDim == 'M_only', my.annotated.ATAC.sexAge$log2FoldChange_M <0)>0,],
       col.v=6,
       type="b",
       B=F, lwd=0.5,
       col=c("deepskyblue"));

dev.off()


#######################
sink(file = paste(Sys.Date(),"_OmicsCircos_ATAC_height_session_Info.txt", sep =""))
sessionInfo()
sink()
