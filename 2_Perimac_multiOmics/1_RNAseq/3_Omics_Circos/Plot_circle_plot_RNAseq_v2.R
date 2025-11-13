setwd('/Volumes/BB_Home_HQ/Immune_sex_dimorphism_Aging/Macrophages/Clean_Macrophage_Data/RNAseq/OmicsCircos')
options(stringsAsFactors = FALSE)

library("OmicCircos")
library('bitops')


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

### 
my.homer.annot <- read.csv('2019-06-12_HOMER_annotated_TSS_mm10.xls',sep="\t",header=T)

######################################################################
my.RNA.sexAge <- read.csv('../DEseq2_analysis/DESeq2//2022-03-30_Peritoneal_Macrophages_AGING_Separated_Sex_Merged_Table_ALL_genes_statistics.txt',sep="\t",header=T)

# add gene name and merge with HOMER info
my.RNA.sexAge$Gene.Name <- rownames(my.RNA.sexAge)
my.RNA.sexAge.annotated <- merge(my.homer.annot[,c(2:4,16)], my.RNA.sexAge, by= "Gene.Name")


#### divergently regulated genes
my.F_notM.tab <- bitAnd(my.RNA.sexAge.annotated$padj_F < 0.05, my.RNA.sexAge.annotated$padj_M > 0.1) > 0
my.M_notF.tab <- bitAnd(my.RNA.sexAge.annotated$padj_M < 0.05, my.RNA.sexAge.annotated$padj_F > 0.1) > 0

my.RNA.sexAge.annotated$SexDimAge <- "No"
my.RNA.sexAge.annotated$SexDimAge[my.F_notM.tab] <- "F_only"
my.RNA.sexAge.annotated$SexDimAge[my.M_notF.tab] <- "M_only"

my.RNA.sexAge.v2 <- my.RNA.sexAge.annotated[my.RNA.sexAge.annotated$SexDimAge != "No",] # 864


my.annotated.RNA.sexAge <- data.frame('chr'                  = my.RNA.sexAge.v2$Chr,
                                      'po'                   = 0.5*(my.RNA.sexAge.v2$Start + my.RNA.sexAge.v2$End),
                                      'Gene.Name'            = my.RNA.sexAge.v2$Gene.Name,
                                      'log2FoldChange_F'     = my.RNA.sexAge.v2$log2FoldChange_F,
                                      'padj_F'               = my.RNA.sexAge.v2$padj_F,
                                      'ABS_log2FoldChange_F' = abs(my.RNA.sexAge.v2$log2FoldChange_F),
                                      'log2FoldChange_M'     = my.RNA.sexAge.v2$log2FoldChange_M,
                                      'padj_M'               = my.RNA.sexAge.v2$padj_M,
                                      'ABS_log2FoldChange_M' = abs(my.RNA.sexAge.v2$log2FoldChange_M),
                                      'sexDim'               = my.RNA.sexAge.v2$SexDimAge
)

# mapping	
# data frame or matrix containing mapping information and values. 
# Column 1: the segment or chromosome ID; 
# column 2: the position; 
# column 3: the position label (optional) or the value and additional columns are the values. 
# such as gene expression and copy number. Missing values are allowed and will be ignored.

#### divergently regulated genes
pdf(paste(Sys.Date(),"Circos_PeriMac_Aging_DE_SexSpecific_RNA_NoDir.pdf", sep = "_"), width = 7, height = 7)
par(mar=c(0, 0, 0, 0))
plot(c(1,800), c(1,800), type="n", axes=FALSE, xlab="", ylab="", main="")
circos(R=400, cir=mm10.db, W=5, type="chr", print.chr.lab=TRUE, scale=FALSE, col = "black")
circos(R=300, cir=mm10.db, W=100,
       mapping=my.annotated.RNA.sexAge[my.annotated.RNA.sexAge$sexDim == 'F_only',],
       col.v=6,
       type="b",
       B=F, lwd=0.5,
       col=c("deeppink"));
circos(R=180, cir=mm10.db, W=100,
       mapping=my.annotated.RNA.sexAge[my.annotated.RNA.sexAge$sexDim == 'M_only',],
       col.v=9,
       type="b",
       B=F, lwd=0.5,
       col=c("deepskyblue"));
dev.off()


pdf(paste(Sys.Date(),"Circos_PeriMac_Aging_DE_SexSpecific_RNA_WithDir.pdf", sep = "_"), width = 7, height = 7)
par(mar=c(0, 0, 0, 0))
plot(c(1,800), c(1,800), type="n", axes=FALSE, xlab="", ylab="", main="")
circos(R=400, cir=mm10.db, W=5, type="chr", print.chr.lab=TRUE, scale=FALSE, col = "black")
circos(R=320, cir=mm10.db, W=80,
       mapping=my.annotated.RNA.sexAge[bitAnd(my.annotated.RNA.sexAge$sexDim == 'F_only', my.annotated.RNA.sexAge$log2FoldChange_F >0)>0,],
       col.v=6,
       type="b",
       B=F, lwd=0.5,
       col=c("deeppink4"));
circos(R=240, cir=mm10.db, W=80,
       mapping=my.annotated.RNA.sexAge[bitAnd(my.annotated.RNA.sexAge$sexDim == 'F_only', my.annotated.RNA.sexAge$log2FoldChange_F <0)>0,],
       col.v=6,
       type="b",
       B=F, lwd=0.5,
       col=c("deeppink"));


circos(R=140, cir=mm10.db, W=80,
       mapping=my.annotated.RNA.sexAge[bitAnd(my.annotated.RNA.sexAge$sexDim == 'M_only', my.annotated.RNA.sexAge$log2FoldChange_M >0)>0,],
       col.v=6,
       type="b",
       B=F, lwd=0.5,
       col=c("deepskyblue4"));
circos(R=50, cir=mm10.db, W=80,
       mapping=my.annotated.RNA.sexAge[bitAnd(my.annotated.RNA.sexAge$sexDim == 'M_only', my.annotated.RNA.sexAge$log2FoldChange_M <0)>0,],
       col.v=6,
       type="b",
       B=F, lwd=0.5,
       col=c("deepskyblue"));

dev.off()

#######################
sink(file = paste(Sys.Date(),"_OmicsCircos_RNA_session_Info.txt", sep =""))
sessionInfo()
sink()
