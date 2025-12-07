setwd('/Users/berenice/Desktop/Irf2_CUTRUN/Omics_Circos')
options(stringsAsFactors = FALSE);
library("OmicCircos")
library('bitops')

# 2025-12-01
# plot Irf2 peaks


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

# Peak positions
irf2.cutrun        <- read.csv("../MSPC/MSPC_Irf2_PeriMac_4samples_ConsensusPeaks.simple.bed", header = F, sep = "\t")
irf2.chip          <- read.csv("../01_ChIPseq/MSPC/MSPC_Irf2_BMDM_ConsensusPeaks.simple.bed" , header = F, sep = "\t")

irf2.cutrun.track <- data.frame('chr'            = irf2.cutrun$V1,
                                'po'             = 0.5*(irf2.cutrun$V2 + irf2.cutrun$V3),
                                'val'            = 1)

irf2.chip.track   <- data.frame('chr'            = irf2.chip$V1,
                                'po'             = 0.5*(irf2.chip$V2 + irf2.chip$V3),
                                'val'            = 1)

# Diff peaks
irf2.cutrun.fdr5   <- read.csv("../DiffBind/MSPC_4samples/2025-11-28_PeriMac_Irf2_stringent_Aging_DiffBind_DESeq2_FDR5_Peaks.txt", header = T, sep = "\t")

irf2.cutrun.track.fdr <- data.frame('chr'            = irf2.cutrun.fdr5$Chr,
                                    'po'             = 0.5*(irf2.cutrun.fdr5$Start + irf2.cutrun.fdr5$End),
                                    'val'            = irf2.cutrun.fdr5$Fold)

# mapping	
# data frame or matrix containing mapping information and values. 
# Column 1: the segment or chromosome ID; 
# column 2: the position; 
# column 3: the position label (optional) or the value and additional columns are the values. 
# such as gene expression and copy number. Missing values are allowed and will be ignored.

#### 
pdf(paste(Sys.Date(),"Circos_BMDM_PeriMac_Irf2.pdf", sep = "_"), width = 4, height = 4)
par(mar=c(1, 1, 1, 1))
plot(c(1,800), c(1,800), type="n", axes=FALSE, xlab="", ylab="", main="")
circos(R=400, cir=mm10.db, W=5, type="chr", print.chr.lab=TRUE, scale=FALSE, col = "black")
circos(R=310, cir=mm10.db, W=80,
       mapping=irf2.chip.track,
       type="b3",
       B=F, lwd=0.5,
       col=c("pink3"));
circos(R=185, cir=mm10.db, W=80,
       mapping=irf2.cutrun.track,
       type="b3",
       B=F, lwd=0.5,
       col=c("deeppink"));
circos(R=120, cir=mm10.db, W=60,
       mapping=irf2.cutrun.track.fdr,
       col.v = 3,
       type="b",
       B=F, lwd=0.5,
       col=c("grey"));
dev.off()

#######################
sink(file = paste(Sys.Date(),"_OmicsCircos_Irf2_session_Info.txt", sep =""))
sessionInfo()
sink()
