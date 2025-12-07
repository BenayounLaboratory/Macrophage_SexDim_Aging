setwd('/Users/berenice/Desktop/Irf2_CUTRUN/DiffBind/MSPC_4samples')
options(stringsAsFactors = F)

library('DiffBind')
library("rtracklayer")
library("pheatmap")
library("beeswarm")

# 2025-11-28
# use diffbind to analyze Irf cut and run
# use stringent peak definition to improve statistical power (peak in all samples)

######################################################################################
####################   Diffbind Analysis : Irf2 binding PeriMac ####################
######################################################################################

################################################################################
##### A. Run Duffbind/DEseq2 analysis using dm6 SpikeIn data

# read design
perimac.aging <- dba(sampleSheet="PeriMac_Irf2_4s.csv",skipLines=1,attributes=c(DBA_ID,DBA_CONDITION))

# remove ENCODE blacklist intervals
perimac.aging <- dba.blacklist(perimac.aging, blacklist = DBA_BLACKLIST_MM10) # Removed: 72 merged (of 1565) and 72 (of 1565) consensus.

# count reads
perimac.aging <- dba.count(perimac.aging)

# normalize to spike in data (bams specified in samplesheet)
perimac.aging <- dba.normalize(perimac.aging, spikein = T)

# set comparisons for differential analysis
perimac.aging <- dba.contrast(perimac.aging, categories=DBA_CONDITION, minMembers=2)

# analyze for differential peaks
perimac.aging <- dba.analyze(perimac.aging, method=DBA_DESEQ2)

perimac.aging
# 4 Samples, 1493 sites in matrix:
#             ID Condition Replicate    Reads FRiP
# 1 PeriMac_YF_1        YF         1 14175150 0.01
# 2 PeriMac_YF_2        YF         2 15146434 0.01
# 3 PeriMac_OF_1        OF         1 11691807 0.01
# 4 PeriMac_OF_2        OF         2 11236872 0.01

# Design: [~Condition] | 1 Contrast:
#   Factor Group Samples Group2 Samples2 DB.DESeq2
# 1 Condition    OF       2     YF        2       541

pdf(paste0(Sys.Date(),"PeriMac_aging_Irf2_clust_heatmap_Stitched_stringent.pdf"))
plot(perimac.aging, colScheme="Reds")
dev.off()

save(perimac.aging, file = paste0(Sys.Date(),"_DiffBind_Irf2_CUTRUN_stringent_object.RData"))

# PCA analysis
pdf(paste0(Sys.Date(),"PeriMac_aging_Irf2_PCA_stringent.pdf"))
dba.plotPCA(perimac.aging, vColors = c("deeppink","deeppink4"))
dev.off()

# Signal in peaks
pdf(paste0(Sys.Date(),"PeriMac_aging_Irf2_boxplot_signal_stringent.pdf"))
dba.plotBox(perimac.aging, vColors = c("deeppink4","deeppink"), bDBIncreased = F, bDBDecreased = F )
dev.off()

pvals.1 <- dba.plotBox(perimac.aging, vColors = c("deeppink4","deeppink"), bDBIncreased = F, bDBDecreased = F )
pvals.1
#          OF.DB    YF.DB
# OF.DB 1.00e+00 2.56e-90
# YF.DB 2.56e-90 1.00e+00

###### Plot volcanos
library(ggplot2) 
theme_set(theme_bw())

pdf(paste0(Sys.Date(),"PeriMac_aging_Irf2_volcanoplot_stringent.pdf"))
dba.plotVolcano(perimac.aging)
dev.off()


###### Plot heatmaps
hmap <- colorRampPalette(c("dodgerblue3", "gray95", "firebrick3"))(100)

pdf(paste0(Sys.Date(),"PeriMac_aging_Irf2_HEATMAP_FDR5.pdf"))
dba.plotHeatmap(perimac.aging, correlations=FALSE, scale="row", colScheme = hmap, ColAttributes = NULL, th = 0.05)
dev.off()
################################################################################

################################################################################
##### B. Extract DEseq2 results and normalized counts

## all peaks results
perimac.aging.DB <- dba.report(perimac.aging, th = 1)
write.table(data.frame(perimac.aging.DB), file=paste0(Sys.Date(),"PeriMac_aging_Irf2_stringent_DESeq2_res_allPeaks.txt")    , sep = "\t", row.names = F, quote = F)

## FDR5 only
perimac.aging.DB.fdr5 <- dba.report(perimac.aging, th = 0.05)
write.table(data.frame(perimac.aging.DB.fdr5), file=paste0(Sys.Date(),"PeriMac_aging_Irf2_stringent_DESeq2_res_FDR5.txt")    , sep = "\t", row.names = F, quote = F)

### Extract normalized counts
my.norm.counts <- dba.peakset(perimac.aging, bRetrieve=TRUE)
my.norm.counts <- data.frame(my.norm.counts)

write.table(my.norm.counts, 
            file = paste0(Sys.Date(),"_PeriMac_Irf2_stringent_Aging_Normalized_count_matrix.txt"),
            quote = F, col.names = T, row.names = F, sep = "\t")

write.table(my.norm.counts[,c(1:3)], 
            file = paste0(Sys.Date(),"_PeriMac_Irf2_stringent_Aging_diffbind_peaks.bed"),
            quote = F, col.names = F, row.names = F, sep = "\t")
################################################################################

################################################################################
##### C. Plot MDS and heatmaps for figures

##### MDS analysis
mds.result <- cmdscale(1-cor(my.norm.counts[,6:9],method="spearman"), k = 2, eig = FALSE, add = FALSE, x.ret = FALSE)
x <- mds.result[, 1]
y <- mds.result[, 2]

pdf(paste0(Sys.Date(),"PeriMac_aging_Irf2_stringent_MDS_plot.pdf"))
plot(x, y, 
     xlab = "MDS dimension 1", ylab = "MDS dimension 2",
     main="Multi-dimensional Scaling", 
     cex=3, col= c(rep("deeppink",2),rep("deeppink4",2)), pch= 16,
     cex.lab = 1.5,
     cex.axis = 1.5,
     xlim = c(-0.15,0.15),
     ylim = c(-0.08,0.08),
     las = 1
)
dev.off()


#### plots heatmaps

# significant peaks
sig.peaks          <- data.frame(perimac.aging.DB.fdr5) ### 541
sig.peaks$PeakName <- paste(sig.peaks$seqnames, sig.peaks$start,sig.peaks$end, sep = "_")

# all peaks
all.peaks          <- data.frame(perimac.aging.DB) ### 1493
all.peaks$PeakName <- paste(all.peaks$seqnames, all.peaks$start,all.peaks$end, sep = "_")

# add name to count table
my.norm.counts.annot          <- my.norm.counts
my.norm.counts.annot$PeakName <- paste(my.norm.counts$seqnames, my.norm.counts$start, my.norm.counts$end, sep = "_")

# get significant peaks data for heatmap
heat.count <- my.norm.counts.annot[my.norm.counts.annot$PeakName %in% sig.peaks$PeakName,6:9] # 423

pdf(paste0(Sys.Date(),"PeriMac_aging_Irf2_FDR5_peaks_Heatmaps.pdf"), onefile = F)
my.heatmap.title <- paste("Aging significant (FDR<5%), ", nrow(heat.count), " peaks",sep="")
pheatmap(heat.count,
         cluster_cols = F,
         cluster_rows = T,
         colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
         show_rownames = F, scale="row",
         main = my.heatmap.title, cellwidth = 15, border_color = NA, cellheight = 0.2)
dev.off()

################################################################################

################################################################################
##### D. Load HOMER annotation

my.peak.annot <- read.csv('HOMER_2025-11-28_PeriMac_Irf2_stringent_Aging_diffbind_peaks.xls', sep = "\t", header = T)
colnames(my.peak.annot)[1] <- "PeakID"
my.peak.annot$PeakName <- paste(my.peak.annot$Chr,my.peak.annot$Start-1,my.peak.annot$End,sep = "_")

## merge annotation
data.frame(perimac.aging.DB)

res.irf2.annot <- merge(my.peak.annot[,c("PeakName","Chr", "Start", "End", "Annotation", "Distance.to.TSS", "Gene.Name","Gene.Description" )],all.peaks, by = "PeakName") # 2121
res.irf2.fdr5  <- merge(my.peak.annot[,c("PeakName","Chr", "Start", "End", "Annotation", "Distance.to.TSS", "Gene.Name","Gene.Description" )],sig.peaks, by = "PeakName") # 423

###
write.table(all.peaks[,c(1:3)], file = paste0(Sys.Date(),"_PeriMac_Irf2_stringent_Aging_diffbind_peaks_BACKGROUND.bed"),quote = F, col.names = F, row.names = F, sep = "\t")
write.table(sig.peaks[,c(1:3)], file = paste0(Sys.Date(),"_PeriMac_Irf2_stringent_Aging_diffbind_peaks_FDR5.bed"),quote = F, col.names = F, row.names = F, sep = "\t")


write.table(res.irf2.annot, 
            file = paste0(Sys.Date(),"_PeriMac_Irf2_stringent_Aging_DiffBind_DESeq2_ALL_Peaks.txt"),
            quote = F, col.names = T, row.names = F, sep = "\t")

write.table(res.irf2.fdr5, 
            file = paste0(Sys.Date(),"_PeriMac_Irf2_stringent_Aging_DiffBind_DESeq2_FDR5_Peaks.txt"),
            quote = F, col.names = T, row.names = F, sep = "\t")

################################################################################

################################################################################
##### E. Zoom on Hk3 promoter

res.irf2.annot[res.irf2.annot$Gene.Name %in% "Hk3",]
# PeakName   Chr    Start      End Annotation Distance.to.TSS Gene.Name
# 403 chr13_55022528_55022928 chr13 55022529 55022928 Intergenic           -1343       Hk3
# Gene.Description seqnames    start      end width strand     Conc  Conc_OF  Conc_YF      Fold
# 403     hexokinase 3    chr13 55022528 55022928   401      * 5.948694 5.097376 6.480484 -1.081851
# p.value        FDR
# 403 0.01137693 0.04528917

hk3.count <- my.norm.counts.annot[my.norm.counts.annot$PeakName %in% "chr13_55022528_55022928",6:9]

pdf(paste0(Sys.Date(),"PeriMac_aging_Irf2_Hk3_peak_Beeswarm.pdf"), width = 3.5, height = 5)
beeswarm(list("YF" = hk3.count[1:2], "OF" = hk3.count[3:4]),
         pch = 16, las = 1, pwcol = c("deeppink","deeppink","deeppink4","deeppink4"),
         main = "chr13_55022528_55022928 (Hk3 promoter)", cex = 2,
         ylim = c(0,150), ylab = "SpikeIn normalized peak signal (A.U.)")
text(1.5,145, "4.5E-2")
dev.off()


################################################################################

#######################
sink(file = paste(Sys.Date(),"Diffbind_Irf2_CUTandRUN_session_Info.txt", sep =""))
sessionInfo()
sink()

