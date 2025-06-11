setwd('/Volumes/BB_HQ_1/Immune_sex_dimorphism_Aging/Macrophages/Clean_Macrophage_Data/ChIP-seq/ChIPseeker')
options(stringsAsFactors = FALSE)

# load libraries for plotting
## loading packages
library(ChIPseeker)
library(clusterProfiler)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(ReactomePA)

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

options(java.parameters = "-Xmx16g" )
require(openxlsx)

# 2025-01-16
# run Irf2 ChIP-seq functional enrichment

######################################################################
### 1. read peak file
peak <- readPeakFile("../MSPC/MSPC_Irf2_BMDM_ConsensusPeaks.noScore.bed")
peak


my.outprefix <- paste0(Sys.Date(),"_Irf2_BMDM_ChIPseq_CHIPSEEKER_analysis")

######################################################################
### 2. general QC plots
# ChIP peaks coverage plot
pdf(paste0(my.outprefix,"_covplot.pdf"))
covplot(peak)
dev.off()

# Profile of ChIP peaks binding to TSS regions
# First of all, for calculating the profile of ChIP peaks binding to TSS regions, we should prepare the TSS regions, which are defined as the flanking sequence of the TSS sites. Then align the peaks that are mapping to these regions, and generate the tagMatrix.
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrix <- getTagMatrix(peak, windows=promoter)

pdf(paste0(my.outprefix,"_tagMatrix.pdf"))
tagHeatmap(tagMatrix, xlim = c(-3000,3000))
dev.off()

pdf(paste0(my.outprefix,"_plotAvgProf.pdf"))
plotAvgProf(tagMatrix, xlim=c(-3000, 3000), conf = 0.95, resample = 1000)
dev.off()

# Peak Annotation
peakAnno <- annotatePeak(peak, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")

pdf(paste0(my.outprefix,"_plotAnnoPie.pdf"))
plotAnnoPie(peakAnno,   ndigit = 1)
dev.off()

pdf(paste0(my.outprefix,"_plotAnnoBar.pdf"))
plotAnnoBar(peakAnno)
dev.off()

######################################################################
### 3. functional enrichment

#REACTOME
pathway1 <- enrichPathway(as.data.frame(peakAnno)$geneId, organism = "mouse")
head(pathway1, 2)

pdf(paste0(my.outprefix,"_REACTOME_enrichment.pdf"))
dotplot(pathway1) + scale_y_discrete(labels = wrap_format(40))
dev.off()

#GO
pathway2 <- enrichGO(as.data.frame(peakAnno)$geneId, ont = "ALL", OrgDb="org.Mm.eg.db")
head(pathway2, 2)

pdf(paste0(my.outprefix,"_GO_ALL_enrichment.pdf"))
dotplot(pathway2) + scale_y_discrete(labels = wrap_format(40))
dev.off()


# export results
chipseeker.res <- list("REACTOME" = data.frame(pathway1),
                       "GO_ALL"   = data.frame(pathway2))
write.xlsx(chipseeker.res, rowNames = F, file = paste0(my.outprefix,"_ClusterProfiler_REACTOME_GO_Results.xlsx"))


################################################################################################

#######################
sink(file = paste(my.outprefix,"_session_Info.txt", sep =""))
sessionInfo()
sink()

