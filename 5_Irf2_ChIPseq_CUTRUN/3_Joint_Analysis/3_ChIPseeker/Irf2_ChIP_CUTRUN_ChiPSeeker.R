setwd('/Users/berenice/Desktop/Irf2_CUTRUN/ChIPSeeker/ChIP_CUTRUN_summary')
options(stringsAsFactors = FALSE)

## loading packages
library(ChIPseeker)
library(clusterProfiler)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(ReactomePA)
library(ggplot2)
library(scales)

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

options(java.parameters = "-Xmx16g" )
require(openxlsx)

# 2025-12-02
# run Irf2 ChIP and CUTandRUN functional enrichment

######################################################################
### 1. read peak files
perimac.peak <- readPeakFile("/Users/berenice/Desktop/Irf2_CUTRUN/MSPC/MSPC_Irf2_PeriMac_4samples_ConsensusPeaks.simple.bed")
bmdm.peak    <- readPeakFile("/Users/berenice/Desktop/Irf2_CUTRUN/01_ChIPseq/MSPC/MSPC_Irf2_BMDM_ConsensusPeaks.simple.bed")

my.outprefix <- paste0(Sys.Date(),"_PeriMac_Irf2_BMDM_Perimac_CHIPSEEKER_analysis")
######################################################################

######################################################################
### 2. general QC plots

# Peak Annotation
perimac.peakAnno <- annotatePeak(perimac.peak , tssRegion=c(-5000, 5000), TxDb=txdb, annoDb="org.Mm.eg.db")
bmdm.peakAnno    <- annotatePeak(bmdm.peak    , tssRegion=c(-5000, 5000), TxDb=txdb, annoDb="org.Mm.eg.db")

peakAnnoList <- list("BMDM"    = bmdm.peakAnno,
                     "PeriMac" = perimac.peakAnno)

pdf(paste0(my.outprefix,"_plotAnnoBar.pdf"))
plotAnnoBar(peakAnnoList)
dev.off()

######################################################################
### 3. functional enrichment

#REACTOME
pm.react <- enrichPathway(as.data.frame(perimac.peakAnno)$geneId, organism = "mouse")
bm.react <- enrichPathway(as.data.frame(bmdm.peakAnno   )$geneId, organism = "mouse")

pdf(paste0(my.outprefix,"_REACTOME_enrichment.pdf"))
dotplot(pm.react, title = "PeriMac") + scale_y_discrete(labels = wrap_format(40))
dotplot(bm.react, title = "BMDM"   ) + scale_y_discrete(labels = wrap_format(40))
dev.off()

#GO
pm.go    <- enrichGO(as.data.frame(perimac.peakAnno)$geneId,  ont = "ALL", OrgDb="org.Mm.eg.db")
bm.go    <- enrichGO(as.data.frame(bmdm.peakAnno   )$geneId,  ont = "ALL", OrgDb="org.Mm.eg.db")

pdf(paste0(my.outprefix,"_GO_ALL_enrichment.pdf"))
dotplot(pm.go, title = "PeriMac" ) + scale_y_discrete(labels = wrap_format(40))
dotplot(bm.go, title = "BMDM"    ) + scale_y_discrete(labels = wrap_format(40))
dev.off()


# export results
chipseeker.res <- list("BMDM_REACTOME"    = data.frame(bm.react),
                       "PeriMac_REACTOME" = data.frame(pm.react),
                       "BMDM_GO_ALL"      = data.frame(bm.go),
                       "PeriMac_GO_ALL"   = data.frame(pm.go) )
write.xlsx(chipseeker.res, rowNames = F, file = paste0(my.outprefix,"_ClusterProfiler_REACTOME_GO_Irf2_BMDM_PeriMac_Results.xlsx"))
################################################################################################

#######################
sink(file = paste(my.outprefix,"_session_Info.txt", sep =""))
sessionInfo()
sink()

