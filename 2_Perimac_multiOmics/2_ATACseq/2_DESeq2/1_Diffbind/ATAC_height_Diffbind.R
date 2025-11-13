setwd('/Volumes/BB_Home_HQ/Immune_sex_dimorphism_Aging/Macrophages/Clean_Macrophage_Data/ATACseq/Diffbind')
options(stringsAsFactors = F)

library('DiffBind')
library("rtracklayer")

# 2021-10-14
# use diffbind to generate normlized count table
# run without low purity sample 1 (< 95%)


################################################################################
####################   Diffbind Analysis : ATAC height  ####################
################################################################################

perimac.aging <- dba(sampleSheet="PeriMac_ATAC_height.csv",skipLines=1,attributes=c(DBA_ID,DBA_CONDITION))
perimac.aging <- dba.count(perimac.aging)

pdf(paste0(Sys.Date(),"PeriMac_aging_ATAC_height_clust_heatmap_Stitched.pdf"))
plot(perimac.aging, colScheme="Reds")
dev.off()

perimac.aging <- dba.contrast(perimac.aging, categories=DBA_CONDITION,minMembers=2)
perimac.aging <- dba.analyze(perimac.aging,method=DBA_ALL_METHODS)
perimac.aging

# 6 Contrasts:
#   Group1 Members1 Group2 Members2 DB.edgeR DB.DESeq2
# 1     YF        4     OF        4     3882      1460
# 2     YF        4     YM        5     4939       684
# 3     YF        4     OM        4     8120      3613
# 4     OF        4     YM        5     6121      5251
# 5     OF        4     OM        4     3253       752
# 6     YM        5     OM        4     3273     30464

#############
# DESeq2 table norm needed for linear modeling
# export normalized counts
my.norm.counts <- dba.peakset(perimac.aging,bRetrieve=TRUE)
my.norm.counts <- data.frame(my.norm.counts)

write.table(my.norm.counts, 
            file = paste0(Sys.Date(),"_PeriMac_ATAC_Sex_Aging_Normalized_count_matrix.txt"),
            quote = F, col.names = T, row.names = F, sep = "\t")

write.table(my.norm.counts[,c(1:3)], 
            file = paste0(Sys.Date(),"_PeriMac_ATAC_Sex_Aging_diffbind_peaks.bed"),
            quote = F, col.names = F, row.names = F, sep = "\t")

################################################################################################

#######################
sink(file = paste(Sys.Date(),"Diffbind_ATAC_session_Info.txt", sep =""))
sessionInfo()
sink()

