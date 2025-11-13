setwd('/Volumes/BB_HQ_1/Immune_sex_dimorphism_Aging/Macrophages/Clean_Macrophage_Data/CUT_and_RUN/H3K4me3_Aging/H3K4me3_height/Diffbind/')
options(stringsAsFactors = F)

library('DiffBind')
library("rtracklayer")

# 2021-09-21
# use diffbind to generate normlized count table

# 2021-09-30
# rerun but without low purity sample 14 and 16 (< 95%)

# 2024-03-28
# Run BMDM data

######################################################################################
####################   Diffbind Analysis : H3K4me3 height PeriMac ####################
######################################################################################

perimac.aging <- dba(sampleSheet="PeriMac_K4me3_height_v2.csv",skipLines=1,attributes=c(DBA_ID,DBA_CONDITION))
perimac.aging <- dba.count(perimac.aging)

pdf(paste0(Sys.Date(),"PeriMac_aging_H3K4me3_height_clust_heatmap_Stitched.pdf"))
plot(perimac.aging, colScheme="Reds")
dev.off()

perimac.aging <- dba.contrast(perimac.aging, categories=DBA_CONDITION,minMembers=2)
perimac.aging <- dba.analyze(perimac.aging,method=DBA_ALL_METHODS)
perimac.aging

# 6 Contrasts:
#   Group1 Members1 Group2 Members2 DB.edgeR DB.DESeq2
# 1     YF        4     YM        5     2686      2300
# 2     YF        4     OF        3     3150      2043
# 3     YF        4     OM        5     4195      2967
# 4     YM        5     OF        3     1480       955
# 5     YM        5     OM        5      631       469
# 6     OF        3     OM        5      353       136


#############
# DESeq2 table norm needed for linear modeling
# export normalized counts
my.norm.counts <- dba.peakset(perimac.aging,bRetrieve=TRUE)
my.norm.counts <- data.frame(my.norm.counts)

write.table(my.norm.counts, 
            file = paste0(Sys.Date(),"_PeriMac_H3K4me3_Sex_Aging_Normalized_count_matrix.txt"),
            quote = F, col.names = T, row.names = F, sep = "\t")

write.table(my.norm.counts[,c(1:3)], 
            file = paste0(Sys.Date(),"_PeriMac_H3K4me3_Sex_Aging_diffbind_peaks.bed"),
            quote = F, col.names = F, row.names = F, sep = "\t")
################################################################################################

#######################
sink(file = paste(Sys.Date(),"Diffbind_H3K4me3_CUTandRUN_session_Info.txt", sep =""))
sessionInfo()
sink()

