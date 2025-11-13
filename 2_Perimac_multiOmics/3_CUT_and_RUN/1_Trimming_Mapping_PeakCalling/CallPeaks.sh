# remove PeriMac OF, animal 14 (90% purity only) and animal 16 (93% only)
# Keep only if > 95%

# Merge reads for meta peak calling
samtools merge POOLED_95PURITY_PeriMac_AC42_H3K4me3_250K.CLEAN.bam  ../BAM_clean/MHK202_PeriMac_AC42_21m_M_20_H3K4me3_250K.mp.srt.rmdup.bam ../BAM_clean/MHK201_PeriMac_AC42_5m_M_18_H3K4me3_250K.mp.srt.rmdup.bam ../BAM_clean/MHK200_PeriMac_AC42_5m_F_17_H3K4me3_250K.mp.srt.rmdup.bam ../BAM_clean/MHK198_PeriMac_AC42_5m_M_15_H3K4me3_250K.mp.srt.rmdup.bam ../BAM_clean/MHK196_PeriMac_AC42_21m_M_13_H3K4me3_250K.mp.srt.rmdup.bam ../BAM_clean/MHK195_PeriMac_AC42_21m_M_12_H3K4me3_250K.mp.srt.rmdup.bam ../BAM_clean/MHK194_PeriMac_AC42_21m_F_11_H3K4me3_250K.mp.srt.rmdup.bam ../BAM_clean/MHK193_PeriMac_AC42_5m_M_10_H3K4me3_250K.mp.srt.rmdup.bam ../BAM_clean/MHK192_PeriMac_AC42_5m_F_9_H3K4me3_250K.mp.srt.rmdup.bam ../BAM_clean/MHK191_PeriMac_AC42_5m_F_8_H3K4me3_250K.mp.srt.rmdup.bam ../BAM_clean/MHK190_PeriMac_AC42_5m_M_7_H3K4me3_250K.mp.srt.rmdup.bam ../BAM_clean/MHK189_PeriMac_AC42_21m_F_6_H3K4me3_250K.mp.srt.rmdup.bam ../BAM_clean/MHK188_PeriMac_AC42_21m_M_5_H3K4me3_250K.mp.srt.rmdup.bam ../BAM_clean/MHK187_PeriMac_AC42_21m_M_4_H3K4me3_250K.mp.srt.rmdup.bam ../BAM_clean/MHK186_PeriMac_AC42_21m_F_3_H3K4me3_250K.mp.srt.rmdup.bam ../BAM_clean/MHK185_PeriMac_AC42_5m_M_2_H3K4me3_250K.mp.srt.rmdup.bam ../BAM_clean/MHK184_PeriMac_AC42_5m_F_1_H3K4me3_250K.mp.srt.rmdup.bam 

# call meta peaks, all age merged, all QC reads
macs2 callpeak -t POOLED_95PURITY_PeriMac_AC42_H3K4me3_250K.CLEAN.bam  -f "BAMPE" -g mm --keep-dup all --broad -n POOLED_95PURITY_PeriMac_AC42_H3K4me3_250K_MACS2  

