# animal #1 for AC43 animals was < 95% purity QC
# exclude animal one from meta-analysis 

# Merge reads for meta peak calling
samtools merge POOLED_PeriMac_AC43_ATACseq.CLEAN.bam  ../BAM_clean/YM_PeriMac7_AC43_ATACseq.mp.srt.rmdup.bam ../BAM_clean/YM_PeriMac2_AC43_ATACseq.mp.srt.rmdup.bam ../BAM_clean/YM_PeriMac18_AC43_ATACseq.mp.srt.rmdup.bam ../BAM_clean/YM_PeriMac15_AC43_ATACseq.mp.srt.rmdup.bam ../BAM_clean/YM_PeriMac10_AC43_ATACseq.mp.srt.rmdup.bam ../BAM_clean/YF_PeriMac9_AC43_ATACseq.mp.srt.rmdup.bam ../BAM_clean/YF_PeriMac8_AC43_ATACseq.mp.srt.rmdup.bam ../BAM_clean/YF_PeriMac17_AC43_ATACseq.mp.srt.rmdup.bam ../BAM_clean/YF_PeriMac16_AC43_ATACseq.mp.srt.rmdup.bam ../BAM_clean/OM_PeriMac5_AC43_ATACseq.mp.srt.rmdup.bam ../BAM_clean/OM_PeriMac4_AC43_ATACseq.mp.srt.rmdup.bam ../BAM_clean/OM_PeriMac13_AC43_ATACseq.mp.srt.rmdup.bam ../BAM_clean/OM_PeriMac12_AC43_ATACseq.mp.srt.rmdup.bam ../BAM_clean/OF_PeriMac6_AC43_ATACseq.mp.srt.rmdup.bam ../BAM_clean/OF_PeriMac3_AC43_ATACseq.mp.srt.rmdup.bam ../BAM_clean/OF_PeriMac14_AC43_ATACseq.mp.srt.rmdup.bam ../BAM_clean/OF_PeriMac11_AC43_ATACseq.mp.srt.rmdup.bam

# call meta peaks, all age merged, all QC reads
macs2 callpeak -t POOLED_PeriMac_AC43_ATACseq.CLEAN.bam  -f "BAMPE" -g mm --keep-dup all --broad -n POOLED_PeriMac_AC43_ATACseq_MACS2  


MACS2==2.2.7.1