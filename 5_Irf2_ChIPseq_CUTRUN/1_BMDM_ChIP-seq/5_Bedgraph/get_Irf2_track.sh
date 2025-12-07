# samtools merge IRF2_YF_BMDM_merged.bam 01_0HW5_02F7USC_BMDM-YF-1_IRF2_mm_i14-2.mp.srt.bam 02_0HW6_02F7USC_BMDM-YF-2_IRF2_mm_i16-2.mp.srt.bam
# 
# samtools sort IRF2_YF_BMDM_merged.bam IRF2_YF_BMDM_merged.srt 
# samtools index IRF2_YF_BMDM_merged.srt.bam

#########

# pair with correct input
makeUCSCfile 01_0HW5_02F7USC_BMDM-YF-1_IRF2_mm_i14-2_FIXSEQ_TAGs -i 00_0I13_02F7USC_Pooled_Input_mm_i29-4_FIXSEQ_TAGs  -o BMDM_YF1_IRF2_ChIPseq.bedgraph  -style chipseq
makeUCSCfile 02_0HW6_02F7USC_BMDM-YF-2_IRF2_mm_i16-2_FIXSEQ_TAGs -i 00_0I13_02F7USC_Pooled_Input_mm_i29-4_FIXSEQ_TAGs  -o BMDM_YF2_IRF2_ChIPseq.bedgraph  -style chipseq


makeTagDirectory IRF2_YF_BMDM_merged_TAGs IRF2_YF_BMDM_merged.srt.bam -genome mm10
makeUCSCfile IRF2_YF_BMDM_merged_TAGs -i 00_0I13_02F7USC_Pooled_Input_mm_i29-4_FIXSEQ_TAGs  -o aNSC_MSI1_RIP_Sample1_1_RIPseq.bedgraph  -style chipseq
