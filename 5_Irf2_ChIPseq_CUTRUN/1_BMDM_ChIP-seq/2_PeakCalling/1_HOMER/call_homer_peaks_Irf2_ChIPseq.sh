findPeaks /Volumes/BB_Travel_2/Mph_Aging/ChIP-seq/BAM_clean/01_0HW5_02F7USC_BMDM-YF-1_IRF2_mm_i14-2_FIXSEQ_TAGs  -o 01_BMDM_IRF2_ChIPseq_YF_1.pos    -style factor -i /Volumes/BB_Travel_2/Mph_Aging/ChIP-seq/BAM_clean/00_0I13_02F7USC_Pooled_Input_mm_i29-4_FIXSEQ_TAGs
findPeaks /Volumes/BB_Travel_2/Mph_Aging/ChIP-seq/BAM_clean/02_0HW6_02F7USC_BMDM-YF-2_IRF2_mm_i16-2_FIXSEQ_TAGs  -o 02_BMDM_IRF2_ChIPseq_YF_2.pos    -style factor -i /Volumes/BB_Travel_2/Mph_Aging/ChIP-seq/BAM_clean/00_0I13_02F7USC_Pooled_Input_mm_i29-4_FIXSEQ_TAGs

parse_pos_for_MSPC_input.pl 01_BMDM_IRF2_ChIPseq_YF_1.pos
parse_pos_for_MSPC_input.pl 02_BMDM_IRF2_ChIPseq_YF_2.pos


#      541 01_BMDM_IRF2_ChIPseq_YF_1.cl.bed
#      949 02_BMDM_IRF2_ChIPseq_YF_2.cl.bed
