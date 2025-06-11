# cat ../BAM_clean/01_0HW5_02F7USC_BMDM-YF-1_IRF2_mm_i14-2.FIXSEQ_CLEANED_READS.bed ../BAM_clean/02_0HW6_02F7USC_BMDM-YF-2_IRF2_mm_i16-2.FIXSEQ_CLEANED_READS.bed > POOLED_BMDM_YF_IRF2_ChIP.FIXSEQ_CLEANED_READS.bed


macs2 callpeak -t POOLED_BMDM_YF_IRF2_ChIP.FIXSEQ_CLEANED_READS.bed -c ../BAM_clean/00_0I13_02F7USC_Pooled_Input_mm_i29-4.FIXSEQ_CLEANED_READS.bed  -f "BED" -g mm --keep-dup all -n BMDM_Pooled_YF_IRF2_ChIPseq_MACS2  
macs2 callpeak -t ../BAM_clean/01_0HW5_02F7USC_BMDM-YF-1_IRF2_mm_i14-2.FIXSEQ_CLEANED_READS.bed -c ../BAM_clean/00_0I13_02F7USC_Pooled_Input_mm_i29-4.FIXSEQ_CLEANED_READS.bed  -f "BED" -g mm --keep-dup all -n BMDM_YF1_IRF2_ChIPseq_MACS2  
macs2 callpeak -t ../BAM_clean/02_0HW6_02F7USC_BMDM-YF-2_IRF2_mm_i16-2.FIXSEQ_CLEANED_READS.bed -c ../BAM_clean/00_0I13_02F7USC_Pooled_Input_mm_i29-4.FIXSEQ_CLEANED_READS.bed  -f "BED" -g mm --keep-dup all -n BMDM_YF2_IRF2_ChIPseq_MACS2  
