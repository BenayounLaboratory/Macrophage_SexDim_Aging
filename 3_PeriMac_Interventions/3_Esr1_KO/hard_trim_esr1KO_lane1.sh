# ### use trimgalore to recapitulate the behavior of fastx_trimmer (doesn't run on newer macos)
# ## cut to 100 to help mapping with STAR

# trim_galore --paired --hardtrim5 100 ../FASTQ/lane1/KO_1_2_CKDL240029243-1A_22FWHJLT4_L1_1.fq.gz   ../FASTQ/lane1/KO_1_2_CKDL240029243-1A_22FWHJLT4_L1_2.fq.gz
# trim_galore --paired --hardtrim5 100 ../FASTQ/lane1/KO_2_4_CKDL240029243-1A_22FWHJLT4_L1_1.fq.gz   ../FASTQ/lane1/KO_2_4_CKDL240029243-1A_22FWHJLT4_L1_2.fq.gz
# trim_galore --paired --hardtrim5 100 ../FASTQ/lane1/KO_3_6_CKDL240029243-1A_22FWHJLT4_L1_1.fq.gz   ../FASTQ/lane1/KO_3_6_CKDL240029243-1A_22FWHJLT4_L1_2.fq.gz
# trim_galore --paired --hardtrim5 100 ../FASTQ/lane1/KO_4_8_CKDL240029243-1A_22FWHJLT4_L1_1.fq.gz   ../FASTQ/lane1/KO_4_8_CKDL240029243-1A_22FWHJLT4_L1_2.fq.gz
# trim_galore --paired --hardtrim5 100 ../FASTQ/lane1/KO_5_10_CKDL240029243-1A_22FWHJLT4_L1_1.fq.gz  ../FASTQ/lane1/KO_5_10_CKDL240029243-1A_22FWHJLT4_L1_2.fq.gz
# trim_galore --paired --hardtrim5 100 ../FASTQ/lane1/KO_6_12_CKDL240029243-1A_22FWHJLT4_L1_1.fq.gz  ../FASTQ/lane1/KO_6_12_CKDL240029243-1A_22FWHJLT4_L1_2.fq.gz
# trim_galore --paired --hardtrim5 100 ../FASTQ/lane1/WT_1_1_CKDL240029243-1A_22FWHJLT4_L1_1.fq.gz   ../FASTQ/lane1/WT_1_1_CKDL240029243-1A_22FWHJLT4_L1_2.fq.gz
# trim_galore --paired --hardtrim5 100 ../FASTQ/lane1/WT_2_3_CKDL240029243-1A_22FWHJLT4_L1_1.fq.gz   ../FASTQ/lane1/WT_2_3_CKDL240029243-1A_22FWHJLT4_L1_2.fq.gz
# trim_galore --paired --hardtrim5 100 ../FASTQ/lane1/WT_3_5_CKDL240029243-1A_22FWHJLT4_L1_1.fq.gz   ../FASTQ/lane1/WT_3_5_CKDL240029243-1A_22FWHJLT4_L1_2.fq.gz
# trim_galore --paired --hardtrim5 100 ../FASTQ/lane1/WT_4_7_CKDL240029243-1A_22FWHJLT4_L1_1.fq.gz   ../FASTQ/lane1/WT_4_7_CKDL240029243-1A_22FWHJLT4_L1_2.fq.gz
# trim_galore --paired --hardtrim5 100 ../FASTQ/lane1/WT_5_9_CKDL240029243-1A_22FWHJLT4_L1_1.fq.gz   ../FASTQ/lane1/WT_5_9_CKDL240029243-1A_22FWHJLT4_L1_2.fq.gz
# trim_galore --paired --hardtrim5 100 ../FASTQ/lane1/WT_6_11_CKDL240029243-1A_22FWHJLT4_L1_1.fq.gz  ../FASTQ/lane1/WT_6_11_CKDL240029243-1A_22FWHJLT4_L1_2.fq.gz
# 
# # ### cut 9 based on fastqc
# trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired KO_1_2_CKDL240029243-1A_22FWHJLT4_L1_1.100bp_5prime.fq.gz   KO_1_2_CKDL240029243-1A_22FWHJLT4_L1_2.100bp_5prime.fq.gz
# trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired KO_2_4_CKDL240029243-1A_22FWHJLT4_L1_1.100bp_5prime.fq.gz   KO_2_4_CKDL240029243-1A_22FWHJLT4_L1_2.100bp_5prime.fq.gz
# trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired KO_3_6_CKDL240029243-1A_22FWHJLT4_L1_1.100bp_5prime.fq.gz   KO_3_6_CKDL240029243-1A_22FWHJLT4_L1_2.100bp_5prime.fq.gz
# trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired KO_4_8_CKDL240029243-1A_22FWHJLT4_L1_1.100bp_5prime.fq.gz   KO_4_8_CKDL240029243-1A_22FWHJLT4_L1_2.100bp_5prime.fq.gz
# trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired KO_5_10_CKDL240029243-1A_22FWHJLT4_L1_1.100bp_5prime.fq.gz  KO_5_10_CKDL240029243-1A_22FWHJLT4_L1_2.100bp_5prime.fq.gz
# trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired KO_6_12_CKDL240029243-1A_22FWHJLT4_L1_1.100bp_5prime.fq.gz  KO_6_12_CKDL240029243-1A_22FWHJLT4_L1_2.100bp_5prime.fq.gz
# trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired WT_1_1_CKDL240029243-1A_22FWHJLT4_L1_1.100bp_5prime.fq.gz   WT_1_1_CKDL240029243-1A_22FWHJLT4_L1_2.100bp_5prime.fq.gz
# trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired WT_2_3_CKDL240029243-1A_22FWHJLT4_L1_1.100bp_5prime.fq.gz   WT_2_3_CKDL240029243-1A_22FWHJLT4_L1_2.100bp_5prime.fq.gz
# trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired WT_3_5_CKDL240029243-1A_22FWHJLT4_L1_1.100bp_5prime.fq.gz   WT_3_5_CKDL240029243-1A_22FWHJLT4_L1_2.100bp_5prime.fq.gz
# trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired WT_4_7_CKDL240029243-1A_22FWHJLT4_L1_1.100bp_5prime.fq.gz   WT_4_7_CKDL240029243-1A_22FWHJLT4_L1_2.100bp_5prime.fq.gz
# trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired WT_5_9_CKDL240029243-1A_22FWHJLT4_L1_1.100bp_5prime.fq.gz   WT_5_9_CKDL240029243-1A_22FWHJLT4_L1_2.100bp_5prime.fq.gz
# trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired WT_6_11_CKDL240029243-1A_22FWHJLT4_L1_1.100bp_5prime.fq.gz  WT_6_11_CKDL240029243-1A_22FWHJLT4_L1_2.100bp_5prime.fq.gz

#### rename output files
mv KO_1_2_CKDL240029243-1A_22FWHJLT4_L1_1.100bp_5prime_val_1.fq.gz      PeriMac_YF_Esr1_KO_RNAseq_rep1_2_1.fastq.gz
mv KO_1_2_CKDL240029243-1A_22FWHJLT4_L1_2.100bp_5prime_val_2.fq.gz      PeriMac_YF_Esr1_KO_RNAseq_rep1_2_2.fastq.gz
mv KO_2_4_CKDL240029243-1A_22FWHJLT4_L1_1.100bp_5prime_val_1.fq.gz      PeriMac_YF_Esr1_KO_RNAseq_rep2_4_1.fastq.gz
mv KO_2_4_CKDL240029243-1A_22FWHJLT4_L1_2.100bp_5prime_val_2.fq.gz      PeriMac_YF_Esr1_KO_RNAseq_rep2_4_2.fastq.gz
mv KO_3_6_CKDL240029243-1A_22FWHJLT4_L1_1.100bp_5prime_val_1.fq.gz      PeriMac_YF_Esr1_KO_RNAseq_rep3_6_1.fastq.gz
mv KO_3_6_CKDL240029243-1A_22FWHJLT4_L1_2.100bp_5prime_val_2.fq.gz      PeriMac_YF_Esr1_KO_RNAseq_rep3_6_2.fastq.gz
mv KO_4_8_CKDL240029243-1A_22FWHJLT4_L1_1.100bp_5prime_val_1.fq.gz      PeriMac_YF_Esr1_KO_RNAseq_rep4_8_1.fastq.gz
mv KO_4_8_CKDL240029243-1A_22FWHJLT4_L1_2.100bp_5prime_val_2.fq.gz      PeriMac_YF_Esr1_KO_RNAseq_rep4_8_2.fastq.gz
mv KO_5_10_CKDL240029243-1A_22FWHJLT4_L1_1.100bp_5prime_val_1.fq.gz     PeriMac_YF_Esr1_KO_RNAseq_rep5_10_1.fastq.gz
mv KO_5_10_CKDL240029243-1A_22FWHJLT4_L1_2.100bp_5prime_val_2.fq.gz     PeriMac_YF_Esr1_KO_RNAseq_rep5_10_2.fastq.gz
mv KO_6_12_CKDL240029243-1A_22FWHJLT4_L1_1.100bp_5prime_val_1.fq.gz     PeriMac_YF_Esr1_KO_RNAseq_rep6_12_1.fastq.gz
mv KO_6_12_CKDL240029243-1A_22FWHJLT4_L1_2.100bp_5prime_val_2.fq.gz     PeriMac_YF_Esr1_KO_RNAseq_rep6_12_2.fastq.gz
mv WT_1_1_CKDL240029243-1A_22FWHJLT4_L1_1.100bp_5prime_val_1.fq.gz      PeriMac_YF_Esr1_WT_RNAseq_rep1_1_1.fastq.gz
mv WT_1_1_CKDL240029243-1A_22FWHJLT4_L1_2.100bp_5prime_val_2.fq.gz      PeriMac_YF_Esr1_WT_RNAseq_rep1_1_2.fastq.gz
mv WT_2_3_CKDL240029243-1A_22FWHJLT4_L1_1.100bp_5prime_val_1.fq.gz      PeriMac_YF_Esr1_WT_RNAseq_rep2_3_1.fastq.gz
mv WT_2_3_CKDL240029243-1A_22FWHJLT4_L1_2.100bp_5prime_val_2.fq.gz      PeriMac_YF_Esr1_WT_RNAseq_rep2_3_2.fastq.gz
mv WT_3_5_CKDL240029243-1A_22FWHJLT4_L1_1.100bp_5prime_val_1.fq.gz      PeriMac_YF_Esr1_WT_RNAseq_rep3_5_1.fastq.gz
mv WT_3_5_CKDL240029243-1A_22FWHJLT4_L1_2.100bp_5prime_val_2.fq.gz      PeriMac_YF_Esr1_WT_RNAseq_rep3_5_2.fastq.gz
mv WT_4_7_CKDL240029243-1A_22FWHJLT4_L1_1.100bp_5prime_val_1.fq.gz      PeriMac_YF_Esr1_WT_RNAseq_rep4_7_1.fastq.gz
mv WT_4_7_CKDL240029243-1A_22FWHJLT4_L1_2.100bp_5prime_val_2.fq.gz      PeriMac_YF_Esr1_WT_RNAseq_rep4_7_2.fastq.gz
mv WT_5_9_CKDL240029243-1A_22FWHJLT4_L1_1.100bp_5prime_val_1.fq.gz      PeriMac_YF_Esr1_WT_RNAseq_rep5_9_1.fastq.gz
mv WT_5_9_CKDL240029243-1A_22FWHJLT4_L1_2.100bp_5prime_val_2.fq.gz      PeriMac_YF_Esr1_WT_RNAseq_rep5_9_2.fastq.gz
mv WT_6_11_CKDL240029243-1A_22FWHJLT4_L1_1.100bp_5prime_val_1.fq.gz     PeriMac_YF_Esr1_WT_RNAseq_rep6_11_1.fastq.gz
mv WT_6_11_CKDL240029243-1A_22FWHJLT4_L1_2.100bp_5prime_val_2.fq.gz     PeriMac_YF_Esr1_WT_RNAseq_rep6_11_2.fastq.gz