# ### use trimgalore to recapitulate the behavior of fastx_trimmer (doesn't run on newer macos)
# ## cut to 100 to help mapping with STAR
trim_galore --paired --hardtrim5 100 ../FASTQ/Lane2/KO_1_2_CKDL240033869-1A_22HVYLLT4_L1_1.fq.gz  ../FASTQ/Lane2/KO_1_2_CKDL240033869-1A_22HVYLLT4_L1_2.fq.gz
trim_galore --paired --hardtrim5 100 ../FASTQ/Lane2/KO_2_4_CKDL240033869-1A_22HVYLLT4_L1_1.fq.gz  ../FASTQ/Lane2/KO_2_4_CKDL240033869-1A_22HVYLLT4_L1_2.fq.gz
trim_galore --paired --hardtrim5 100 ../FASTQ/Lane2/KO_3_6_CKDL240033870-1A_22HVYLLT4_L1_1.fq.gz  ../FASTQ/Lane2/KO_3_6_CKDL240033870-1A_22HVYLLT4_L1_2.fq.gz
trim_galore --paired --hardtrim5 100 ../FASTQ/Lane2/WT_1_1_CKDL240033869-1A_22HVYLLT4_L1_1.fq.gz  ../FASTQ/Lane2/WT_1_1_CKDL240033869-1A_22HVYLLT4_L1_2.fq.gz
trim_galore --paired --hardtrim5 100 ../FASTQ/Lane2/WT_2_3_CKDL240033869-1A_22HVYLLT4_L1_1.fq.gz  ../FASTQ/Lane2/WT_2_3_CKDL240033869-1A_22HVYLLT4_L1_2.fq.gz
trim_galore --paired --hardtrim5 100 ../FASTQ/Lane2/WT_4_7_CKDL240033870-1A_22HVYLLT4_L1_1.fq.gz  ../FASTQ/Lane2/WT_4_7_CKDL240033870-1A_22HVYLLT4_L1_2.fq.gz
trim_galore --paired --hardtrim5 100 ../FASTQ/Lane2/WT_5_9_CKDL240033870-1A_22HVYLLT4_L1_1.fq.gz  ../FASTQ/Lane2/WT_5_9_CKDL240033870-1A_22HVYLLT4_L1_2.fq.gz

# # ### cut 9 based on fastqc
trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired KO_1_2_CKDL240033869-1A_22HVYLLT4_L1_1.100bp_5prime.fq.gz  KO_1_2_CKDL240033869-1A_22HVYLLT4_L1_2.100bp_5prime.fq.gz
trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired KO_2_4_CKDL240033869-1A_22HVYLLT4_L1_1.100bp_5prime.fq.gz  KO_2_4_CKDL240033869-1A_22HVYLLT4_L1_2.100bp_5prime.fq.gz
trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired KO_3_6_CKDL240033870-1A_22HVYLLT4_L1_1.100bp_5prime.fq.gz  KO_3_6_CKDL240033870-1A_22HVYLLT4_L1_2.100bp_5prime.fq.gz
trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired WT_1_1_CKDL240033869-1A_22HVYLLT4_L1_1.100bp_5prime.fq.gz  WT_1_1_CKDL240033869-1A_22HVYLLT4_L1_2.100bp_5prime.fq.gz
trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired WT_2_3_CKDL240033869-1A_22HVYLLT4_L1_1.100bp_5prime.fq.gz  WT_2_3_CKDL240033869-1A_22HVYLLT4_L1_2.100bp_5prime.fq.gz
trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired WT_4_7_CKDL240033870-1A_22HVYLLT4_L1_1.100bp_5prime.fq.gz  WT_4_7_CKDL240033870-1A_22HVYLLT4_L1_2.100bp_5prime.fq.gz
trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired WT_5_9_CKDL240033870-1A_22HVYLLT4_L1_1.100bp_5prime.fq.gz  WT_5_9_CKDL240033870-1A_22HVYLLT4_L1_2.100bp_5prime.fq.gz



#### rename output files
mv KO_1_2_CKDL240033869-1A_22HVYLLT4_L1_1.100bp_5prime_val_1.fq.gz   PeriMac_YF_Esr1_KO_RNAseq_DEEPER_rep1_2_1.fastq.gz
mv KO_1_2_CKDL240033869-1A_22HVYLLT4_L1_2.100bp_5prime_val_2.fq.gz   PeriMac_YF_Esr1_KO_RNAseq_DEEPER_rep1_2_2.fastq.gz
mv KO_2_4_CKDL240033869-1A_22HVYLLT4_L1_1.100bp_5prime_val_1.fq.gz   PeriMac_YF_Esr1_KO_RNAseq_DEEPER_rep2_4_1.fastq.gz
mv KO_2_4_CKDL240033869-1A_22HVYLLT4_L1_2.100bp_5prime_val_2.fq.gz   PeriMac_YF_Esr1_KO_RNAseq_DEEPER_rep2_4_2.fastq.gz
mv KO_3_6_CKDL240033870-1A_22HVYLLT4_L1_1.100bp_5prime_val_1.fq.gz   PeriMac_YF_Esr1_KO_RNAseq_DEEPER_rep3_6_1.fastq.gz
mv KO_3_6_CKDL240033870-1A_22HVYLLT4_L1_2.100bp_5prime_val_2.fq.gz   PeriMac_YF_Esr1_KO_RNAseq_DEEPER_rep3_6_2.fastq.gz
mv WT_1_1_CKDL240033869-1A_22HVYLLT4_L1_1.100bp_5prime_val_1.fq.gz   PeriMac_YF_Esr1_WT_RNAseq_DEEPER_rep1_1_1.fastq.gz
mv WT_1_1_CKDL240033869-1A_22HVYLLT4_L1_2.100bp_5prime_val_2.fq.gz   PeriMac_YF_Esr1_WT_RNAseq_DEEPER_rep1_1_2.fastq.gz
mv WT_2_3_CKDL240033869-1A_22HVYLLT4_L1_1.100bp_5prime_val_1.fq.gz   PeriMac_YF_Esr1_WT_RNAseq_DEEPER_rep2_3_1.fastq.gz
mv WT_2_3_CKDL240033869-1A_22HVYLLT4_L1_2.100bp_5prime_val_2.fq.gz   PeriMac_YF_Esr1_WT_RNAseq_DEEPER_rep2_3_2.fastq.gz
mv WT_4_7_CKDL240033870-1A_22HVYLLT4_L1_1.100bp_5prime_val_1.fq.gz   PeriMac_YF_Esr1_WT_RNAseq_DEEPER_rep4_7_1.fastq.gz
mv WT_4_7_CKDL240033870-1A_22HVYLLT4_L1_2.100bp_5prime_val_2.fq.gz   PeriMac_YF_Esr1_WT_RNAseq_DEEPER_rep4_7_2.fastq.gz
mv WT_5_9_CKDL240033870-1A_22HVYLLT4_L1_1.100bp_5prime_val_1.fq.gz   PeriMac_YF_Esr1_WT_RNAseq_DEEPER_rep5_9_1.fastq.gz
mv WT_5_9_CKDL240033870-1A_22HVYLLT4_L1_2.100bp_5prime_val_2.fq.gz   PeriMac_YF_Esr1_WT_RNAseq_DEEPER_rep5_9_2.fastq.gz
