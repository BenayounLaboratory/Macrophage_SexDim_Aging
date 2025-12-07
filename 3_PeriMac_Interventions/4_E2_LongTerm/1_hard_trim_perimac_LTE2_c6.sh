# ## ## cut to 100 to help mapping with STAR
# trim_galore --paired --hardtrim5 100 ../FASTQ/LTE2_C6_1_CKDL250031352-1A_2325CNLT4_L2_1.fq.gz     ../FASTQ/LTE2_C6_1_CKDL250031352-1A_2325CNLT4_L2_2.fq.gz
# trim_galore --paired --hardtrim5 100 ../FASTQ/LTE2_C6_3_CKDL250031352-1A_2325CNLT4_L2_1.fq.gz     ../FASTQ/LTE2_C6_3_CKDL250031352-1A_2325CNLT4_L2_2.fq.gz
# trim_galore --paired --hardtrim5 100 ../FASTQ/LTE2_C6_5_CKDL250031352-1A_2325CNLT4_L2_1.fq.gz     ../FASTQ/LTE2_C6_5_CKDL250031352-1A_2325CNLT4_L2_2.fq.gz
# trim_galore --paired --hardtrim5 100 ../FASTQ/LTE2_C6_6_CKDL250031352-1A_2325CNLT4_L2_1.fq.gz     ../FASTQ/LTE2_C6_6_CKDL250031352-1A_2325CNLT4_L2_2.fq.gz
# trim_galore --paired --hardtrim5 100 ../FASTQ/LTE2_C6_7_CKDL250031352-1A_2325CNLT4_L2_1.fq.gz     ../FASTQ/LTE2_C6_7_CKDL250031352-1A_2325CNLT4_L2_2.fq.gz
# trim_galore --paired --hardtrim5 100 ../FASTQ/LTE2_C6_8_CKDL250031352-1A_2325CNLT4_L2_1.fq.gz     ../FASTQ/LTE2_C6_8_CKDL250031352-1A_2325CNLT4_L2_2.fq.gz
# trim_galore --paired --hardtrim5 100 ../FASTQ/LTE2_C6_9_CKDL250031352-1A_2325CNLT4_L2_1.fq.gz     ../FASTQ/LTE2_C6_9_CKDL250031352-1A_2325CNLT4_L2_2.fq.gz
# trim_galore --paired --hardtrim5 100 ../FASTQ/LTE2_C6_10_CKDL250031352-1A_2325CNLT4_L2_1.fq.gz    ../FASTQ/LTE2_C6_10_CKDL250031352-1A_2325CNLT4_L2_2.fq.gz
# trim_galore --paired --hardtrim5 100 ../FASTQ/LTE2_C6_11_CKDL250031352-1A_2325CNLT4_L2_1.fq.gz    ../FASTQ/LTE2_C6_11_CKDL250031352-1A_2325CNLT4_L2_2.fq.gz
# trim_galore --paired --hardtrim5 100 ../FASTQ/LTE2_C6_12_CKDL250031352-1A_2325CNLT4_L2_1.fq.gz    ../FASTQ/LTE2_C6_12_CKDL250031352-1A_2325CNLT4_L2_2.fq.gz
# trim_galore --paired --hardtrim5 100 ../FASTQ/LTE2_C6_13_CKDL250031352-1A_2325CNLT4_L2_1.fq.gz    ../FASTQ/LTE2_C6_13_CKDL250031352-1A_2325CNLT4_L2_2.fq.gz
# trim_galore --paired --hardtrim5 100 ../FASTQ/LTE2_C6_14_CKDL250031352-1A_2325CNLT4_L2_1.fq.gz    ../FASTQ/LTE2_C6_14_CKDL250031352-1A_2325CNLT4_L2_2.fq.gz
# trim_galore --paired --hardtrim5 100 ../FASTQ/LTE2_C6_15_CKDL250031352-1A_2325CNLT4_L2_1.fq.gz    ../FASTQ/LTE2_C6_15_CKDL250031352-1A_2325CNLT4_L2_2.fq.gz

## ### cut 9 based on fastqc
trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired LTE2_C6_1_CKDL250031352-1A_2325CNLT4_L2_1.100bp_5prime.fq.gz     LTE2_C6_1_CKDL250031352-1A_2325CNLT4_L2_2.100bp_5prime.fq.gz
trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired LTE2_C6_3_CKDL250031352-1A_2325CNLT4_L2_1.100bp_5prime.fq.gz     LTE2_C6_3_CKDL250031352-1A_2325CNLT4_L2_2.100bp_5prime.fq.gz
trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired LTE2_C6_5_CKDL250031352-1A_2325CNLT4_L2_1.100bp_5prime.fq.gz     LTE2_C6_5_CKDL250031352-1A_2325CNLT4_L2_2.100bp_5prime.fq.gz
trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired LTE2_C6_6_CKDL250031352-1A_2325CNLT4_L2_1.100bp_5prime.fq.gz     LTE2_C6_6_CKDL250031352-1A_2325CNLT4_L2_2.100bp_5prime.fq.gz
trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired LTE2_C6_7_CKDL250031352-1A_2325CNLT4_L2_1.100bp_5prime.fq.gz     LTE2_C6_7_CKDL250031352-1A_2325CNLT4_L2_2.100bp_5prime.fq.gz
trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired LTE2_C6_8_CKDL250031352-1A_2325CNLT4_L2_1.100bp_5prime.fq.gz     LTE2_C6_8_CKDL250031352-1A_2325CNLT4_L2_2.100bp_5prime.fq.gz
trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired LTE2_C6_9_CKDL250031352-1A_2325CNLT4_L2_1.100bp_5prime.fq.gz     LTE2_C6_9_CKDL250031352-1A_2325CNLT4_L2_2.100bp_5prime.fq.gz
trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired LTE2_C6_10_CKDL250031352-1A_2325CNLT4_L2_1.100bp_5prime.fq.gz    LTE2_C6_10_CKDL250031352-1A_2325CNLT4_L2_2.100bp_5prime.fq.gz
trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired LTE2_C6_11_CKDL250031352-1A_2325CNLT4_L2_1.100bp_5prime.fq.gz    LTE2_C6_11_CKDL250031352-1A_2325CNLT4_L2_2.100bp_5prime.fq.gz
trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired LTE2_C6_12_CKDL250031352-1A_2325CNLT4_L2_1.100bp_5prime.fq.gz    LTE2_C6_12_CKDL250031352-1A_2325CNLT4_L2_2.100bp_5prime.fq.gz
trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired LTE2_C6_13_CKDL250031352-1A_2325CNLT4_L2_1.100bp_5prime.fq.gz    LTE2_C6_13_CKDL250031352-1A_2325CNLT4_L2_2.100bp_5prime.fq.gz
trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired LTE2_C6_14_CKDL250031352-1A_2325CNLT4_L2_1.100bp_5prime.fq.gz    LTE2_C6_14_CKDL250031352-1A_2325CNLT4_L2_2.100bp_5prime.fq.gz
trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired LTE2_C6_15_CKDL250031352-1A_2325CNLT4_L2_1.100bp_5prime.fq.gz    LTE2_C6_15_CKDL250031352-1A_2325CNLT4_L2_2.100bp_5prime.fq.gz


#### # rename output files
mv LTE2_C6_1_CKDL250031352-1A_2325CNLT4_L2_1.100bp_5prime_val_1.fq.gz   PeriMac_LTE2_01_O_E2_1.fastq.gz 
mv LTE2_C6_3_CKDL250031352-1A_2325CNLT4_L2_1.100bp_5prime_val_1.fq.gz   PeriMac_LTE2_03_Y_Ve_1.fastq.gz 
mv LTE2_C6_5_CKDL250031352-1A_2325CNLT4_L2_1.100bp_5prime_val_1.fq.gz   PeriMac_LTE2_05_O_Ve_1.fastq.gz 
mv LTE2_C6_6_CKDL250031352-1A_2325CNLT4_L2_1.100bp_5prime_val_1.fq.gz   PeriMac_LTE2_06_Y_Ve_1.fastq.gz 
mv LTE2_C6_7_CKDL250031352-1A_2325CNLT4_L2_1.100bp_5prime_val_1.fq.gz   PeriMac_LTE2_07_O_E2_1.fastq.gz 
mv LTE2_C6_8_CKDL250031352-1A_2325CNLT4_L2_1.100bp_5prime_val_1.fq.gz   PeriMac_LTE2_08_O_Ve_1.fastq.gz 
mv LTE2_C6_9_CKDL250031352-1A_2325CNLT4_L2_1.100bp_5prime_val_1.fq.gz   PeriMac_LTE2_09_Y_Ve_1.fastq.gz 
mv LTE2_C6_10_CKDL250031352-1A_2325CNLT4_L2_1.100bp_5prime_val_1.fq.gz  PeriMac_LTE2_10_O_E2_1.fastq.gz 
mv LTE2_C6_11_CKDL250031352-1A_2325CNLT4_L2_1.100bp_5prime_val_1.fq.gz  PeriMac_LTE2_11_O_Ve_1.fastq.gz 
mv LTE2_C6_12_CKDL250031352-1A_2325CNLT4_L2_1.100bp_5prime_val_1.fq.gz  PeriMac_LTE2_12_Y_Ve_1.fastq.gz 
mv LTE2_C6_13_CKDL250031352-1A_2325CNLT4_L2_1.100bp_5prime_val_1.fq.gz  PeriMac_LTE2_13_O_E2_1.fastq.gz 
mv LTE2_C6_14_CKDL250031352-1A_2325CNLT4_L2_1.100bp_5prime_val_1.fq.gz  PeriMac_LTE2_14_O_Ve_1.fastq.gz 
mv LTE2_C6_15_CKDL250031352-1A_2325CNLT4_L2_1.100bp_5prime_val_1.fq.gz  PeriMac_LTE2_15_Y_Ve_1.fastq.gz 
  
mv LTE2_C6_1_CKDL250031352-1A_2325CNLT4_L2_2.100bp_5prime_val_2.fq.gz   PeriMac_LTE2_01_O_E2_2.fastq.gz 
mv LTE2_C6_3_CKDL250031352-1A_2325CNLT4_L2_2.100bp_5prime_val_2.fq.gz   PeriMac_LTE2_03_Y_Ve_2.fastq.gz 
mv LTE2_C6_5_CKDL250031352-1A_2325CNLT4_L2_2.100bp_5prime_val_2.fq.gz   PeriMac_LTE2_05_O_Ve_2.fastq.gz 
mv LTE2_C6_6_CKDL250031352-1A_2325CNLT4_L2_2.100bp_5prime_val_2.fq.gz   PeriMac_LTE2_06_Y_Ve_2.fastq.gz 
mv LTE2_C6_7_CKDL250031352-1A_2325CNLT4_L2_2.100bp_5prime_val_2.fq.gz   PeriMac_LTE2_07_O_E2_2.fastq.gz 
mv LTE2_C6_8_CKDL250031352-1A_2325CNLT4_L2_2.100bp_5prime_val_2.fq.gz   PeriMac_LTE2_08_O_Ve_2.fastq.gz 
mv LTE2_C6_9_CKDL250031352-1A_2325CNLT4_L2_2.100bp_5prime_val_2.fq.gz   PeriMac_LTE2_09_Y_Ve_2.fastq.gz 
mv LTE2_C6_10_CKDL250031352-1A_2325CNLT4_L2_2.100bp_5prime_val_2.fq.gz  PeriMac_LTE2_10_O_E2_2.fastq.gz 
mv LTE2_C6_11_CKDL250031352-1A_2325CNLT4_L2_2.100bp_5prime_val_2.fq.gz  PeriMac_LTE2_11_O_Ve_2.fastq.gz 
mv LTE2_C6_12_CKDL250031352-1A_2325CNLT4_L2_2.100bp_5prime_val_2.fq.gz  PeriMac_LTE2_12_Y_Ve_2.fastq.gz 
mv LTE2_C6_13_CKDL250031352-1A_2325CNLT4_L2_2.100bp_5prime_val_2.fq.gz  PeriMac_LTE2_13_O_E2_2.fastq.gz 
mv LTE2_C6_14_CKDL250031352-1A_2325CNLT4_L2_2.100bp_5prime_val_2.fq.gz  PeriMac_LTE2_14_O_Ve_2.fastq.gz 
mv LTE2_C6_15_CKDL250031352-1A_2325CNLT4_L2_2.100bp_5prime_val_2.fq.gz  PeriMac_LTE2_15_Y_Ve_2.fastq.gz 
