## cut to 100 to help mapping with STAR
trim_galore --paired --hardtrim5 100 ../FASTQ/OldE2_1_CKDL230030340-1A_22CWNNLT3_L1_1.fq.gz       ../FASTQ/OldE2_1_CKDL230030340-1A_22CWNNLT3_L1_2.fq.gz
trim_galore --paired --hardtrim5 100 ../FASTQ/OldE2_2_CKDL230030340-1A_22CWNNLT3_L1_1.fq.gz       ../FASTQ/OldE2_2_CKDL230030340-1A_22CWNNLT3_L1_2.fq.gz
trim_galore --paired --hardtrim5 100 ../FASTQ/OldE2_3_CKDL230030340-1A_22CWNNLT3_L1_1.fq.gz       ../FASTQ/OldE2_3_CKDL230030340-1A_22CWNNLT3_L1_2.fq.gz
trim_galore --paired --hardtrim5 100 ../FASTQ/OldE2_4_CKDL230030340-1A_22CWNNLT3_L1_1.fq.gz       ../FASTQ/OldE2_4_CKDL230030340-1A_22CWNNLT3_L1_2.fq.gz
trim_galore --paired --hardtrim5 100 ../FASTQ/OldE2_5_CKDL230030340-1A_22CWNNLT3_L1_1.fq.gz       ../FASTQ/OldE2_5_CKDL230030340-1A_22CWNNLT3_L1_2.fq.gz
trim_galore --paired --hardtrim5 100 ../FASTQ/OldFemale1_CKDL230030340-1A_22CWNNLT3_L1_1.fq.gz    ../FASTQ/OldFemale1_CKDL230030340-1A_22CWNNLT3_L1_2.fq.gz
trim_galore --paired --hardtrim5 100 ../FASTQ/OldFemale2_CKDL230030340-1A_22CWNNLT3_L1_1.fq.gz    ../FASTQ/OldFemale2_CKDL230030340-1A_22CWNNLT3_L1_2.fq.gz
trim_galore --paired --hardtrim5 100 ../FASTQ/OldFemale3_CKDL230030340-1A_22CWNNLT3_L1_1.fq.gz    ../FASTQ/OldFemale3_CKDL230030340-1A_22CWNNLT3_L1_2.fq.gz
trim_galore --paired --hardtrim5 100 ../FASTQ/OldFemale4_CKDL230030340-1A_22CWNNLT3_L1_1.fq.gz    ../FASTQ/OldFemale4_CKDL230030340-1A_22CWNNLT3_L1_2.fq.gz
trim_galore --paired --hardtrim5 100 ../FASTQ/OldFemale5_CKDL230030340-1A_22CWNNLT3_L1_1.fq.gz    ../FASTQ/OldFemale5_CKDL230030340-1A_22CWNNLT3_L1_2.fq.gz
trim_galore --paired --hardtrim5 100 ../FASTQ/YoungFemale1_CKDL230030340-1A_22CWNNLT3_L1_1.fq.gz  ../FASTQ/YoungFemale1_CKDL230030340-1A_22CWNNLT3_L1_2.fq.gz
trim_galore --paired --hardtrim5 100 ../FASTQ/YoungFemale2_CKDL230030340-1A_22CWNNLT3_L1_1.fq.gz  ../FASTQ/YoungFemale2_CKDL230030340-1A_22CWNNLT3_L1_2.fq.gz
trim_galore --paired --hardtrim5 100 ../FASTQ/YoungFemale3_CKDL230030340-1A_22CWNNLT3_L1_1.fq.gz  ../FASTQ/YoungFemale3_CKDL230030340-1A_22CWNNLT3_L1_2.fq.gz
trim_galore --paired --hardtrim5 100 ../FASTQ/YoungFemale4_CKDL230030340-1A_22CWNNLT3_L1_1.fq.gz  ../FASTQ/YoungFemale4_CKDL230030340-1A_22CWNNLT3_L1_2.fq.gz
trim_galore --paired --hardtrim5 100 ../FASTQ/YoungFemale5_CKDL230030340-1A_22CWNNLT3_L1_1.fq.gz  ../FASTQ/YoungFemale5_CKDL230030340-1A_22CWNNLT3_L1_2.fq.gz

## cut 9 based on fastqc
trim_galore --stringency 5 --length 50 --clip_R1 9 --clip_R2 9 --paired OldE2_1_CKDL230030340-1A_22CWNNLT3_L1_1.100bp_5prime.fq.gz       OldE2_1_CKDL230030340-1A_22CWNNLT3_L1_2.100bp_5prime.fq.gz
trim_galore --stringency 5 --length 50 --clip_R1 9 --clip_R2 9 --paired OldE2_2_CKDL230030340-1A_22CWNNLT3_L1_1.100bp_5prime.fq.gz       OldE2_2_CKDL230030340-1A_22CWNNLT3_L1_2.100bp_5prime.fq.gz
trim_galore --stringency 5 --length 50 --clip_R1 9 --clip_R2 9 --paired OldE2_3_CKDL230030340-1A_22CWNNLT3_L1_1.100bp_5prime.fq.gz       OldE2_3_CKDL230030340-1A_22CWNNLT3_L1_2.100bp_5prime.fq.gz
trim_galore --stringency 5 --length 50 --clip_R1 9 --clip_R2 9 --paired OldE2_4_CKDL230030340-1A_22CWNNLT3_L1_1.100bp_5prime.fq.gz       OldE2_4_CKDL230030340-1A_22CWNNLT3_L1_2.100bp_5prime.fq.gz
trim_galore --stringency 5 --length 50 --clip_R1 9 --clip_R2 9 --paired OldE2_5_CKDL230030340-1A_22CWNNLT3_L1_1.100bp_5prime.fq.gz       OldE2_5_CKDL230030340-1A_22CWNNLT3_L1_2.100bp_5prime.fq.gz
trim_galore --stringency 5 --length 50 --clip_R1 9 --clip_R2 9 --paired OldFemale1_CKDL230030340-1A_22CWNNLT3_L1_1.100bp_5prime.fq.gz    OldFemale1_CKDL230030340-1A_22CWNNLT3_L1_2.100bp_5prime.fq.gz
trim_galore --stringency 5 --length 50 --clip_R1 9 --clip_R2 9 --paired OldFemale2_CKDL230030340-1A_22CWNNLT3_L1_1.100bp_5prime.fq.gz    OldFemale2_CKDL230030340-1A_22CWNNLT3_L1_2.100bp_5prime.fq.gz
trim_galore --stringency 5 --length 50 --clip_R1 9 --clip_R2 9 --paired OldFemale3_CKDL230030340-1A_22CWNNLT3_L1_1.100bp_5prime.fq.gz    OldFemale3_CKDL230030340-1A_22CWNNLT3_L1_2.100bp_5prime.fq.gz
trim_galore --stringency 5 --length 50 --clip_R1 9 --clip_R2 9 --paired OldFemale4_CKDL230030340-1A_22CWNNLT3_L1_1.100bp_5prime.fq.gz    OldFemale4_CKDL230030340-1A_22CWNNLT3_L1_2.100bp_5prime.fq.gz
trim_galore --stringency 5 --length 50 --clip_R1 9 --clip_R2 9 --paired OldFemale5_CKDL230030340-1A_22CWNNLT3_L1_1.100bp_5prime.fq.gz    OldFemale5_CKDL230030340-1A_22CWNNLT3_L1_2.100bp_5prime.fq.gz
trim_galore --stringency 5 --length 50 --clip_R1 9 --clip_R2 9 --paired YoungFemale1_CKDL230030340-1A_22CWNNLT3_L1_1.100bp_5prime.fq.gz  YoungFemale1_CKDL230030340-1A_22CWNNLT3_L1_2.100bp_5prime.fq.gz
trim_galore --stringency 5 --length 50 --clip_R1 9 --clip_R2 9 --paired YoungFemale2_CKDL230030340-1A_22CWNNLT3_L1_1.100bp_5prime.fq.gz  YoungFemale2_CKDL230030340-1A_22CWNNLT3_L1_2.100bp_5prime.fq.gz
trim_galore --stringency 5 --length 50 --clip_R1 9 --clip_R2 9 --paired YoungFemale3_CKDL230030340-1A_22CWNNLT3_L1_1.100bp_5prime.fq.gz  YoungFemale3_CKDL230030340-1A_22CWNNLT3_L1_2.100bp_5prime.fq.gz
trim_galore --stringency 5 --length 50 --clip_R1 9 --clip_R2 9 --paired YoungFemale4_CKDL230030340-1A_22CWNNLT3_L1_1.100bp_5prime.fq.gz  YoungFemale4_CKDL230030340-1A_22CWNNLT3_L1_2.100bp_5prime.fq.gz
trim_galore --stringency 5 --length 50 --clip_R1 9 --clip_R2 9 --paired YoungFemale5_CKDL230030340-1A_22CWNNLT3_L1_1.100bp_5prime.fq.gz  YoungFemale5_CKDL230030340-1A_22CWNNLT3_L1_2.100bp_5prime.fq.gz
