# ### use trimgalore to recapitulate the behavior of fastx_trimmer (doesn't run on newer macos)
# ## cut to 100 to help mapping with STAR

trim_galore --paired --hardtrim5 100 ../FASTQ/SET4/Luc1_1.fq.gz ../FASTQ/SET4/Luc1_2.fq.gz
trim_galore --paired --hardtrim5 100 ../FASTQ/SET4/Luc2_1.fq.gz ../FASTQ/SET4/Luc2_2.fq.gz
trim_galore --paired --hardtrim5 100 ../FASTQ/SET4/Luc3_1.fq.gz ../FASTQ/SET4/Luc3_2.fq.gz
trim_galore --paired --hardtrim5 100 ../FASTQ/SET4/Mef2c_69_1_1.fq.gz ../FASTQ/SET4/Mef2c_69_1_2.fq.gz
trim_galore --paired --hardtrim5 100 ../FASTQ/SET4/Mef2c_69_2_1.fq.gz ../FASTQ/SET4/Mef2c_69_2_2.fq.gz
trim_galore --paired --hardtrim5 100 ../FASTQ/SET4/Mef2c_69_3_1.fq.gz ../FASTQ/SET4/Mef2c_69_3_2.fq.gz
trim_galore --paired --hardtrim5 100 ../FASTQ/SET4/Mef2c_72_1_1.fq.gz ../FASTQ/SET4/Mef2c_72_1_2.fq.gz
trim_galore --paired --hardtrim5 100 ../FASTQ/SET4/Mef2c_72_2_1.fq.gz ../FASTQ/SET4/Mef2c_72_2_2.fq.gz
trim_galore --paired --hardtrim5 100 ../FASTQ/SET4/Mef2c_72_3_1.fq.gz ../FASTQ/SET4/Mef2c_72_3_2.fq.gz
trim_galore --paired --hardtrim5 100 ../FASTQ/SET4/Non1_1.fq.gz ../FASTQ/SET4/Non1_2.fq.gz
trim_galore --paired --hardtrim5 100 ../FASTQ/SET4/Non2_1.fq.gz ../FASTQ/SET4/Non2_2.fq.gz
trim_galore --paired --hardtrim5 100 ../FASTQ/SET4/Non3_1.fq.gz ../FASTQ/SET4/Non3_2.fq.gz


# ### cut 9 based on fastqc
trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired Luc1_1.100bp_5prime.fq.gz Luc1_2.100bp_5prime.fq.gz
trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired Luc2_1.100bp_5prime.fq.gz Luc2_2.100bp_5prime.fq.gz
trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired Luc3_1.100bp_5prime.fq.gz Luc3_2.100bp_5prime.fq.gz
trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired Mef2c_69_1_1.100bp_5prime.fq.gz Mef2c_69_1_2.100bp_5prime.fq.gz
trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired Mef2c_69_2_1.100bp_5prime.fq.gz Mef2c_69_2_2.100bp_5prime.fq.gz
trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired Mef2c_69_3_1.100bp_5prime.fq.gz Mef2c_69_3_2.100bp_5prime.fq.gz
trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired Mef2c_72_1_1.100bp_5prime.fq.gz Mef2c_72_1_2.100bp_5prime.fq.gz
trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired Mef2c_72_2_1.100bp_5prime.fq.gz Mef2c_72_2_2.100bp_5prime.fq.gz
trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired Mef2c_72_3_1.100bp_5prime.fq.gz Mef2c_72_3_2.100bp_5prime.fq.gz
trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired Non1_1.100bp_5prime.fq.gz Non1_2.100bp_5prime.fq.gz
trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired Non2_1.100bp_5prime.fq.gz Non2_2.100bp_5prime.fq.gz
trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired Non3_1.100bp_5prime.fq.gz Non3_2.100bp_5prime.fq.gz


# rename output files
mv Luc1_1.100bp_5prime_val_1.fq.gz        shLuc_rep10_RNAseq_SET4_1.fastq.gz
mv Luc1_2.100bp_5prime_val_2.fq.gz        shLuc_rep10_RNAseq_SET4_2.fastq.gz
mv Luc2_1.100bp_5prime_val_1.fq.gz        shLuc_rep11_RNAseq_SET4_1.fastq.gz
mv Luc2_2.100bp_5prime_val_2.fq.gz        shLuc_rep11_RNAseq_SET4_2.fastq.gz
mv Luc3_1.100bp_5prime_val_1.fq.gz        shLuc_rep12_RNAseq_SET4_1.fastq.gz
mv Luc3_2.100bp_5prime_val_2.fq.gz        shLuc_rep12_RNAseq_SET4_2.fastq.gz
mv Mef2c_69_1_1.100bp_5prime_val_1.fq.gz   shMef2c_69_rep1_RNAseq_SET4_1.fastq.gz
mv Mef2c_69_1_2.100bp_5prime_val_2.fq.gz   shMef2c_69_rep1_RNAseq_SET4_2.fastq.gz
mv Mef2c_69_2_1.100bp_5prime_val_1.fq.gz   shMef2c_69_rep2_RNAseq_SET4_1.fastq.gz
mv Mef2c_69_2_2.100bp_5prime_val_2.fq.gz   shMef2c_69_rep2_RNAseq_SET4_2.fastq.gz
mv Mef2c_69_3_1.100bp_5prime_val_1.fq.gz   shMef2c_69_rep3_RNAseq_SET4_1.fastq.gz
mv Mef2c_69_3_2.100bp_5prime_val_2.fq.gz   shMef2c_69_rep3_RNAseq_SET4_2.fastq.gz
mv Mef2c_72_1_1.100bp_5prime_val_1.fq.gz   shMef2c_72_rep1_RNAseq_SET4_1.fastq.gz
mv Mef2c_72_1_2.100bp_5prime_val_2.fq.gz   shMef2c_72_rep1_RNAseq_SET4_2.fastq.gz
mv Mef2c_72_2_1.100bp_5prime_val_1.fq.gz   shMef2c_72_rep2_RNAseq_SET4_1.fastq.gz
mv Mef2c_72_2_2.100bp_5prime_val_2.fq.gz   shMef2c_72_rep2_RNAseq_SET4_2.fastq.gz
mv Mef2c_72_3_1.100bp_5prime_val_1.fq.gz   shMef2c_72_rep3_RNAseq_SET4_1.fastq.gz
mv Mef2c_72_3_2.100bp_5prime_val_2.fq.gz   shMef2c_72_rep3_RNAseq_SET4_2.fastq.gz
mv Non1_1.100bp_5prime_val_1.fq.gz        shNon_rep10_RNAseq_SET4_1.fastq.gz
mv Non1_2.100bp_5prime_val_2.fq.gz        shNon_rep10_RNAseq_SET4_2.fastq.gz
mv Non2_1.100bp_5prime_val_1.fq.gz        shNon_rep11_RNAseq_SET4_1.fastq.gz
mv Non2_2.100bp_5prime_val_2.fq.gz        shNon_rep11_RNAseq_SET4_2.fastq.gz
mv Non3_1.100bp_5prime_val_1.fq.gz        shNon_rep12_RNAseq_SET4_1.fastq.gz
mv Non3_2.100bp_5prime_val_2.fq.gz        shNon_rep12_RNAseq_SET4_2.fastq.gz