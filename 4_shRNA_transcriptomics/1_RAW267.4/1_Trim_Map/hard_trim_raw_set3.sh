# ### use trimgalore to recapitulate the behavior of fastx_trimmer (doesn't run on newer macos)
# ## cut to 100 to help mapping with STAR
# trim_galore --paired --hardtrim5 100 ../FASTQ/SET3/Luc7_1.fq.gz ../FASTQ/SET3/Luc7_2.fq.gz
# trim_galore --paired --hardtrim5 100 ../FASTQ/SET3/Luc8_1.fq.gz ../FASTQ/SET3/Luc8_2.fq.gz
# trim_galore --paired --hardtrim5 100 ../FASTQ/SET3/Luc9_1.fq.gz ../FASTQ/SET3/Luc9_2.fq.gz
# trim_galore --paired --hardtrim5 100 ../FASTQ/SET3/Mef2c69_1_1.fq.gz ../FASTQ/SET3/Mef2c69_1_2.fq.gz
# trim_galore --paired --hardtrim5 100 ../FASTQ/SET3/Mef2c69_2_1.fq.gz ../FASTQ/SET3/Mef2c69_2_2.fq.gz
# trim_galore --paired --hardtrim5 100 ../FASTQ/SET3/Mef2c69_3_1.fq.gz ../FASTQ/SET3/Mef2c69_3_2.fq.gz
# trim_galore --paired --hardtrim5 100 ../FASTQ/SET3/Mef2c72_1_1.fq.gz ../FASTQ/SET3/Mef2c72_1_2.fq.gz
# trim_galore --paired --hardtrim5 100 ../FASTQ/SET3/Mef2c72_2_1.fq.gz ../FASTQ/SET3/Mef2c72_2_2.fq.gz
# trim_galore --paired --hardtrim5 100 ../FASTQ/SET3/Mef2c72_3_1.fq.gz ../FASTQ/SET3/Mef2c72_3_2.fq.gz
# trim_galore --paired --hardtrim5 100 ../FASTQ/SET3/Meis1cds1_1.fq.gz ../FASTQ/SET3/Meis1cds1_2.fq.gz
# trim_galore --paired --hardtrim5 100 ../FASTQ/SET3/Meis1cds2_1.fq.gz ../FASTQ/SET3/Meis1cds2_2.fq.gz
# trim_galore --paired --hardtrim5 100 ../FASTQ/SET3/Meis1cds3_1.fq.gz ../FASTQ/SET3/Meis1cds3_2.fq.gz
# trim_galore --paired --hardtrim5 100 ../FASTQ/SET3/Meis13utr1_1.fq.gz ../FASTQ/SET3/Meis13utr1_2.fq.gz
# trim_galore --paired --hardtrim5 100 ../FASTQ/SET3/Meis13utr2_1.fq.gz ../FASTQ/SET3/Meis13utr2_2.fq.gz
# trim_galore --paired --hardtrim5 100 ../FASTQ/SET3/Meis13utr3_1.fq.gz ../FASTQ/SET3/Meis13utr3_2.fq.gz
# trim_galore --paired --hardtrim5 100 ../FASTQ/SET3/Non7_1.fq.gz ../FASTQ/SET3/Non7_2.fq.gz
# trim_galore --paired --hardtrim5 100 ../FASTQ/SET3/Non8_1.fq.gz ../FASTQ/SET3/Non8_2.fq.gz
# trim_galore --paired --hardtrim5 100 ../FASTQ/SET3/Non9_1.fq.gz ../FASTQ/SET3/Non9_2.fq.gz
# 
# ### cut 9 based on fastqc
# trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired Luc7_1.100bp_5prime.fq.gz Luc7_2.100bp_5prime.fq.gz
# trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired Luc8_1.100bp_5prime.fq.gz Luc8_2.100bp_5prime.fq.gz
# trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired Luc9_1.100bp_5prime.fq.gz Luc9_2.100bp_5prime.fq.gz
# trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired Mef2c69_1_1.100bp_5prime.fq.gz Mef2c69_1_2.100bp_5prime.fq.gz
# trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired Mef2c69_2_1.100bp_5prime.fq.gz Mef2c69_2_2.100bp_5prime.fq.gz
# trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired Mef2c69_3_1.100bp_5prime.fq.gz Mef2c69_3_2.100bp_5prime.fq.gz
# trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired Mef2c72_1_1.100bp_5prime.fq.gz Mef2c72_1_2.100bp_5prime.fq.gz
# trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired Mef2c72_2_1.100bp_5prime.fq.gz Mef2c72_2_2.100bp_5prime.fq.gz
# trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired Mef2c72_3_1.100bp_5prime.fq.gz Mef2c72_3_2.100bp_5prime.fq.gz
# trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired Meis1cds1_1.100bp_5prime.fq.gz Meis1cds1_2.100bp_5prime.fq.gz
# trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired Meis1cds2_1.100bp_5prime.fq.gz Meis1cds2_2.100bp_5prime.fq.gz
# trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired Meis1cds3_1.100bp_5prime.fq.gz Meis1cds3_2.100bp_5prime.fq.gz
# trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired Meis13utr1_1.100bp_5prime.fq.gz Meis13utr1_2.100bp_5prime.fq.gz
# trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired Meis13utr2_1.100bp_5prime.fq.gz Meis13utr2_2.100bp_5prime.fq.gz
# trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired Meis13utr3_1.100bp_5prime.fq.gz Meis13utr3_2.100bp_5prime.fq.gz
# trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired Non7_1.100bp_5prime.fq.gz Non7_2.100bp_5prime.fq.gz
# trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired Non8_1.100bp_5prime.fq.gz Non8_2.100bp_5prime.fq.gz
# trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired Non9_1.100bp_5prime.fq.gz Non9_2.100bp_5prime.fq.gz

# rename output files
mv Luc7_1.100bp_5prime_val_1.fq.gz        shLuc_rep7_RNAseq_SET3_1.fastq.gz
mv Luc7_2.100bp_5prime_val_2.fq.gz        shLuc_rep7_RNAseq_SET3_2.fastq.gz
mv Luc8_1.100bp_5prime_val_1.fq.gz        shLuc_rep8_RNAseq_SET3_1.fastq.gz
mv Luc8_2.100bp_5prime_val_2.fq.gz        shLuc_rep8_RNAseq_SET3_2.fastq.gz
mv Luc9_1.100bp_5prime_val_1.fq.gz        shLuc_rep9_RNAseq_SET3_1.fastq.gz
mv Luc9_2.100bp_5prime_val_2.fq.gz        shLuc_rep9_RNAseq_SET3_2.fastq.gz
mv Mef2c69_1_1.100bp_5prime_val_1.fq.gz   shMef2c69_rep1_RNAseq_SET3_1.fastq.gz
mv Mef2c69_1_2.100bp_5prime_val_2.fq.gz   shMef2c69_rep1_RNAseq_SET3_2.fastq.gz
mv Mef2c69_2_1.100bp_5prime_val_1.fq.gz   shMef2c69_rep2_RNAseq_SET3_1.fastq.gz
mv Mef2c69_2_2.100bp_5prime_val_2.fq.gz   shMef2c69_rep2_RNAseq_SET3_2.fastq.gz
mv Mef2c69_3_1.100bp_5prime_val_1.fq.gz   shMef2c69_rep3_RNAseq_SET3_1.fastq.gz
mv Mef2c69_3_2.100bp_5prime_val_2.fq.gz   shMef2c69_rep3_RNAseq_SET3_2.fastq.gz
mv Mef2c72_1_1.100bp_5prime_val_1.fq.gz   shMef2c72_rep1_RNAseq_SET3_1.fastq.gz
mv Mef2c72_1_2.100bp_5prime_val_2.fq.gz   shMef2c72_rep1_RNAseq_SET3_2.fastq.gz
mv Mef2c72_2_1.100bp_5prime_val_1.fq.gz   shMef2c72_rep2_RNAseq_SET3_1.fastq.gz
mv Mef2c72_2_2.100bp_5prime_val_2.fq.gz   shMef2c72_rep2_RNAseq_SET3_2.fastq.gz
mv Mef2c72_3_1.100bp_5prime_val_1.fq.gz   shMef2c72_rep3_RNAseq_SET3_1.fastq.gz
mv Mef2c72_3_2.100bp_5prime_val_2.fq.gz   shMef2c72_rep3_RNAseq_SET3_2.fastq.gz
mv Meis1cds1_1.100bp_5prime_val_1.fq.gz   shMeis1cds_rep1_RNAseq_SET3_1.fastq.gz
mv Meis1cds1_2.100bp_5prime_val_2.fq.gz   shMeis1cds_rep1_RNAseq_SET3_2.fastq.gz
mv Meis1cds2_1.100bp_5prime_val_1.fq.gz   shMeis1cds_rep2_RNAseq_SET3_1.fastq.gz
mv Meis1cds2_2.100bp_5prime_val_2.fq.gz   shMeis1cds_rep2_RNAseq_SET3_2.fastq.gz
mv Meis1cds3_1.100bp_5prime_val_1.fq.gz   shMeis1cds_rep3_RNAseq_SET3_1.fastq.gz
mv Meis1cds3_2.100bp_5prime_val_2.fq.gz   shMeis1cds_rep3_RNAseq_SET3_2.fastq.gz
mv Meis13utr1_1.100bp_5prime_val_1.fq.gz  shMeis13utr_rep1_RNAseq_SET3_1.fastq.gz
mv Meis13utr1_2.100bp_5prime_val_2.fq.gz  shMeis13utr_rep1_RNAseq_SET3_2.fastq.gz
mv Meis13utr2_1.100bp_5prime_val_1.fq.gz  shMeis13utr_rep2_RNAseq_SET3_1.fastq.gz
mv Meis13utr2_2.100bp_5prime_val_2.fq.gz  shMeis13utr_rep2_RNAseq_SET3_2.fastq.gz
mv Meis13utr3_1.100bp_5prime_val_1.fq.gz  shMeis13utr_rep3_RNAseq_SET3_1.fastq.gz
mv Meis13utr3_2.100bp_5prime_val_2.fq.gz  shMeis13utr_rep3_RNAseq_SET3_2.fastq.gz
mv Non7_1.100bp_5prime_val_1.fq.gz        shNon_rep7_RNAseq_SET3_1.fastq.gz
mv Non7_2.100bp_5prime_val_2.fq.gz        shNon_rep7_RNAseq_SET3_2.fastq.gz
mv Non8_1.100bp_5prime_val_1.fq.gz        shNon_rep8_RNAseq_SET3_1.fastq.gz
mv Non8_2.100bp_5prime_val_2.fq.gz        shNon_rep8_RNAseq_SET3_2.fastq.gz
mv Non9_1.100bp_5prime_val_1.fq.gz        shNon_rep9_RNAseq_SET3_1.fastq.gz
mv Non9_2.100bp_5prime_val_2.fq.gz        shNon_rep9_RNAseq_SET3_2.fastq.gz