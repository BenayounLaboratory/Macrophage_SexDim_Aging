# ### use trimgalore to recapitulate the behavior of fastx_trimmer (doesn't run on newer macos)
# ## cut to 100 to help mapping with STAR
trim_galore --paired --hardtrim5 100 ../FASTQ/Irf2_3utr1_1.fq.gz   ../FASTQ/Irf2_3utr1_2.fq.gz
trim_galore --paired --hardtrim5 100 ../FASTQ/Irf2_3utr2_1.fq.gz   ../FASTQ/Irf2_3utr2_2.fq.gz
trim_galore --paired --hardtrim5 100 ../FASTQ/Irf2_3utr3_1.fq.gz   ../FASTQ/Irf2_3utr3_2.fq.gz
trim_galore --paired --hardtrim5 100 ../FASTQ/Irf2_cds1_1.fq.gz    ../FASTQ/Irf2_cds1_2.fq.gz
trim_galore --paired --hardtrim5 100 ../FASTQ/Irf2_cds2_1.fq.gz    ../FASTQ/Irf2_cds2_2.fq.gz
trim_galore --paired --hardtrim5 100 ../FASTQ/Irf2_cds3_1.fq.gz    ../FASTQ/Irf2_cds3_2.fq.gz
trim_galore --paired --hardtrim5 100 ../FASTQ/Luc1_1.fq.gz         ../FASTQ/Luc1_2.fq.gz
trim_galore --paired --hardtrim5 100 ../FASTQ/Luc2_1.fq.gz         ../FASTQ/Luc2_2.fq.gz
trim_galore --paired --hardtrim5 100 ../FASTQ/Luc3_1.fq.gz         ../FASTQ/Luc3_2.fq.gz
trim_galore --paired --hardtrim5 100 ../FASTQ/Non1_1.fq.gz         ../FASTQ/Non1_2.fq.gz
trim_galore --paired --hardtrim5 100 ../FASTQ/Non2_1.fq.gz         ../FASTQ/Non2_2.fq.gz
trim_galore --paired --hardtrim5 100 ../FASTQ/Non3_1.fq.gz         ../FASTQ/Non3_2.fq.gz


# ### cut 9 based on fastqc
trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired Irf2_3utr1_1.100bp_5prime.fq.gz   Irf2_3utr1_2.100bp_5prime.fq.gz
trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired Irf2_3utr2_1.100bp_5prime.fq.gz   Irf2_3utr2_2.100bp_5prime.fq.gz
trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired Irf2_3utr3_1.100bp_5prime.fq.gz   Irf2_3utr3_2.100bp_5prime.fq.gz
trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired Irf2_cds1_1.100bp_5prime.fq.gz    Irf2_cds1_2.100bp_5prime.fq.gz
trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired Irf2_cds2_1.100bp_5prime.fq.gz    Irf2_cds2_2.100bp_5prime.fq.gz
trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired Irf2_cds3_1.100bp_5prime.fq.gz    Irf2_cds3_2.100bp_5prime.fq.gz
trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired Luc1_1.100bp_5prime.fq.gz         Luc1_2.100bp_5prime.fq.gz
trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired Luc2_1.100bp_5prime.fq.gz         Luc2_2.100bp_5prime.fq.gz
trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired Luc3_1.100bp_5prime.fq.gz         Luc3_2.100bp_5prime.fq.gz
trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired Non1_1.100bp_5prime.fq.gz         Non1_2.100bp_5prime.fq.gz
trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired Non2_1.100bp_5prime.fq.gz         Non2_2.100bp_5prime.fq.gz
trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired Non3_1.100bp_5prime.fq.gz         Non3_2.100bp_5prime.fq.gz


# rename output files
mv Irf2_3utr1_1.100bp_5prime_val_1.fq.gz   J774A1_shIrf2_3utr_rep1_1.fastq.gz      
mv Irf2_3utr2_1.100bp_5prime_val_1.fq.gz   J774A1_shIrf2_3utr_rep2_1.fastq.gz      
mv Irf2_3utr3_1.100bp_5prime_val_1.fq.gz   J774A1_shIrf2_3utr_rep3_1.fastq.gz      
mv Irf2_cds1_1.100bp_5prime_val_1.fq.gz    J774A1_shIrf2_cds_rep1_1.fastq.gz       
mv Irf2_cds2_1.100bp_5prime_val_1.fq.gz    J774A1_shIrf2_cds_rep2_1.fastq.gz       
mv Irf2_cds3_1.100bp_5prime_val_1.fq.gz    J774A1_shIrf2_cds_rep3_1.fastq.gz       
mv Luc1_1.100bp_5prime_val_1.fq.gz         J774A1_shLuc_rep1_1.fastq.gz            
mv Luc2_1.100bp_5prime_val_1.fq.gz         J774A1_shLuc_rep2_1.fastq.gz            
mv Luc3_1.100bp_5prime_val_1.fq.gz         J774A1_shLuc_rep3_1.fastq.gz            
mv Non1_1.100bp_5prime_val_1.fq.gz         J774A1_shNon1_1.fastq.gz                
mv Non2_1.100bp_5prime_val_1.fq.gz         J774A1_shNon2_1.fastq.gz                
mv Non3_1.100bp_5prime_val_1.fq.gz         J774A1_shNon3_1.fastq.gz                
mv Irf2_3utr1_2.100bp_5prime_val_2.fq.gz   J774A1_shIrf2_3utr_rep1_2.fastq.gz      
mv Irf2_3utr2_2.100bp_5prime_val_2.fq.gz   J774A1_shIrf2_3utr_rep2_2.fastq.gz      
mv Irf2_3utr3_2.100bp_5prime_val_2.fq.gz   J774A1_shIrf2_3utr_rep3_2.fastq.gz      
mv Irf2_cds1_2.100bp_5prime_val_2.fq.gz    J774A1_shIrf2_cds_rep1_2.fastq.gz       
mv Irf2_cds2_2.100bp_5prime_val_2.fq.gz    J774A1_shIrf2_cds_rep2_2.fastq.gz       
mv Irf2_cds3_2.100bp_5prime_val_2.fq.gz    J774A1_shIrf2_cds_rep3_2.fastq.gz       
mv Luc1_2.100bp_5prime_val_2.fq.gz         J774A1_shLuc_rep1_2.fastq.gz           
mv Luc2_2.100bp_5prime_val_2.fq.gz         J774A1_shLuc_rep2_2.fastq.gz           
mv Luc3_2.100bp_5prime_val_2.fq.gz         J774A1_shLuc_rep3_2.fastq.gz           
mv Non1_2.100bp_5prime_val_2.fq.gz         J774A1_shNon_rep1_2.fastq.gz           
mv Non2_2.100bp_5prime_val_2.fq.gz         J774A1_shNon_rep2_2.fastq.gz           
mv Non3_2.100bp_5prime_val_2.fq.gz         J774A1_shNon_rep3_2.fastq.gz           

