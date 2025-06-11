#### gzcat ../FASTQ/Irf2cds1_1.fq.gz           |  fastx_trimmer -f 9 -l 100 -z -Q33 -i - -o shIrf2_cds_rep1_RNAseq_SET2_1.fastq.gz      

### use trimgalore to recapitulate the behavior of fastx_trimmer (doesn't run on newer macos)
## cut to 100 to help mapping with STAR
#trim_galore --paired --hardtrim5 100 ../FASTQ/SET2/Luc4_1.fq.gz ../FASTQ/SET2/Luc4_2.fq.gz
#trim_galore --paired --hardtrim5 100 ../FASTQ/SET2/Luc5_1.fq.gz ../FASTQ/SET2/Luc5_2.fq.gz
#trim_galore --paired --hardtrim5 100 ../FASTQ/SET2/Luc6_1.fq.gz ../FASTQ/SET2/Luc6_2.fq.gz
#trim_galore --paired --hardtrim5 100 ../FASTQ/SET2/Non4_1.fq.gz ../FASTQ/SET2/Non4_2.fq.gz
#trim_galore --paired --hardtrim5 100 ../FASTQ/SET2/Non5_1.fq.gz ../FASTQ/SET2/Non5_2.fq.gz
#trim_galore --paired --hardtrim5 100 ../FASTQ/SET2/Non6_1.fq.gz ../FASTQ/SET2/Non6_2.fq.gz
#trim_galore --paired --hardtrim5 100 ../FASTQ/SET2/Tal1CDS1_1.fq.gz ../FASTQ/SET2/Tal1CDS1_2.fq.gz
#trim_galore --paired --hardtrim5 100 ../FASTQ/SET2/Tal1CDS2_1.fq.gz ../FASTQ/SET2/Tal1CDS2_2.fq.gz
#trim_galore --paired --hardtrim5 100 ../FASTQ/SET2/Tal1CDS3_1.fq.gz ../FASTQ/SET2/Tal1CDS3_2.fq.gz
#trim_galore --paired --hardtrim5 100 ../FASTQ/SET2/TAL13UTR1_1.fq.gz ../FASTQ/SET2/TAL13UTR1_2.fq.gz
#trim_galore --paired --hardtrim5 100 ../FASTQ/SET2/TAL13UTR2_1.fq.gz ../FASTQ/SET2/TAL13UTR2_2.fq.gz
#trim_galore --paired --hardtrim5 100 ../FASTQ/SET2/TAL13UTR3_1.fq.gz ../FASTQ/SET2/TAL13UTR3_2.fq.gz
#
### cut 9 based on fastqc
#trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired Luc4_1.100bp_5prime.fq.gz Luc4_2.100bp_5prime.fq.gz
#trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired Luc5_1.100bp_5prime.fq.gz Luc5_2.100bp_5prime.fq.gz
#trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired Luc6_1.100bp_5prime.fq.gz Luc6_2.100bp_5prime.fq.gz
#trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired Non4_1.100bp_5prime.fq.gz Non4_2.100bp_5prime.fq.gz
#trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired Non5_1.100bp_5prime.fq.gz Non5_2.100bp_5prime.fq.gz
#trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired Non6_1.100bp_5prime.fq.gz Non6_2.100bp_5prime.fq.gz
#trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired Tal1CDS1_1.100bp_5prime.fq.gz Tal1CDS1_2.100bp_5prime.fq.gz
#trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired Tal1CDS2_1.100bp_5prime.fq.gz Tal1CDS2_2.100bp_5prime.fq.gz
#trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired Tal1CDS3_1.100bp_5prime.fq.gz Tal1CDS3_2.100bp_5prime.fq.gz
#trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired TAL13UTR1_1.100bp_5prime.fq.gz TAL13UTR1_2.100bp_5prime.fq.gz
#trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired TAL13UTR2_1.100bp_5prime.fq.gz TAL13UTR2_2.100bp_5prime.fq.gz
#trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired TAL13UTR3_1.100bp_5prime.fq.gz TAL13UTR3_2.100bp_5prime.fq.gz

# rename output files

#mv Luc4_1.100bp_5prime_val_1.fq.gz        shLuc_rep4_RNAseq_SET2_1.100bp_5prime_val_1.fastq.gz
#mv Luc4_2.100bp_5prime_val_2.fq.gz        shLuc_rep4_RNAseq_SET2_2.100bp_5prime_val_2.fastq.gz
#mv Luc5_1.100bp_5prime_val_1.fq.gz        shLuc_rep5_RNAseq_SET2_1.100bp_5prime_val_1.fastq.gz
#mv Luc5_2.100bp_5prime_val_2.fq.gz        shLuc_rep5_RNAseq_SET2_2.100bp_5prime_val_2.fastq.gz
#mv Luc6_1.100bp_5prime_val_1.fq.gz        shLuc_rep6_RNAseq_SET2_1.100bp_5prime_val_1.fastq.gz
#mv Luc6_2.100bp_5prime_val_2.fq.gz        shLuc_rep6_RNAseq_SET2_2.100bp_5prime_val_2.fastq.gz
#mv Non4_1.100bp_5prime_val_1.fq.gz        shNon_rep4_RNAseq_SET2_1.100bp_5prime_val_1.fastq.gz
#mv Non4_2.100bp_5prime_val_2.fq.gz        shNon_rep4_RNAseq_SET2_2.100bp_5prime_val_2.fastq.gz
#mv Non5_1.100bp_5prime_val_1.fq.gz        shNon_rep5_RNAseq_SET2_1.100bp_5prime_val_1.fastq.gz
#mv Non5_2.100bp_5prime_val_2.fq.gz        shNon_rep5_RNAseq_SET2_2.100bp_5prime_val_2.fastq.gz
#mv Non6_1.100bp_5prime_val_1.fq.gz        shNon_rep6_RNAseq_SET2_1.100bp_5prime_val_1.fastq.gz
#mv Non6_2.100bp_5prime_val_2.fq.gz        shNon_rep6_RNAseq_SET2_2.100bp_5prime_val_2.fastq.gz
#mv Tal1CDS1_1.100bp_5prime_val_1.fq.gz    shTal1_CDS_rep1_RNAseq_SET2_1.100bp_5prime_val_1.fastq.gz
#mv Tal1CDS1_2.100bp_5prime_val_2.fq.gz    shTal1_CDS_rep1_RNAseq_SET2_2.100bp_5prime_val_2.fastq.gz
#mv Tal1CDS2_1.100bp_5prime_val_1.fq.gz    shTal1_CDS_rep2_RNAseq_SET2_1.100bp_5prime_val_1.fastq.gz
#mv Tal1CDS2_2.100bp_5prime_val_2.fq.gz    shTal1_CDS_rep2_RNAseq_SET2_2.100bp_5prime_val_2.fastq.gz
#mv Tal1CDS3_1.100bp_5prime_val_1.fq.gz    shTal1_CDS_rep3_RNAseq_SET2_1.100bp_5prime_val_1.fastq.gz
#mv Tal1CDS3_2.100bp_5prime_val_2.fq.gz    shTal1_CDS_rep3_RNAseq_SET2_2.100bp_5prime_val_2.fastq.gz
#mv TAL13UTR1_1.100bp_5prime_val_1.fq.gz   shTal1_3UTR_rep1_RNAseq_SET2_1.100bp_5prime_val_1.fastq.gz
#mv TAL13UTR1_2.100bp_5prime_val_2.fq.gz   shTal1_3UTR_rep1_RNAseq_SET2_2.100bp_5prime_val_2.fastq.gz
#mv TAL13UTR2_1.100bp_5prime_val_1.fq.gz   shTal1_3UTR_rep2_RNAseq_SET2_1.100bp_5prime_val_1.fastq.gz
#mv TAL13UTR2_2.100bp_5prime_val_2.fq.gz   shTal1_3UTR_rep2_RNAseq_SET2_2.100bp_5prime_val_2.fastq.gz
#mv TAL13UTR3_1.100bp_5prime_val_1.fq.gz   shTal1_3UTR_rep3_RNAseq_SET2_1.100bp_5prime_val_1.fastq.gz
#mv TAL13UTR3_2.100bp_5prime_val_2.fq.gz   shTal1_3UTR_rep3_RNAseq_SET2_2.100bp_5prime_val_2.fastq.gz

mv shTal1_3UTR_rep3_RNAseq_SET2_1.100bp_5prime_val_1.fastq.gz   shTal1_3UTR_rep3_RNAseq_SET2_1.fastq.gz
mv shTal1_3UTR_rep3_RNAseq_SET2_2.100bp_5prime_val_2.fastq.gz   shTal1_3UTR_rep3_RNAseq_SET2_2.fastq.gz
mv shTal1_3UTR_rep2_RNAseq_SET2_1.100bp_5prime_val_1.fastq.gz   shTal1_3UTR_rep2_RNAseq_SET2_1.fastq.gz
mv shTal1_3UTR_rep2_RNAseq_SET2_2.100bp_5prime_val_2.fastq.gz   shTal1_3UTR_rep2_RNAseq_SET2_2.fastq.gz
mv shTal1_3UTR_rep1_RNAseq_SET2_2.100bp_5prime_val_2.fastq.gz   shTal1_3UTR_rep1_RNAseq_SET2_2.fastq.gz
mv shTal1_3UTR_rep1_RNAseq_SET2_1.100bp_5prime_val_1.fastq.gz   shTal1_3UTR_rep1_RNAseq_SET2_1.fastq.gz
mv shTal1_CDS_rep3_RNAseq_SET2_2.100bp_5prime_val_2.fastq.gz    shTal1_CDS_rep3_RNAseq_SET2_2.fastq.gz
mv shTal1_CDS_rep3_RNAseq_SET2_1.100bp_5prime_val_1.fastq.gz    shTal1_CDS_rep3_RNAseq_SET2_1.fastq.gz
mv shTal1_CDS_rep2_RNAseq_SET2_1.100bp_5prime_val_1.fastq.gz    shTal1_CDS_rep2_RNAseq_SET2_1.fastq.gz
mv shTal1_CDS_rep2_RNAseq_SET2_2.100bp_5prime_val_2.fastq.gz    shTal1_CDS_rep2_RNAseq_SET2_2.fastq.gz
mv shTal1_CDS_rep1_RNAseq_SET2_1.100bp_5prime_val_1.fastq.gz    shTal1_CDS_rep1_RNAseq_SET2_1.fastq.gz
mv shTal1_CDS_rep1_RNAseq_SET2_2.100bp_5prime_val_2.fastq.gz    shTal1_CDS_rep1_RNAseq_SET2_2.fastq.gz
mv shNon_rep6_RNAseq_SET2_1.100bp_5prime_val_1.fastq.gz         shNon_rep6_RNAseq_SET2_1.fastq.gz
mv shNon_rep6_RNAseq_SET2_2.100bp_5prime_val_2.fastq.gz         shNon_rep6_RNAseq_SET2_2.fastq.gz
mv shNon_rep5_RNAseq_SET2_1.100bp_5prime_val_1.fastq.gz         shNon_rep5_RNAseq_SET2_1.fastq.gz
mv shNon_rep5_RNAseq_SET2_2.100bp_5prime_val_2.fastq.gz         shNon_rep5_RNAseq_SET2_2.fastq.gz
mv shNon_rep4_RNAseq_SET2_2.100bp_5prime_val_2.fastq.gz         shNon_rep4_RNAseq_SET2_2.fastq.gz
mv shNon_rep4_RNAseq_SET2_1.100bp_5prime_val_1.fastq.gz         shNon_rep4_RNAseq_SET2_1.fastq.gz
mv shLuc_rep6_RNAseq_SET2_1.100bp_5prime_val_1.fastq.gz         shLuc_rep6_RNAseq_SET2_1.fastq.gz
mv shLuc_rep6_RNAseq_SET2_2.100bp_5prime_val_2.fastq.gz         shLuc_rep6_RNAseq_SET2_2.fastq.gz
mv shLuc_rep5_RNAseq_SET2_1.100bp_5prime_val_1.fastq.gz         shLuc_rep5_RNAseq_SET2_1.fastq.gz
mv shLuc_rep5_RNAseq_SET2_2.100bp_5prime_val_2.fastq.gz         shLuc_rep5_RNAseq_SET2_2.fastq.gz
mv shLuc_rep4_RNAseq_SET2_1.100bp_5prime_val_1.fastq.gz         shLuc_rep4_RNAseq_SET2_1.fastq.gz
mv shLuc_rep4_RNAseq_SET2_2.100bp_5prime_val_2.fastq.gz         shLuc_rep4_RNAseq_SET2_2.fastq.gz
