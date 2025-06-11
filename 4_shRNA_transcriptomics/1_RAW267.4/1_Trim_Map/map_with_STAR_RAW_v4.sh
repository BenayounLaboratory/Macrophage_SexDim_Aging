# rerum counting for set 1 (Novogene mRNAseq is unstranded)

export GENIN="/Users/benayoun/Softwares/STAR-2.7.7a/genomes/MM10"

######################################### 
### set 1
#STAR --genomeDir $GENIN --readFilesIn ../FASTQ_HardTrim/shIrf2_3utr_rep1_RNAseq_1.fastq.gz     ../FASTQ_HardTrim/shIrf2_3utr_rep1_RNAseq_2.fastq.gz       --outFilterScoreMinOverLread 0.4 --outFilterMatchNminOverLread 0.4 --readFilesCommand gzcat --runThreadN 1 --outFilterMultimapNmax 50 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate  --alignEndsProtrude 15 ConcordantPair --outFileNamePrefix shIrf2_3utr_rep1_RNAseq_STAR_    
#STAR --genomeDir $GENIN --readFilesIn ../FASTQ_HardTrim/shIrf2_3utr_rep2_RNAseq_1.fastq.gz     ../FASTQ_HardTrim/shIrf2_3utr_rep2_RNAseq_2.fastq.gz       --outFilterScoreMinOverLread 0.4 --outFilterMatchNminOverLread 0.4 --readFilesCommand gzcat --runThreadN 1 --outFilterMultimapNmax 50 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate  --alignEndsProtrude 15 ConcordantPair --outFileNamePrefix shIrf2_3utr_rep2_RNAseq_STAR_    
#STAR --genomeDir $GENIN --readFilesIn ../FASTQ_HardTrim/shIrf2_3utr_rep3_RNAseq_1.fastq.gz     ../FASTQ_HardTrim/shIrf2_3utr_rep3_RNAseq_2.fastq.gz       --outFilterScoreMinOverLread 0.4 --outFilterMatchNminOverLread 0.4 --readFilesCommand gzcat --runThreadN 1 --outFilterMultimapNmax 50 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate  --alignEndsProtrude 15 ConcordantPair --outFileNamePrefix shIrf2_3utr_rep3_RNAseq_STAR_    
#STAR --genomeDir $GENIN --readFilesIn ../FASTQ_HardTrim/shIrf2_cds_rep1_RNAseq_1.fastq.gz      ../FASTQ_HardTrim/shIrf2_cds_rep1_RNAseq_2.fastq.gz        --outFilterScoreMinOverLread 0.4 --outFilterMatchNminOverLread 0.4 --readFilesCommand gzcat --runThreadN 1 --outFilterMultimapNmax 50 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate  --alignEndsProtrude 15 ConcordantPair --outFileNamePrefix shIrf2_cds_rep1_RNAseq_STAR_     
#STAR --genomeDir $GENIN --readFilesIn ../FASTQ_HardTrim/shIrf2_cds_rep2_RNAseq_1.fastq.gz      ../FASTQ_HardTrim/shIrf2_cds_rep2_RNAseq_2.fastq.gz        --outFilterScoreMinOverLread 0.4 --outFilterMatchNminOverLread 0.4 --readFilesCommand gzcat --runThreadN 1 --outFilterMultimapNmax 50 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate  --alignEndsProtrude 15 ConcordantPair --outFileNamePrefix shIrf2_cds_rep2_RNAseq_STAR_     
#STAR --genomeDir $GENIN --readFilesIn ../FASTQ_HardTrim/shIrf2_cds_rep3_RNAseq_1.fastq.gz      ../FASTQ_HardTrim/shIrf2_cds_rep3_RNAseq_2.fastq.gz        --outFilterScoreMinOverLread 0.4 --outFilterMatchNminOverLread 0.4 --readFilesCommand gzcat --runThreadN 1 --outFilterMultimapNmax 50 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate  --alignEndsProtrude 15 ConcordantPair --outFileNamePrefix shIrf2_cds_rep3_RNAseq_STAR_     
#STAR --genomeDir $GENIN --readFilesIn ../FASTQ_HardTrim/shLuc_rep1_RNAseq_1.fastq.gz           ../FASTQ_HardTrim/shLuc_rep1_RNAseq_2.fastq.gz             --outFilterScoreMinOverLread 0.4 --outFilterMatchNminOverLread 0.4 --readFilesCommand gzcat --runThreadN 1 --outFilterMultimapNmax 50 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate  --alignEndsProtrude 15 ConcordantPair --outFileNamePrefix shLuc_rep1_RNAseq_STAR_          
#STAR --genomeDir $GENIN --readFilesIn ../FASTQ_HardTrim/shLuc_rep2_RNAseq_1.fastq.gz           ../FASTQ_HardTrim/shLuc_rep2_RNAseq_2.fastq.gz             --outFilterScoreMinOverLread 0.4 --outFilterMatchNminOverLread 0.4 --readFilesCommand gzcat --runThreadN 1 --outFilterMultimapNmax 50 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate  --alignEndsProtrude 15 ConcordantPair --outFileNamePrefix shLuc_rep2_RNAseq_STAR_          
#STAR --genomeDir $GENIN --readFilesIn ../FASTQ_HardTrim/shLuc_rep3_RNAseq_1.fastq.gz           ../FASTQ_HardTrim/shLuc_rep3_RNAseq_2.fastq.gz             --outFilterScoreMinOverLread 0.4 --outFilterMatchNminOverLread 0.4 --readFilesCommand gzcat --runThreadN 1 --outFilterMultimapNmax 50 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate  --alignEndsProtrude 15 ConcordantPair --outFileNamePrefix shLuc_rep3_RNAseq_STAR_          
#STAR --genomeDir $GENIN --readFilesIn ../FASTQ_HardTrim/shNon_rep1_RNAseq_1.fastq.gz           ../FASTQ_HardTrim/shNon_rep1_RNAseq_2.fastq.gz             --outFilterScoreMinOverLread 0.4 --outFilterMatchNminOverLread 0.4 --readFilesCommand gzcat --runThreadN 1 --outFilterMultimapNmax 50 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate  --alignEndsProtrude 15 ConcordantPair --outFileNamePrefix shNon_rep1_RNAseq_STAR_          
#STAR --genomeDir $GENIN --readFilesIn ../FASTQ_HardTrim/shNon_rep2_RNAseq_1.fastq.gz           ../FASTQ_HardTrim/shNon_rep2_RNAseq_2.fastq.gz             --outFilterScoreMinOverLread 0.4 --outFilterMatchNminOverLread 0.4 --readFilesCommand gzcat --runThreadN 1 --outFilterMultimapNmax 50 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate  --alignEndsProtrude 15 ConcordantPair --outFileNamePrefix shNon_rep2_RNAseq_STAR_          
#STAR --genomeDir $GENIN --readFilesIn ../FASTQ_HardTrim/shNon_rep3_RNAseq_1.fastq.gz           ../FASTQ_HardTrim/shNon_rep3_RNAseq_2.fastq.gz             --outFilterScoreMinOverLread 0.4 --outFilterMatchNminOverLread 0.4 --readFilesCommand gzcat --runThreadN 1 --outFilterMultimapNmax 50 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate  --alignEndsProtrude 15 ConcordantPair --outFileNamePrefix shNon_rep3_RNAseq_STAR_          
#STAR --genomeDir $GENIN --readFilesIn ../FASTQ_HardTrim/shTbl1xr1_3UTR_rep1_RNAseq_1.fastq.gz  ../FASTQ_HardTrim/shTbl1xr1_3UTR_rep1_RNAseq_2.fastq.gz    --outFilterScoreMinOverLread 0.4 --outFilterMatchNminOverLread 0.4 --readFilesCommand gzcat --runThreadN 1 --outFilterMultimapNmax 50 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate  --alignEndsProtrude 15 ConcordantPair --outFileNamePrefix shTbl1xr1_3UTR_rep1_RNAseq_STAR_ 
#STAR --genomeDir $GENIN --readFilesIn ../FASTQ_HardTrim/shTbl1xr1_3UTR_rep2_RNAseq_1.fastq.gz  ../FASTQ_HardTrim/shTbl1xr1_3UTR_rep2_RNAseq_2.fastq.gz    --outFilterScoreMinOverLread 0.4 --outFilterMatchNminOverLread 0.4 --readFilesCommand gzcat --runThreadN 1 --outFilterMultimapNmax 50 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate  --alignEndsProtrude 15 ConcordantPair --outFileNamePrefix shTbl1xr1_3UTR_rep2_RNAseq_STAR_ 
#STAR --genomeDir $GENIN --readFilesIn ../FASTQ_HardTrim/shTbl1xr1_3UTR_rep3_RNAseq_1.fastq.gz  ../FASTQ_HardTrim/shTbl1xr1_3UTR_rep3_RNAseq_2.fastq.gz    --outFilterScoreMinOverLread 0.4 --outFilterMatchNminOverLread 0.4 --readFilesCommand gzcat --runThreadN 1 --outFilterMultimapNmax 50 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate  --alignEndsProtrude 15 ConcordantPair --outFileNamePrefix shTbl1xr1_3UTR_rep3_RNAseq_STAR_ 
#STAR --genomeDir $GENIN --readFilesIn ../FASTQ_HardTrim/shTbl1xr1_cds_rep1_RNAseq_1.fastq.gz   ../FASTQ_HardTrim/shTbl1xr1_cds_rep1_RNAseq_2.fastq.gz     --outFilterScoreMinOverLread 0.4 --outFilterMatchNminOverLread 0.4 --readFilesCommand gzcat --runThreadN 1 --outFilterMultimapNmax 50 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate  --alignEndsProtrude 15 ConcordantPair --outFileNamePrefix shTbl1xr1_cds_rep1_RNAseq_STAR_  
#STAR --genomeDir $GENIN --readFilesIn ../FASTQ_HardTrim/shTbl1xr1_cds_rep3_RNAseq_1.fastq.gz   ../FASTQ_HardTrim/shTbl1xr1_cds_rep3_RNAseq_2.fastq.gz     --outFilterScoreMinOverLread 0.4 --outFilterMatchNminOverLread 0.4 --readFilesCommand gzcat --runThreadN 1 --outFilterMultimapNmax 50 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate  --alignEndsProtrude 15 ConcordantPair --outFileNamePrefix shTbl1xr1_cds_rep3_RNAseq_STAR_  
          
# ### set 2
# STAR --genomeDir $GENIN --readFilesIn ../FASTQ_HardTrim/shTal1_3UTR_rep3_RNAseq_SET2_1.fastq.gz ../FASTQ_HardTrim/shTal1_3UTR_rep3_RNAseq_SET2_2.fastq.gz  --outFilterScoreMinOverLread 0.4 --outFilterMatchNminOverLread 0.4 --readFilesCommand gzcat --runThreadN 1 --outFilterMultimapNmax 50 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate  --alignEndsProtrude 15 ConcordantPair --outFileNamePrefix shTal1_3UTR_rep3_RNAseq_SET2_STAR_    
# STAR --genomeDir $GENIN --readFilesIn ../FASTQ_HardTrim/shTal1_3UTR_rep2_RNAseq_SET2_1.fastq.gz ../FASTQ_HardTrim/shTal1_3UTR_rep2_RNAseq_SET2_2.fastq.gz  --outFilterScoreMinOverLread 0.4 --outFilterMatchNminOverLread 0.4 --readFilesCommand gzcat --runThreadN 1 --outFilterMultimapNmax 50 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate  --alignEndsProtrude 15 ConcordantPair --outFileNamePrefix shTal1_3UTR_rep2_RNAseq_SET2_STAR_    
# STAR --genomeDir $GENIN --readFilesIn ../FASTQ_HardTrim/shTal1_3UTR_rep1_RNAseq_SET2_1.fastq.gz ../FASTQ_HardTrim/shTal1_3UTR_rep1_RNAseq_SET2_2.fastq.gz  --outFilterScoreMinOverLread 0.4 --outFilterMatchNminOverLread 0.4 --readFilesCommand gzcat --runThreadN 1 --outFilterMultimapNmax 50 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate  --alignEndsProtrude 15 ConcordantPair --outFileNamePrefix shTal1_3UTR_rep1_RNAseq_SET2_STAR_    
# STAR --genomeDir $GENIN --readFilesIn ../FASTQ_HardTrim/shTal1_CDS_rep3_RNAseq_SET2_1.fastq.gz  ../FASTQ_HardTrim/shTal1_CDS_rep3_RNAseq_SET2_2.fastq.gz   --outFilterScoreMinOverLread 0.4 --outFilterMatchNminOverLread 0.4 --readFilesCommand gzcat --runThreadN 1 --outFilterMultimapNmax 50 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate  --alignEndsProtrude 15 ConcordantPair --outFileNamePrefix shTal1_CDS_rep3_RNAseq_SET2_STAR_    
# STAR --genomeDir $GENIN --readFilesIn ../FASTQ_HardTrim/shTal1_CDS_rep2_RNAseq_SET2_1.fastq.gz  ../FASTQ_HardTrim/shTal1_CDS_rep2_RNAseq_SET2_2.fastq.gz   --outFilterScoreMinOverLread 0.4 --outFilterMatchNminOverLread 0.4 --readFilesCommand gzcat --runThreadN 1 --outFilterMultimapNmax 50 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate  --alignEndsProtrude 15 ConcordantPair --outFileNamePrefix shTal1_CDS_rep2_RNAseq_SET2_STAR_    
# STAR --genomeDir $GENIN --readFilesIn ../FASTQ_HardTrim/shTal1_CDS_rep1_RNAseq_SET2_1.fastq.gz  ../FASTQ_HardTrim/shTal1_CDS_rep1_RNAseq_SET2_2.fastq.gz   --outFilterScoreMinOverLread 0.4 --outFilterMatchNminOverLread 0.4 --readFilesCommand gzcat --runThreadN 1 --outFilterMultimapNmax 50 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate  --alignEndsProtrude 15 ConcordantPair --outFileNamePrefix shTal1_CDS_rep1_RNAseq_SET2_STAR_    
# STAR --genomeDir $GENIN --readFilesIn ../FASTQ_HardTrim/shNon_rep6_RNAseq_SET2_1.fastq.gz       ../FASTQ_HardTrim/shNon_rep6_RNAseq_SET2_2.fastq.gz        --outFilterScoreMinOverLread 0.4 --outFilterMatchNminOverLread 0.4 --readFilesCommand gzcat --runThreadN 1 --outFilterMultimapNmax 50 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate  --alignEndsProtrude 15 ConcordantPair --outFileNamePrefix shNon_rep6_RNAseq_SET2_2_STAR_    
# STAR --genomeDir $GENIN --readFilesIn ../FASTQ_HardTrim/shNon_rep5_RNAseq_SET2_1.fastq.gz       ../FASTQ_HardTrim/shNon_rep5_RNAseq_SET2_2.fastq.gz        --outFilterScoreMinOverLread 0.4 --outFilterMatchNminOverLread 0.4 --readFilesCommand gzcat --runThreadN 1 --outFilterMultimapNmax 50 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate  --alignEndsProtrude 15 ConcordantPair --outFileNamePrefix shNon_rep5_RNAseq_SET2_2_STAR_    
# STAR --genomeDir $GENIN --readFilesIn ../FASTQ_HardTrim/shNon_rep4_RNAseq_SET2_1.fastq.gz       ../FASTQ_HardTrim/shNon_rep4_RNAseq_SET2_2.fastq.gz        --outFilterScoreMinOverLread 0.4 --outFilterMatchNminOverLread 0.4 --readFilesCommand gzcat --runThreadN 1 --outFilterMultimapNmax 50 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate  --alignEndsProtrude 15 ConcordantPair --outFileNamePrefix shNon_rep4_RNAseq_SET2_2_STAR_    
# STAR --genomeDir $GENIN --readFilesIn ../FASTQ_HardTrim/shLuc_rep6_RNAseq_SET2_1.fastq.gz       ../FASTQ_HardTrim/shLuc_rep6_RNAseq_SET2_2.fastq.gz        --outFilterScoreMinOverLread 0.4 --outFilterMatchNminOverLread 0.4 --readFilesCommand gzcat --runThreadN 1 --outFilterMultimapNmax 50 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate  --alignEndsProtrude 15 ConcordantPair --outFileNamePrefix shLuc_rep6_RNAseq_SET2_2_STAR_    
# STAR --genomeDir $GENIN --readFilesIn ../FASTQ_HardTrim/shLuc_rep5_RNAseq_SET2_1.fastq.gz       ../FASTQ_HardTrim/shLuc_rep5_RNAseq_SET2_2.fastq.gz        --outFilterScoreMinOverLread 0.4 --outFilterMatchNminOverLread 0.4 --readFilesCommand gzcat --runThreadN 1 --outFilterMultimapNmax 50 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate  --alignEndsProtrude 15 ConcordantPair --outFileNamePrefix shLuc_rep5_RNAseq_SET2_2_STAR_    
# STAR --genomeDir $GENIN --readFilesIn ../FASTQ_HardTrim/shLuc_rep4_RNAseq_SET2_1.fastq.gz       ../FASTQ_HardTrim/shLuc_rep4_RNAseq_SET2_2.fastq.gz        --outFilterScoreMinOverLread 0.4 --outFilterMatchNminOverLread 0.4 --readFilesCommand gzcat --runThreadN 1 --outFilterMultimapNmax 50 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate  --alignEndsProtrude 15 ConcordantPair --outFileNamePrefix shLuc_rep4_RNAseq_SET2_2_STAR_    

### set 3
#STAR --genomeDir $GENIN --readFilesIn ../FASTQ_HardTrim/shNon_rep9_RNAseq_SET3_1.fastq.gz ../FASTQ_HardTrim/shNon_rep9_RNAseq_SET3_2.fastq.gz              --outFilterScoreMinOverLread 0.4 --outFilterMatchNminOverLread 0.4 --readFilesCommand gzcat --runThreadN 1 --outFilterMultimapNmax 50 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate  --alignEndsProtrude 15 ConcordantPair --outFileNamePrefix shNon_rep9_RNAseq_SET3_STAR_ 
#STAR --genomeDir $GENIN --readFilesIn ../FASTQ_HardTrim/shNon_rep8_RNAseq_SET3_1.fastq.gz ../FASTQ_HardTrim/shNon_rep8_RNAseq_SET3_2.fastq.gz              --outFilterScoreMinOverLread 0.4 --outFilterMatchNminOverLread 0.4 --readFilesCommand gzcat --runThreadN 1 --outFilterMultimapNmax 50 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate  --alignEndsProtrude 15 ConcordantPair --outFileNamePrefix shNon_rep8_RNAseq_SET3_STAR_ 
#STAR --genomeDir $GENIN --readFilesIn ../FASTQ_HardTrim/shNon_rep7_RNAseq_SET3_1.fastq.gz ../FASTQ_HardTrim/shNon_rep7_RNAseq_SET3_2.fastq.gz              --outFilterScoreMinOverLread 0.4 --outFilterMatchNminOverLread 0.4 --readFilesCommand gzcat --runThreadN 1 --outFilterMultimapNmax 50 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate  --alignEndsProtrude 15 ConcordantPair --outFileNamePrefix shNon_rep7_RNAseq_SET3_STAR_ 
#STAR --genomeDir $GENIN --readFilesIn ../FASTQ_HardTrim/shMeis13utr_rep3_RNAseq_SET3_1.fastq.gz ../FASTQ_HardTrim/shMeis13utr_rep3_RNAseq_SET3_2.fastq.gz  --outFilterScoreMinOverLread 0.4 --outFilterMatchNminOverLread 0.4 --readFilesCommand gzcat --runThreadN 1 --outFilterMultimapNmax 50 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate  --alignEndsProtrude 15 ConcordantPair --outFileNamePrefix shMeis13utr_rep3_RNAseq_SET3_STAR_    
#STAR --genomeDir $GENIN --readFilesIn ../FASTQ_HardTrim/shMeis13utr_rep2_RNAseq_SET3_1.fastq.gz ../FASTQ_HardTrim/shMeis13utr_rep2_RNAseq_SET3_2.fastq.gz  --outFilterScoreMinOverLread 0.4 --outFilterMatchNminOverLread 0.4 --readFilesCommand gzcat --runThreadN 1 --outFilterMultimapNmax 50 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate  --alignEndsProtrude 15 ConcordantPair --outFileNamePrefix shMeis13utr_rep2_RNAseq_SET3_STAR_    
#STAR --genomeDir $GENIN --readFilesIn ../FASTQ_HardTrim/shMeis13utr_rep1_RNAseq_SET3_1.fastq.gz ../FASTQ_HardTrim/shMeis13utr_rep1_RNAseq_SET3_2.fastq.gz  --outFilterScoreMinOverLread 0.4 --outFilterMatchNminOverLread 0.4 --readFilesCommand gzcat --runThreadN 1 --outFilterMultimapNmax 50 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate  --alignEndsProtrude 15 ConcordantPair --outFileNamePrefix shMeis13utr_rep1_RNAseq_SET3_STAR_    
#STAR --genomeDir $GENIN --readFilesIn ../FASTQ_HardTrim/shMeis1cds_rep3_RNAseq_SET3_1.fastq.gz ../FASTQ_HardTrim/shMeis1cds_rep3_RNAseq_SET3_2.fastq.gz    --outFilterScoreMinOverLread 0.4 --outFilterMatchNminOverLread 0.4 --readFilesCommand gzcat --runThreadN 1 --outFilterMultimapNmax 50 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate  --alignEndsProtrude 15 ConcordantPair --outFileNamePrefix shMeis1cds_rep3_RNAseq_SET3_STAR_  
#STAR --genomeDir $GENIN --readFilesIn ../FASTQ_HardTrim/shMeis1cds_rep2_RNAseq_SET3_1.fastq.gz ../FASTQ_HardTrim/shMeis1cds_rep2_RNAseq_SET3_2.fastq.gz    --outFilterScoreMinOverLread 0.4 --outFilterMatchNminOverLread 0.4 --readFilesCommand gzcat --runThreadN 1 --outFilterMultimapNmax 50 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate  --alignEndsProtrude 15 ConcordantPair --outFileNamePrefix shMeis1cds_rep2_RNAseq_SET3_STAR_  
#STAR --genomeDir $GENIN --readFilesIn ../FASTQ_HardTrim/shMeis1cds_rep1_RNAseq_SET3_1.fastq.gz ../FASTQ_HardTrim/shMeis1cds_rep1_RNAseq_SET3_2.fastq.gz    --outFilterScoreMinOverLread 0.4 --outFilterMatchNminOverLread 0.4 --readFilesCommand gzcat --runThreadN 1 --outFilterMultimapNmax 50 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate  --alignEndsProtrude 15 ConcordantPair --outFileNamePrefix shMeis1cds_rep1_RNAseq_SET3_STAR_  
#STAR --genomeDir $GENIN --readFilesIn ../FASTQ_HardTrim/shLuc_rep9_RNAseq_SET3_1.fastq.gz ../FASTQ_HardTrim/shLuc_rep9_RNAseq_SET3_2.fastq.gz              --outFilterScoreMinOverLread 0.4 --outFilterMatchNminOverLread 0.4 --readFilesCommand gzcat --runThreadN 1 --outFilterMultimapNmax 50 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate  --alignEndsProtrude 15 ConcordantPair --outFileNamePrefix shLuc_rep9_RNAseq_SET3_1_STAR_ 
#STAR --genomeDir $GENIN --readFilesIn ../FASTQ_HardTrim/shLuc_rep8_RNAseq_SET3_1.fastq.gz ../FASTQ_HardTrim/shLuc_rep8_RNAseq_SET3_2.fastq.gz              --outFilterScoreMinOverLread 0.4 --outFilterMatchNminOverLread 0.4 --readFilesCommand gzcat --runThreadN 1 --outFilterMultimapNmax 50 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate  --alignEndsProtrude 15 ConcordantPair --outFileNamePrefix shLuc_rep8_RNAseq_SET3_1_STAR_ 
#STAR --genomeDir $GENIN --readFilesIn ../FASTQ_HardTrim/shLuc_rep7_RNAseq_SET3_1.fastq.gz ../FASTQ_HardTrim/shLuc_rep7_RNAseq_SET3_2.fastq.gz              --outFilterScoreMinOverLread 0.4 --outFilterMatchNminOverLread 0.4 --readFilesCommand gzcat --runThreadN 1 --outFilterMultimapNmax 50 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate  --alignEndsProtrude 15 ConcordantPair --outFileNamePrefix shLuc_rep7_RNAseq_SET3_1_STAR_ 

### set 4
STAR --genomeDir $GENIN --readFilesIn ../FASTQ_HardTrim/shLuc_rep10_RNAseq_SET4_1.fastq.gz ../FASTQ_HardTrim/shLuc_rep10_RNAseq_SET4_2.fastq.gz                  --outFilterScoreMinOverLread 0.4 --outFilterMatchNminOverLread 0.4 --readFilesCommand gzcat --runThreadN 1 --outFilterMultimapNmax 50 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate  --alignEndsProtrude 15 ConcordantPair --outFileNamePrefix shLuc_rep10_RNAseq_SET4_STAR
STAR --genomeDir $GENIN --readFilesIn ../FASTQ_HardTrim/shLuc_rep11_RNAseq_SET4_1.fastq.gz ../FASTQ_HardTrim/shLuc_rep11_RNAseq_SET4_2.fastq.gz                  --outFilterScoreMinOverLread 0.4 --outFilterMatchNminOverLread 0.4 --readFilesCommand gzcat --runThreadN 1 --outFilterMultimapNmax 50 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate  --alignEndsProtrude 15 ConcordantPair --outFileNamePrefix shLuc_rep11_RNAseq_SET4_STAR
STAR --genomeDir $GENIN --readFilesIn ../FASTQ_HardTrim/shLuc_rep12_RNAseq_SET4_1.fastq.gz ../FASTQ_HardTrim/shLuc_rep12_RNAseq_SET4_2.fastq.gz                  --outFilterScoreMinOverLread 0.4 --outFilterMatchNminOverLread 0.4 --readFilesCommand gzcat --runThreadN 1 --outFilterMultimapNmax 50 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate  --alignEndsProtrude 15 ConcordantPair --outFileNamePrefix shLuc_rep12_RNAseq_SET4_STAR
STAR --genomeDir $GENIN --readFilesIn ../FASTQ_HardTrim/shMef2c_69_rep1_RNAseq_SET4_1.fastq.gz ../FASTQ_HardTrim/shMef2c_69_rep1_RNAseq_SET4_2.fastq.gz          --outFilterScoreMinOverLread 0.4 --outFilterMatchNminOverLread 0.4 --readFilesCommand gzcat --runThreadN 1 --outFilterMultimapNmax 50 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate  --alignEndsProtrude 15 ConcordantPair --outFileNamePrefix shMef2c_69_rep1_RNAseq_SET4_STAR
STAR --genomeDir $GENIN --readFilesIn ../FASTQ_HardTrim/shMef2c_69_rep2_RNAseq_SET4_1.fastq.gz ../FASTQ_HardTrim/shMef2c_69_rep2_RNAseq_SET4_2.fastq.gz          --outFilterScoreMinOverLread 0.4 --outFilterMatchNminOverLread 0.4 --readFilesCommand gzcat --runThreadN 1 --outFilterMultimapNmax 50 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate  --alignEndsProtrude 15 ConcordantPair --outFileNamePrefix shMef2c_69_rep2_RNAseq_SET4_STAR
STAR --genomeDir $GENIN --readFilesIn ../FASTQ_HardTrim/shMef2c_69_rep3_RNAseq_SET4_1.fastq.gz ../FASTQ_HardTrim/shMef2c_69_rep3_RNAseq_SET4_2.fastq.gz          --outFilterScoreMinOverLread 0.4 --outFilterMatchNminOverLread 0.4 --readFilesCommand gzcat --runThreadN 1 --outFilterMultimapNmax 50 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate  --alignEndsProtrude 15 ConcordantPair --outFileNamePrefix shMef2c_69_rep3_RNAseq_SET4_STAR
STAR --genomeDir $GENIN --readFilesIn ../FASTQ_HardTrim/shMef2c_72_rep1_RNAseq_SET4_1.fastq.gz ../FASTQ_HardTrim/shMef2c_72_rep1_RNAseq_SET4_2.fastq.gz          --outFilterScoreMinOverLread 0.4 --outFilterMatchNminOverLread 0.4 --readFilesCommand gzcat --runThreadN 1 --outFilterMultimapNmax 50 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate  --alignEndsProtrude 15 ConcordantPair --outFileNamePrefix shMef2c_72_rep1_RNAseq_SET4_STAR
STAR --genomeDir $GENIN --readFilesIn ../FASTQ_HardTrim/shMef2c_72_rep2_RNAseq_SET4_1.fastq.gz ../FASTQ_HardTrim/shMef2c_72_rep2_RNAseq_SET4_2.fastq.gz          --outFilterScoreMinOverLread 0.4 --outFilterMatchNminOverLread 0.4 --readFilesCommand gzcat --runThreadN 1 --outFilterMultimapNmax 50 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate  --alignEndsProtrude 15 ConcordantPair --outFileNamePrefix shMef2c_72_rep2_RNAseq_SET4_STAR
STAR --genomeDir $GENIN --readFilesIn ../FASTQ_HardTrim/shMef2c_72_rep3_RNAseq_SET4_1.fastq.gz ../FASTQ_HardTrim/shMef2c_72_rep3_RNAseq_SET4_2.fastq.gz          --outFilterScoreMinOverLread 0.4 --outFilterMatchNminOverLread 0.4 --readFilesCommand gzcat --runThreadN 1 --outFilterMultimapNmax 50 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate  --alignEndsProtrude 15 ConcordantPair --outFileNamePrefix shMef2c_72_rep3_RNAseq_SET4_STAR
STAR --genomeDir $GENIN --readFilesIn ../FASTQ_HardTrim/shNon_rep10_RNAseq_SET4_1.fastq.gz ../FASTQ_HardTrim/shNon_rep10_RNAseq_SET4_2.fastq.gz                  --outFilterScoreMinOverLread 0.4 --outFilterMatchNminOverLread 0.4 --readFilesCommand gzcat --runThreadN 1 --outFilterMultimapNmax 50 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate  --alignEndsProtrude 15 ConcordantPair --outFileNamePrefix shNon_rep10_RNAseq_SET4_STAR
STAR --genomeDir $GENIN --readFilesIn ../FASTQ_HardTrim/shNon_rep11_RNAseq_SET4_1.fastq.gz ../FASTQ_HardTrim/shNon_rep11_RNAseq_SET4_2.fastq.gz                  --outFilterScoreMinOverLread 0.4 --outFilterMatchNminOverLread 0.4 --readFilesCommand gzcat --runThreadN 1 --outFilterMultimapNmax 50 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate  --alignEndsProtrude 15 ConcordantPair --outFileNamePrefix shNon_rep11_RNAseq_SET4_STAR
STAR --genomeDir $GENIN --readFilesIn ../FASTQ_HardTrim/shNon_rep12_RNAseq_SET4_1.fastq.gz ../FASTQ_HardTrim/shNon_rep12_RNAseq_SET4_2.fastq.gz                  --outFilterScoreMinOverLread 0.4 --outFilterMatchNminOverLread 0.4 --readFilesCommand gzcat --runThreadN 1 --outFilterMultimapNmax 50 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate  --alignEndsProtrude 15 ConcordantPair --outFileNamePrefix shNon_rep12_RNAseq_SET4_STAR
                      
# ######################################### 
# one file per KD TF
### set 1
featureCounts -t exon -D 1500 -p --primary -T 3 -s 0 -a /Users/benayoun/Softwares/Genomes/mm10_genes.gtf -o 2024-08-15_RAW_shIRF2_CLEAN_counts_unstranded.txt \
     shNon_rep1_RNAseq_STAR_Aligned.sortedByCoord.out.bam          \
     shNon_rep2_RNAseq_STAR_Aligned.sortedByCoord.out.bam          \
     shNon_rep3_RNAseq_STAR_Aligned.sortedByCoord.out.bam          \
     shLuc_rep1_RNAseq_STAR_Aligned.sortedByCoord.out.bam          \
     shLuc_rep2_RNAseq_STAR_Aligned.sortedByCoord.out.bam          \
     shLuc_rep3_RNAseq_STAR_Aligned.sortedByCoord.out.bam          \
     shIrf2_3utr_rep1_RNAseq_STAR_Aligned.sortedByCoord.out.bam    \
     shIrf2_3utr_rep2_RNAseq_STAR_Aligned.sortedByCoord.out.bam    \
     shIrf2_3utr_rep3_RNAseq_STAR_Aligned.sortedByCoord.out.bam    \
     shIrf2_cds_rep1_RNAseq_STAR_Aligned.sortedByCoord.out.bam     \
     shIrf2_cds_rep2_RNAseq_STAR_Aligned.sortedByCoord.out.bam     \
     shIrf2_cds_rep3_RNAseq_STAR_Aligned.sortedByCoord.out.bam     

featureCounts -t exon -D 1500 -p --primary -T 3 -s 0 -a /Users/benayoun/Softwares/Genomes/mm10_genes.gtf -o 2024-08-15_RAW_shTBL1XR1_CLEAN_counts_unstranded.txt \
     shNon_rep1_RNAseq_STAR_Aligned.sortedByCoord.out.bam          \
     shNon_rep2_RNAseq_STAR_Aligned.sortedByCoord.out.bam          \
     shNon_rep3_RNAseq_STAR_Aligned.sortedByCoord.out.bam          \
     shLuc_rep1_RNAseq_STAR_Aligned.sortedByCoord.out.bam          \
     shLuc_rep2_RNAseq_STAR_Aligned.sortedByCoord.out.bam          \
     shLuc_rep3_RNAseq_STAR_Aligned.sortedByCoord.out.bam          \
     shTbl1xr1_3UTR_rep1_RNAseq_STAR_Aligned.sortedByCoord.out.bam \
     shTbl1xr1_3UTR_rep2_RNAseq_STAR_Aligned.sortedByCoord.out.bam \
     shTbl1xr1_3UTR_rep3_RNAseq_STAR_Aligned.sortedByCoord.out.bam \
     shTbl1xr1_cds_rep1_RNAseq_STAR_Aligned.sortedByCoord.out.bam  \
     shTbl1xr1_cds_rep3_RNAseq_STAR_Aligned.sortedByCoord.out.bam
     
### set 2
featureCounts -t exon -D 1500 -p --primary -T 3 -s 0 -a /Users/benayoun/Softwares/Genomes/mm10_genes.gtf -o 2024-08-15_RAW_shTAL1_CLEAN_counts_unstranded.txt \
     shNon_rep4_RNAseq_SET2_2_STAR_Aligned.sortedByCoord.out.bam       \
     shNon_rep5_RNAseq_SET2_2_STAR_Aligned.sortedByCoord.out.bam       \
     shNon_rep6_RNAseq_SET2_2_STAR_Aligned.sortedByCoord.out.bam       \
     shLuc_rep4_RNAseq_SET2_2_STAR_Aligned.sortedByCoord.out.bam       \
     shLuc_rep5_RNAseq_SET2_2_STAR_Aligned.sortedByCoord.out.bam       \
     shLuc_rep6_RNAseq_SET2_2_STAR_Aligned.sortedByCoord.out.bam       \
     shTal1_3UTR_rep1_RNAseq_SET2_STAR_Aligned.sortedByCoord.out.bam   \
     shTal1_3UTR_rep2_RNAseq_SET2_STAR_Aligned.sortedByCoord.out.bam   \
     shTal1_3UTR_rep3_RNAseq_SET2_STAR_Aligned.sortedByCoord.out.bam   \
     shTal1_CDS_rep1_RNAseq_SET2_STAR_Aligned.sortedByCoord.out.bam    \
     shTal1_CDS_rep2_RNAseq_SET2_STAR_Aligned.sortedByCoord.out.bam    \
     shTal1_CDS_rep3_RNAseq_SET2_STAR_Aligned.sortedByCoord.out.bam   
     
          
### set 3
featureCounts -t exon -D 1500 -p --primary -T 3 -s 0 -a /Users/benayoun/Softwares/Genomes/mm10_genes.gtf -o 22024-08-15_RAW_shMeis1_CLEAN_counts_unstranded.txt \
	shNon_rep9_RNAseq_SET3_STAR_Aligned.sortedByCoord.out.bam       \
	shNon_rep8_RNAseq_SET3_STAR_Aligned.sortedByCoord.out.bam       \
	shNon_rep7_RNAseq_SET3_STAR_Aligned.sortedByCoord.out.bam       \
	shLuc_rep9_RNAseq_SET3_1_STAR_Aligned.sortedByCoord.out.bam     \
	shLuc_rep8_RNAseq_SET3_1_STAR_Aligned.sortedByCoord.out.bam     \
	shLuc_rep7_RNAseq_SET3_1_STAR_Aligned.sortedByCoord.out.bam     \
	shMeis13utr_rep3_RNAseq_SET3_STAR_Aligned.sortedByCoord.out.bam \
	shMeis13utr_rep2_RNAseq_SET3_STAR_Aligned.sortedByCoord.out.bam \
	shMeis13utr_rep1_RNAseq_SET3_STAR_Aligned.sortedByCoord.out.bam \
	shMeis1cds_rep3_RNAseq_SET3_STAR_Aligned.sortedByCoord.out.bam  \
	shMeis1cds_rep2_RNAseq_SET3_STAR_Aligned.sortedByCoord.out.bam  \
	shMeis1cds_rep1_RNAseq_SET3_STAR_Aligned.sortedByCoord.out.bam         

### set 4
featureCounts -t exon -D 1500 -p --primary -T 3 -s 0 -a /Users/benayoun/Softwares/Genomes/mm10_genes.gtf -o 22024-08-15_RAW_shMef2c_CLEAN_counts_unstranded.txt \
	shLuc_rep10_RNAseq_SET4_STARAligned.sortedByCoord.out.bam      \
	shLuc_rep11_RNAseq_SET4_STARAligned.sortedByCoord.out.bam      \
	shLuc_rep12_RNAseq_SET4_STARAligned.sortedByCoord.out.bam      \
	shMef2c_69_rep1_RNAseq_SET4_STARAligned.sortedByCoord.out.bam  \
	shMef2c_69_rep2_RNAseq_SET4_STARAligned.sortedByCoord.out.bam  \
	shMef2c_69_rep3_RNAseq_SET4_STARAligned.sortedByCoord.out.bam  \
	shMef2c_72_rep1_RNAseq_SET4_STARAligned.sortedByCoord.out.bam  \
	shMef2c_72_rep2_RNAseq_SET4_STARAligned.sortedByCoord.out.bam  \
	shMef2c_72_rep3_RNAseq_SET4_STARAligned.sortedByCoord.out.bam  \
	shNon_rep10_RNAseq_SET4_STARAligned.sortedByCoord.out.bam      \
	shNon_rep11_RNAseq_SET4_STARAligned.sortedByCoord.out.bam      \
	shNon_rep12_RNAseq_SET4_STARAligned.sortedByCoord.out.bam   














