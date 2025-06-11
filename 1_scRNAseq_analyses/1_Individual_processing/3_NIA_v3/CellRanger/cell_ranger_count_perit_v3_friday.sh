cellranger count --id perit_pool_1_CKDL230020266-1A_H755CDSX7 \
                 --transcriptome /home/benayoun/Softwares/cellranger-6.0.2/Reference/refdata-gex-mm10-2020-A \
                 --fastqs /home/benayoun/Projects/2023-07-18_NIA_peritv3_cohort/FASTQ \
                 --sample NIA_pool1 \
                 --expect-cells 30000 \
                 --localmem 32 \
                 --localcores 6

