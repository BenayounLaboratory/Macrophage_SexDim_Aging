~/Softwares/cellranger-6.1.2/bin/cellranger count --id JAX_perit_pool1 \
                 --transcriptome /home/benayoun/Softwares/cellranger-6.0.2/Reference/refdata-gex-mm10-2020-A \
                 --fastqs ~/Projects/10xGenomics/2021-12-21_JAX_Peritoneal_cohort/FASTQ \
                 --sample JAX_perit_pool1_SI_GA_B7_Trimmed  \
                 --expect-cells 14000 \
                 --localmem 128 \
                 --localcores 16

~/Softwares/cellranger-6.1.2/bin/cellranger count --id JAX_perit_pool2 \
                 --transcriptome /home/benayoun/Softwares/cellranger-6.0.2/Reference/refdata-gex-mm10-2020-A \
                 --fastqs ~/Projects/10xGenomics/2021-12-21_JAX_Peritoneal_cohort/FASTQ \
                 --sample JAX_perit_pool2_SI_GA_B8_Trimmed  \
                 --expect-cells 14000 \
                 --localmem 128 \
                 --localcores 16

~/Softwares/cellranger-6.1.2/bin/cellranger count --id JAX_perit_pool3 \
                 --transcriptome /home/benayoun/Softwares/cellranger-6.0.2/Reference/refdata-gex-mm10-2020-A \
                 --fastqs ~/Projects/10xGenomics/2021-12-21_JAX_Peritoneal_cohort/FASTQ \
                 --sample JAX_perit_pool3_SI_GA_B9_Trimmed  \
                 --expect-cells 14000 \
                 --localmem 128 \
                 --localcores 16
