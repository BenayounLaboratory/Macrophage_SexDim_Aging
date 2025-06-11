#cellranger count --id 4m_Female_perit_pool_1 \
#                 --transcriptome /home/benayoun/Softwares/cellranger-6.0.2/Reference/refdata-gex-mm10-2020-A \
#                 --fastqs /home/benayoun/Projects/10xGenomics/2018-05-08/180508_Peritoneal_cells/ \
#                 --sample 4m_Female \
#                 --expect-cells 4000 \
#                 --localmem 64 \
#                 --localcores 12

#cellranger count --id 20m_Female_perit_pool_1 \
#                 --transcriptome /home/benayoun/Softwares/cellranger-6.0.2/Reference/refdata-gex-mm10-2020-A \
#                 --fastqs /home/benayoun/Projects/10xGenomics/2018-05-08/180508_Peritoneal_cells/ \
#                 --sample 20m_Female \
#                 --expect-cells 4000 \
#                 --localmem 64 \
#                 --localcores 12

#cellranger count --id 4m_Male_perit_pool_1 \
#                 --transcriptome /home/benayoun/Softwares/cellranger-6.0.2/Reference/refdata-gex-mm10-2020-A \
#                 --fastqs /home/benayoun/Projects/10xGenomics/2018-05-08/180508_Peritoneal_cells/ \
#                 --sample 4m_Male \
#                 --expect-cells 4000 \
#                 --localmem 64 \
#                 --localcores 12

#cellranger count --id 20m_Male_perit_pool_1 \
#                 --transcriptome /home/benayoun/Softwares/cellranger-6.0.2/Reference/refdata-gex-mm10-2020-A \
#                 --fastqs /home/benayoun/Projects/10xGenomics/2018-05-08/180508_Peritoneal_cells/ \
#                 --sample 20m_Male \
#                 --expect-cells 4000 \
#                 --localmem 64 \
#                 --localcores 12
                 
                 
#cellranger count --id 4m_Female_perit_pool_2 \
#                 --transcriptome /home/benayoun/Softwares/cellranger-6.0.2/Reference/refdata-gex-mm10-2020-A \
#                 --fastqs /home/benayoun/Projects/10xGenomics/2018-09-13/180911_Peritoneal_and_liver_cells/ \
#                 --sample 4m_Female_peritoneal \
#                 --expect-cells 3000 \
#                 --localmem 64 \
#                 --localcores 12

#cellranger count --id 20m_Female_perit_pool_2 \
#                 --transcriptome /home/benayoun/Softwares/cellranger-6.0.2/Reference/refdata-gex-mm10-2020-A \
#                 --fastqs /home/benayoun/Projects/10xGenomics/2018-09-13/180911_Peritoneal_and_liver_cells/ \
#                 --sample 20m_Female_peritoneal \
#                 --expect-cells 3000 \
#                 --localmem 64 \
#                 --localcores 12

cellranger count --id 4m_Male_perit_pool_2 \
                 --transcriptome /home/benayoun/Softwares/cellranger-6.0.2/Reference/refdata-gex-mm10-2020-A \
                 --fastqs /home/benayoun/Projects/10xGenomics/2018-09-13/180911_Peritoneal_and_liver_cells/ \
                 --sample 4m_Male_peritoneal \
                 --expect-cells 3000 \
                 --localmem 64 \
                 --localcores 12

#cellranger count --id 20m_Male_perit_pool_2 \
#                 --transcriptome /home/benayoun/Softwares/cellranger-6.0.2/Reference/refdata-gex-mm10-2020-A \
#                 --fastqs /home/benayoun/Projects/10xGenomics/2018-09-13/180911_Peritoneal_and_liver_cells/ \
#                 --sample 20m_Male_peritoneal \
#                 --expect-cells 3000 \
#                 --localmem 64 \
#                 --localcores 12
                 
cellranger aggr --id=Peritoneal_scRNAseq_reanalysis \
                 --csv=peritoneal_cells_ALL_experiments.csv \
                 --normalize=none

