###########################################################
##############    README -- 5_Irf2_ChIPseq   ##############
###########################################################


##### 1_Map                   : scripts to map ChIP-seq FASTQ reads to the mouse genome #####

##### 2_PeakCalling           : processing of bulk peritoneal RNA-seq dataset #####
		++ 1_MACS2            : script to perform peak calling in each ChIP-seq replicate with MACS2
		++ 2_MSPC             : script to obtain robust peaks using MSPC

##### 3_Motif                 : script to perform motif enrichment analysis in Irf2 ChIP-seq MSPC robust peaks #####

##### 4_Functional_Enrichment : processing of bulk peritoneal RNA-seq dataset #####
		++ 1_GREAT            : scripts for parsing and plotting GREAT functional enrichment results associated to Irf2 ChIP-seq MSPC robust peaks
		++ 2_ChIPseeeker      : scripts for running and plotting ChIPseeker functional enrichment results associated to Irf2 ChIP-seq MSPC robust peaks

##### 5_ChIP_Comparison_RNA   : scripts to compare Irf2 direct target genes with genes regulated during female macrophage aging #####
