##################################################################
##############    README -- 5_Irf2_ChIPseq_CUTRUN   ##############
##################################################################

##############################################
## &&&&& ##    1_BMDM_ChIP-seq     ## &&&&& ##
##############################################

##### 1_Map                   : scripts to map ChIP-seq FASTQ reads to the mouse genome #####

##### 2_PeakCalling           : peak calling analysis #####
		++ 1_HOMER            : scripts to perform peak calling in each ChIP-seq replicate with HOMER
		++ 2_MSPC             : script to obtain robust peaks using MSPC

##### 3_Motif                 : script to perform motif enrichment analysis in Irf2 ChIP-seq MSPC robust peaks #####

##### 4_Functional_Enrichment : functional enrichment analysis #####
		++ 1_GREAT            : scripts for parsing and plotting GREAT functional enrichment results associated to Irf2 ChIP-seq MSPC robust peaks
		++ 2_ChIPseeeker      : scripts for running and plotting ChIPseeker functional enrichment results associated to Irf2 ChIP-seq MSPC robust peaks

##### 5_Bedgraph              : scripts to generate Irf2 ChIP-seq bedgraphs for visualization #####



###############################################
## &&&&& ##    2_PeriMac_CUTRUN     ## &&&&& ##
###############################################

##### 1_Map                   : scripts to map CUT&RUN FASTQ reads to the mouse genome #####

##### 2_PeakCalling           : peak calling analysis #####
		++ 1_HOMER            : script to perform peak calling in each CUT&RUN sample with HOMER
		++ 2_MSPC             : script to obtain robust peaks using MSPC

##### 3_Motif                 : script to perform motif enrichment analysis in Irf2 CUT&RUN MSPC robust peaks #####

##### 4_Functional_Enrichment : unctional enrichment analysis #####
		++ 1_GREAT            : scripts for parsing and plotting GREAT functional enrichment results associated to Irf2 CUT&RUN MSPC robust peaks
		++ 2_ChIPseeeker      : scripts for running and plotting ChIPseeker functional enrichment results associated to Irf2 CUT&RUN MSPC robust peaks

##### 5_Bedgraph              : scripts to generate Irf2 CUT&RUN bedgraphs, including spike-in scaling, for visualization #####

##### 6_Diffbind              : scripts to perform differential Irf2 CUT&RUN analysis #####


###############################################
## &&&&& ##    3_Joint_Analysis     ## &&&&& ##
###############################################

##### 1_GREAT_parsing          : scripts to parse GREAT enrichments for ChIP and CUT&RUN for plotting #####

##### 2_OmicsCircos            : script to plot positions of ChIP and CUT&RUN peaks on mm10 mouse genome #####

##### 3_ChIPseeker             : script to annotate ChIP and CUT&RUN peaks relative to genes #####

##### 4_Comparison_RNA         : scripts to compare ChIP and CUT&RUN peaks, associated genes and regulated genes in RNAseq #####
