#############################################################
##############  README -- 2_Perimac_multiOmics ##############
#############################################################


##### 1_RNAseq: processing of bulk peritoneal RNA-seq dataset #####
		++ 1_Trimming_Mapping                 : scripts to quality trim FASTQ reads and map them to the mouse genome
		++ 2_DESeq2                           : script to perform differential gene expression analysis with DESeq2
		++ 3_Omics_Circos                     : scripts to plot DEGs on a circular genome plot
		++ 4_Network_Analyst                  : input data and parameters to run NetworkAnalyst analysis
		++ 5_Functional_Enrichment            : scripts to perform functional enrichment analysis on RNA-seq dataset
			- 1_GSEA                          : scripts to run GSEA enrichment analysis of female vs. male transcriptional aging
			- 2_X_chromosome                  : script to perform GSEA with curated X-related gene lists (and input gene lists)
			- 3_Top_Processes                 : scripts to test top candidate processes with curated gene lists by GSEA
				* 1_Phago                     : GSEA analysis for CRISPR/Cas9 hits in phagocytosis screens gene sets
				* 2_Metab                     : GSEA analysis for MSigDB Hallmarks metabolism-related gene sets
				* 3_Polarization              : GSE103958 analysis for M1/M2 DEGs, and GSEA analysis for polarization gene sets
			- 4_TF_analysis                   : scripts to identify candidate misregulated genes (and downstream candidate genes) 
				* 1_Candidates                : script to identify female-specific regulated TFs with aging using AnimalTFBD3.0 data
				* 2_Gene_expression_plotting  : scripts to plot gene expression from RNA-seq of genes of interest


##### 2_ATACseq: processing of bulk peritoneal ATAC-seq dataset #####
		++ 1_Trimming_Mapping_PeakCalling     : scripts to quality trim FASTQ reads, map them to the mouse genome and call meta ATAC-seq peaks
		++ 2_DESeq2                           : script to perform differential accessibility analysis
			- 1_Diffbind                      : script to run DiffBind to obtained normalized peak count matrices
			- 2_DESeq2                        : script to perform differential accessibility analysis with DEseq2
		++ 3_Omics_Circos                     : scripts to plot DAPs on a circular genome plot
		++ 5_Functional_Enrichment            : scripts to run GSEA enrichment analysis of female vs. male chromatin accessibility changes with aging
		++ 6_Motif_analysis                   : scripts to perform motif enrichment analysis in ATAC-seq peak sets


##### 3_CUT_and_RUN: processing of bulk peritoneal H3K4me3 CUT&RUN dataset #####
		++ 1_Trimming_Mapping_PeakCalling     : scripts to quality trim FASTQ reads, map them to the mouse genome and call meta H3K4me3 CUT&RUN peaks
		++ 2_DESeq2                           : script to perform differential H3K4me3 signal analysis
			- 1_Diffbind                      : script to run DiffBind to obtained normalized peak count matrices
			- 2_DESeq2                        : script to perform differential H3K4me3 signal analysis with DEseq2
		++ 3_Omics_Circos                     : scripts to plot differentially H3K4me3 modified peaks on a circular genome plot
		++ 5_Functional_Enrichment            : scripts to run GSEA enrichment analysis of female vs. male on differential H3K4me3 levels with aging
		++ 6_Motif_analysis                   : scripts to perform motif enrichment analysis in H3K4me3 CUT&RUN peak sets
