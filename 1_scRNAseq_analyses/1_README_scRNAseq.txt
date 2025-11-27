############################################################
##############  README -- 1_scRNAseq_analyses ##############
############################################################


##### 1_Individual_processing: preprocessing of each independent scRNAseq dataset #####
		++ 1_NIA_v2: 10xGenomics v2 dataset on cells from NIA mice
			- CellRanger     : cellranger scripts to obtain gene count matrices
			- Seurat         : scripts for data processing in Seurat (doublet removal, cell annotation)
			- AUGUR_CellProp : scripts to perform Augur analysis as well as cell proportion analysis with scProportionTest
		
		++ 2_JAX_v3: 10xGenomics v3.1 dataset on cells from JAX mice; includes HTO processing
			- CellRanger     : cellranger scripts to obtain gene count matrices and perform HTO quantification
			- Seurat         : scripts for data processing in Seurat (HTO demultiplexing, doublet removal, cell annotation)
			- AUGUR_CellProp : scripts to perform Augur analysis as well as cell proportion analysis with scProportionTest

		++ 3_NIA_v3: 10xGenomics v3.1 dataset on cells from NIA mice; includes HTO processing
			- CellRanger     : cellranger scripts to obtain gene count matrices and perform HTO quantification
			- Seurat         : scripts for data processing in Seurat (HTO demultiplexing, doublet removal, cell annotation)
			- AUGUR_CellProp : scripts to perform Augur analysis as well as cell proportion analysis with scProportionTest


##### 2_Joint_analyses: joint analysis of the 3 scRNAseq datasets #####
		++ AUGUR_CellProp : scripts to perform Augur analysis as well as cell proportion analysis with scProportionTest on merged data
		++ Mph_Only       : scripts to plot gene expression from scRNAseq for top candidate genes
		++ SCENIC         : scripts to perform TF activity prediction using SCENIC