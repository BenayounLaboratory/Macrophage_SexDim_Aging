setwd('/Volumes/BB_HQ_1//Immune_sex_dimorphism_Aging/Macrophages/Clean_Macrophage_Data/RNAseq/Functional_enrichment/GSEA_v2/')
options(stringsAsFactors = F)

# 2022-12-20
# prep MSigDB mouse gmts for fast loading for GSEA with phenotest
# TF targets from TFlink

library(DESeq2)
library(phenoTest)
library(qusage)


#########
Sym.m2.reactome        <- read.gmt("/Volumes/BB_HQ_1/PATHWAY_ANNOT/MSigDB/m2.cp.reactome.v2022.1.Mm.symbols.gmt"                                              )
Sym.m3.gtrd            <- read.gmt("/Volumes/BB_HQ_1/PATHWAY_ANNOT/MSigDB/m3.gtrd.v2022.1.Mm.symbols.gmt"                                                     )
Sym.mh.all             <- read.gmt("/Volumes/BB_HQ_1/PATHWAY_ANNOT/MSigDB/mh.all.v2022.1.Mm.symbols.gmt"                                                      )

Sym.TFlink             <- read.gmt("/Volumes/BB_HQ_1/PATHWAY_ANNOT/TF_targets/TFLink_Mus_musculus_interactions_All_GMT_proteinName_v1.0.gmt"                                                       )

# Use parsed clean GO from ENSEMBL (instead of MsigDB)
# Sym.m5.go              <- read.gmt("/Volumes/BB_HQ_1/PATHWAY_ANNOT/MSigDB/m5.go.v2022.1.Mm.symbols.gmt"                                                     )
Sym.ENS.GO.ALL         <- read.gmt("/Volumes/BB_HQ_1/PATHWAY_ANNOT/ENSEMBL/2022-12-21_mouse_Ens108_GO_ALL.gmt"                                                 ) 

save(Sym.m2.reactome   ,
     Sym.ENS.GO.ALL    ,
     Sym.m3.gtrd       ,
     Sym.mh.all        ,
     Sym.TFlink        ,
     file = paste0(Sys.Date(),"_GeneSetCollections_for_Phenotest_GSEA.RData"))
