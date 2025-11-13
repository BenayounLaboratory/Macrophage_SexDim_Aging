setwd('/Volumes/BB_HQ_1/Immune_sex_dimorphism_Aging/Macrophages/Clean_Macrophage_Data/ATACseq/Functional_Enrichment_Analysis/GSEA_v2/')
options(stringsAsFactors = F)

source('Functions_GSEA_ATAC_v4.R')

# 2021-11-11
# Analyze PeriMac ATAC-seq based on GSEA

# 2022-04-08
# Use most signficant and 10kb threshold instead

# 2022-04-15
# rerun to output all (even not significant) to facilitate male/female aging comaparison

# 2023-11-01
# rerun using same pathways definitions as RNAseq

# 2024-07-11
# rerun GSEA with more conservative thresholds (min 20 genes, max 5000)
# to improve power and robustness

# load gene sets as prepped for RNAseq
load('../../../RNAseq/Functional_enrichment/GSEA_v2/2022-12-21_GeneSetCollections_for_Phenotest_GSEA.RData')


######################## A. Load DEseq2 results for analysis ######################## 
# load DEseq2 results, with annotated peaks
my.mph.atac.F_aging    <- read.csv("../../DEseq2/MACS2/2021-11-05_PeriMac_ATACseq_DESeq2_Analysis_AGING_FEMALE_all_genes_statistics_PeakAnnot.txt", header = T, sep = "\t")
my.mph.atac.M_aging    <- read.csv("../../DEseq2/MACS2/2021-11-05_PeriMac_ATACseq_DESeq2_Analysis_AGING_MALE_all_genes_statistics_PeakAnnot.txt"  , header = T, sep = "\t")

######################## B. Prepare GeneLists using DEseq2 t-statistic to rank genes and selecting 1 peak per gene ######################## 
# need to transfrom ATAC data, assign to closest gene if gene is within 10k
# if more than one peak within 10k, get the most significant one

glist.F_aging    <- prepare_atac_glist(my.mph.atac.F_aging    )
glist.M_aging    <- prepare_atac_glist(my.mph.atac.M_aging    )


######################## C. Gene Set Enrichment Analysis ######################## 
mph.aging.F.m2.reactome   <-  run_enrich(glist.F_aging, "ATAC_PeriMac_Aging_Females",  my.fdr = 0.05,  my.ontology =  Sym.m2.reactome       , my.ontology.name = "MSigDB_Reactome"          )  
mph.aging.F.go            <-  run_enrich(glist.F_aging, "ATAC_PeriMac_Aging_Females",  my.fdr = 0.05,  my.ontology =  Sym.ENS.GO.ALL        , my.ontology.name = "EMSEMBL_GO"               )  
mph.aging.F.m3.gtrd       <-  run_enrich(glist.F_aging, "ATAC_PeriMac_Aging_Females",  my.fdr = 0.05,  my.ontology =  Sym.m3.gtrd           , my.ontology.name = "MSigDB_GTRD"              )  
mph.aging.F.mh.all        <-  run_enrich(glist.F_aging, "ATAC_PeriMac_Aging_Females",  my.fdr = 0.05,  my.ontology =  Sym.mh.all            , my.ontology.name = "MSigDB_Hallmarks"         )  
mph.aging.F.TFlink        <-  run_enrich(glist.F_aging, "ATAC_PeriMac_Aging_Females",  my.fdr = 0.05,  my.ontology =  Sym.TFlink            , my.ontology.name = "TFlink_mouse"             )  

mph.aging.M.m2.reactome   <-  run_enrich(glist.M_aging, "ATAC_PeriMac_Aging_Males",  my.fdr = 0.05,  my.ontology =  Sym.m2.reactome       , my.ontology.name = "MSigDB_Reactome"          )  
mph.aging.M.go            <-  run_enrich(glist.M_aging, "ATAC_PeriMac_Aging_Males",  my.fdr = 0.05,  my.ontology =  Sym.ENS.GO.ALL        , my.ontology.name = "EMSEMBL_GO"               )  
mph.aging.M.m3.gtrd       <-  run_enrich(glist.M_aging, "ATAC_PeriMac_Aging_Males",  my.fdr = 0.05,  my.ontology =  Sym.m3.gtrd           , my.ontology.name = "MSigDB_GTRD"              )  
mph.aging.M.mh.all        <-  run_enrich(glist.M_aging, "ATAC_PeriMac_Aging_Males",  my.fdr = 0.05,  my.ontology =  Sym.mh.all            , my.ontology.name = "MSigDB_Hallmarks"         )  
mph.aging.M.TFlink        <-  run_enrich(glist.M_aging, "ATAC_PeriMac_Aging_Males",  my.fdr = 0.05,  my.ontology =  Sym.TFlink            , my.ontology.name = "TFlink_mouse"             )  


#######################
sink(file = paste(Sys.Date(),"_GSEA_ATAC_Sex_session_Info.txt", sep =""))
sessionInfo()
sink()

