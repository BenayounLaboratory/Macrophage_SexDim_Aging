setwd('/Volumes/BB_HQ_1//Immune_sex_dimorphism_Aging/Macrophages/Clean_Macrophage_Data/RNAseq/Functional_enrichment/GSEA_v2/')
options(stringsAsFactors = F)


source('Functions_GSEA_v3.R')

load('2022-12-21_GeneSetCollections_for_Phenotest_GSEA.RData')

######################## A. Load DEseq2 results for analysis ######################## 
# load DEseq2 results
load('../../DEseq2_analysis/DEseq2/2022-03-30_Peritoneal_Macrophages_RNA_seq_results.RData')
my.mph1.RNAseq.process
# [1] "Sex"             "Aging"           "Female_Aging"    "Male_Aging"      "FemaleMaleAging"

my.mph.aging.F           <- data.frame(my.mph1.RNAseq.process$Female_Aging)
my.mph.aging.F$GeneName  <- rownames(my.mph.aging.F )

my.mph.aging.M           <- data.frame(my.mph1.RNAseq.process$Male_Aging)
my.mph.aging.M$GeneName  <- rownames(my.mph.aging.M )

######################## B. Prepare GeneLists using DEseq2 t-statistic to rank genes ######################## 
mph.aging.F.geneList         <- my.mph.aging.F$stat
names(mph.aging.F.geneList)  <- my.mph.aging.F$GeneName
mph.aging.F.geneList         <- sort(mph.aging.F.geneList , decreasing = TRUE)

mph.aging.M.geneList         <- my.mph.aging.M$stat
names(mph.aging.M.geneList)  <- my.mph.aging.M$GeneName
mph.aging.M.geneList         <- sort(mph.aging.M.geneList , decreasing = TRUE)

######################## C. Gene Set Enrichment Analysis ######################## 
mph.aging.F.m2.reactome   <-  run_enrich(mph.aging.F.geneList, "PeriMac_Aging_Females",  my.fdr = 0.05,  my.ontology =  Sym.m2.reactome       , my.ontology.name = "MSigDB_Reactome"          )
mph.aging.F.m5.go         <-  run_enrich(mph.aging.F.geneList, "PeriMac_Aging_Females",  my.fdr = 0.05,  my.ontology =  Sym.ENS.GO.ALL        , my.ontology.name = "EMSEMBL_GO"               )
mph.aging.F.m3.gtrd       <-  run_enrich(mph.aging.F.geneList, "PeriMac_Aging_Females",  my.fdr = 0.05,  my.ontology =  Sym.m3.gtrd           , my.ontology.name = "MSigDB_GTRD"              )


mph.aging.M.m2.reactome   <-  run_enrich(mph.aging.M.geneList, "PeriMac_Aging_Males"  ,  my.fdr = 0.05,  my.ontology =  Sym.m2.reactome       , my.ontology.name = "MSigDB_Reactome"          )
mph.aging.M.m5.go         <-  run_enrich(mph.aging.M.geneList, "PeriMac_Aging_Males"  ,  my.fdr = 0.05,  my.ontology =  Sym.ENS.GO.ALL        , my.ontology.name = "EMSEMBL_GO"               )
mph.aging.M.m3.gtrd       <-  run_enrich(mph.aging.M.geneList, "PeriMac_Aging_Males"  ,  my.fdr = 0.05,  my.ontology =  Sym.m3.gtrd           , my.ontology.name = "MSigDB_GTRD"              )


#######################
sink(file = paste0(Sys.Date(),"_GSEA_Sex_session_Info_STRINGENT.txt"))
sessionInfo()
sink()

