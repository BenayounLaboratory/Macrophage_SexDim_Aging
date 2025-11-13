library(DESeq2)
library(phenoTest)
library(qusage)   

#################################################################
### DEBUG
# my.matrix          <-   mph.aging.F.geneList
# data.name          <-  "PeriMac_Aging_Females"
# my.fdr             <-   0.05
# my.ontology        <-   Sym.c2.cp.kegg
# my.ontology.name   <-  "MSigDB_KEGG"

run_enrich <- function(my.matrix, data.name, my.fdr = 0.05, my.ontology, my.ontology.name) {
  
  # set seed to stabilize output
  set.seed(123456789)
  
  # run phenotest GSEA
  gsea.data <- gsea( x         =  my.matrix    , 
                     gsets     =  my.ontology , 
                     mc.cores  =  2            , 
                     logScale  =  FALSE        , 
                     B         =  10000         , 
                     minGenes  =  20            , 
                     maxGenes  =  5000         )
  
  my.summary <- data.frame(summary(gsea.data))
  my.sig.path.num <- sum(my.summary$fdr < my.fdr )
  
  # write results to file
  my.outfile <- paste(Sys.Date(),data.name, my.ontology.name, "FDR", 100*my.fdr, "Phenotest_GSEA_Analysis_table", my.sig.path.num,"significant_OUTPUT_ALL.txt", sep = "_")
  write.table(my.summary, file = my.outfile, quote = F, sep = "\t")
  
  return(my.summary)
}

