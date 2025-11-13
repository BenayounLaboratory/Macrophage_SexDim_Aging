library(DESeq2)
library(phenoTest)
library(qusage)   

#################################################################
# need to transfrom H3K4me3 height data, assign to closest gene since TSS mark
# only keep if within 10kb of annotated gene

# my.peak.annot <- my.mph.H3K4me3_height.F_aging
# my.abs.dist = 1e4

prepare_H3K4me3_height_glist <- function(my.peak.annot, my.abs.dist = 1e4) {
  
  # ignore peaks with no gene assignment
  my.peak.annot <- my.peak.annot[my.peak.annot$Gene.Name != "",]
  
  # keep only annotations within 10k of assigned gene
  my.peak.annot.2 <- my.peak.annot[abs(my.peak.annot$Distance.to.TSS) < my.abs.dist, ]
  
  # get list of genes
  my.genes <- unique(my.peak.annot.2$Gene.Name)
  
  # create a summary data.frame with only one entry per gene
  my.H3K4me3_height.sum <- data.frame(matrix(0,length(my.genes),4))
  rownames(my.H3K4me3_height.sum) <- my.genes
  colnames(my.H3K4me3_height.sum) <- c("GeneName","logFC","stat","padj")
  my.H3K4me3_height.sum$GeneName  <- my.genes
  
  for (i in 1:length(my.genes)) {
    
    # get H3K4me3_height peak corresponding to gene with closest TSS
    my.H3K4me3_height.peaks <- my.peak.annot.2[my.peak.annot.2$Gene.Name %in% my.genes[i],]
    my.closest <- which.min(abs(my.H3K4me3_height.peaks$Distance.to.TSS))
    
    # collect data
    my.H3K4me3_height.sum$logFC[i] <- my.H3K4me3_height.peaks[my.closest,]$log2FoldChange
    my.H3K4me3_height.sum$stat[i]  <- my.H3K4me3_height.peaks[my.closest,]$stat
    my.H3K4me3_height.sum$padj[i]  <- my.H3K4me3_height.peaks[my.closest,]$padj
    
  }
  
  H3K4me3_height.geneList             <- my.H3K4me3_height.sum$stat
  names(H3K4me3_height.geneList)      <- my.H3K4me3_height.sum$GeneName
  H3K4me3_height.geneList             <- sort(H3K4me3_height.geneList , decreasing = TRUE)
  
  return(H3K4me3_height.geneList)
}



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

