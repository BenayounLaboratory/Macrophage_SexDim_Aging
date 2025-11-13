library(DESeq2)
library(phenoTest)
library(qusage)   

#################################################################
# need to transfrom ATAC data, assign to closest gene if gene is within 10kb
# if more than one peak within 10k, get the most significant one

# my.peak.annot <- my.mph.atac.F_aging
# my.abs.dist = 1e4

prepare_atac_glist <- function(my.peak.annot, my.abs.dist = 1e4) {
  
  # ignore peaks with no gene assignment
  my.peak.annot <- my.peak.annot[my.peak.annot$Gene.Name != "",]
  
  # keep only annotations within 10k of assigned gene
  my.peak.annot.2 <- my.peak.annot[abs(my.peak.annot$Distance.to.TSS) < my.abs.dist, ]
  
  # get list of genes
  my.genes <- unique(my.peak.annot.2$Gene.Name)
  
  # create a summary data.frame with only one entry per gene
  my.atac.sum <- data.frame(matrix(0,length(my.genes),4))
  rownames(my.atac.sum) <- my.genes
  colnames(my.atac.sum) <- c("GeneName","logFC","stat","padj")
  my.atac.sum$GeneName  <- my.genes
  
  for (i in 1:length(my.genes)) {
    
    # get atac peak corresponding to gene with minimize p adjust
    my.atac.peaks <- my.peak.annot.2[my.peak.annot.2$Gene.Name %in% my.genes[i],]
    my.min.p      <- which(my.atac.peaks$padj == min(my.atac.peaks$padj))
    
    # collect data
    # if more than one peak at min p, get largest absolute FC
    
    if (length(my.min.p) == 1) {
      my.atac.sum$logFC[i] <- my.atac.peaks[my.min.p,]$log2FoldChange
      my.atac.sum$stat[i]  <- my.atac.peaks[my.min.p,]$stat
      my.atac.sum$padj[i]  <- my.atac.peaks[my.min.p,]$padj
      
    } else {
      min.p.data <- my.atac.peaks[my.min.p,]
      my.max.stat  <- which.max(abs(min.p.data$stat))
      
      my.atac.sum$logFC[i] <- min.p.data[my.max.stat,]$log2FoldChange
      my.atac.sum$stat[i]  <- min.p.data[my.max.stat,]$stat
      my.atac.sum$padj[i]  <- min.p.data[my.max.stat,]$padj
    }
    
  }
  
  atac.geneList             <- my.atac.sum$stat
  names(atac.geneList)      <- my.atac.sum$GeneName
  atac.geneList             <- sort(atac.geneList , decreasing = TRUE)
  
  return(atac.geneList)
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

