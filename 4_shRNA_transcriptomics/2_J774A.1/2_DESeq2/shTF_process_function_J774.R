
# ### DEBUG
# my.raw.data <- data.irf2
# my.TF       <- "Irf2"
# shNon = 1:3
# shLuc = 4:6
# shTF1 = 7:9
# shTF2 = 10:12


process_shTF <- function (my.raw.data, my.TF, shNon = 1:3, shLuc = 4:6, shTF1 = 7:9, shTF2 = 10:12) {
  
  # remove non essential columns
  my.raw.data <- my.raw.data[,-c(2:6)]
  
  # Clean up names
  rownames(my.raw.data) <- my.raw.data$Geneid
  
  ## address all possible variants in bam file names
  colnames(my.raw.data) <- gsub("_RNAseq_STAR_Aligned.sortedByCoord.out.bam","",   colnames(my.raw.data))
  
  colnames(my.raw.data)[1] <- c("Gene_Symbol")
  
  # only keep genes with at least 1 read in at least half of the samples
  my.good.data <- which(apply(my.raw.data[,-1]>0, 1, sum) >= ncol(my.raw.data[,-1])/2 ) 
  my.matrix3   <- my.raw.data[my.good.data,-1] # 14858 genes for Irf2 example
  
  # get output file prefixes
  my.outprefix <- paste(Sys.Date(),"J774_sh", my.TF ,"data_DESeq2",sep="_")
  
  # covariate building
  my.treat        <- rep(NA,ncol(my.matrix3))
  my.treat[shNon] <- "CTL"
  my.treat[shLuc] <- "CTL"
  my.treat[shTF1] <- "sh1"
  my.treat[shTF2] <- "sh2"
  
  # design matrix
  dataDesign = data.frame( row.names  = colnames(my.matrix3), 
                           treat      = my.treat)
  
  # get matrix using age and genotype as a modeling covariate
  dds <- DESeqDataSetFromMatrix(countData = my.matrix3,
                                colData   = dataDesign,
                                design    = ~ treat)
  
  # run DESeq normalizations and export results
  dds.deseq <- DESeq(dds)
  
  # plot dispersion
  my.disp.out <- paste(my.outprefix,"_dispersion_plot.pdf")
  
  pdf(my.disp.out)
  plotDispEsts(dds.deseq)
  dev.off()
  
  # normalized expression value
  norm.cts <- getVarianceStabilizedData(dds.deseq)
  
  # output result tables of combined analysis to text files
  my.out.ct.mat <- paste0(my.outprefix,"_log2_counts_matrix.txt")
  write.table(norm.cts, file = my.out.ct.mat , sep = "\t" , row.names = T, quote=F)
  
  # color-code 
  my.colors       <-  c(rep("orange",length(shNon)), rep("peru",length(shLuc)), rep("seagreen",length(shTF1)), rep("seagreen3",length(shTF2)))
  
  # expression range
  pdf(paste(my.outprefix,"_Normalized_counts_boxplot.pdf"))
  boxplot(norm.cts,col=my.colors,cex=0.5,ylab="Log2 DESeq2 Normalized counts", las = 2)  
  dev.off()
  
  # do MDS analysis
  mds.result <- cmdscale(1-cor(norm.cts,method="spearman"), k = 2, eig = FALSE, add = FALSE, x.ret = FALSE)
  x <- mds.result[, 1]
  y <- mds.result[, 2]
  
  my.mds.out <- paste0(my.outprefix,"MDS_plot.pdf")
  pdf(my.mds.out, width = 5, height = 5)
  plot(x, y,
       xlab = "MDS dimension 1", ylab = "MDS dimension 2",
       main=paste0("MDS sh",my.TF),
       cex=2, pch = 15, col = my.colors,
       las =1)
  points(x, y, cex = 2, pch = 0)
  legend("topright",c("shNT","shLuc",paste0("sh",my.TF,"_1"),paste0("sh",my.TF,"_2")), col = c("orange","peru","seagreen","seagreen3"), pch = 15, cex = 0.75, bty = 'n')
  dev.off()
  
  # plot KD gene expression
  TF.exp <- list("shControl"   = norm.cts[my.TF, c(shNon, shLuc) ],
                 "sh_1"        = norm.cts[my.TF, shTF1 ]          ,
                 "sh_2"        = norm.cts[my.TF, shTF2 ]          ) 
  
  my.max.exp <- max (norm.cts[my.TF,])
  my.min.exp <- min (norm.cts[my.TF,])
  
  ymax <- ceiling(my.max.exp * 1.1)
  ymin <- floor(my.min.exp * 0.9)
  
  pdf(paste0(my.outprefix,"_Normalized_counts_boxplot_", my.TF,".pdf"), height = 5, width = 3)
  boxplot(TF.exp, outline = F, las = 1, ylim = c(ymin,ymax) , main = my.TF, ylab = "DEseq2 VST normalized log2 counts", cex.axis = 0.75)
  beeswarm(TF.exp, add = T, pch = 15, pwcol = my.colors)
  dev.off()
  
  
  ################################################################
  # contrast KD genotypes
  
  res.shTF1 <- results(dds.deseq, contrast = c("treat", "sh1", "CTL"))
  res.shTF2 <- results(dds.deseq, contrast = c("treat", "sh2", "CTL"))
  
  # exclude any NA values
  res.shTF1   <- res.shTF1[!is.na(res.shTF1$padj),]
  res.shTF2   <- res.shTF2[!is.na(res.shTF2$padj),]
  
  genes.shTF1 <- rownames(res.shTF1)[res.shTF1$padj < 0.05] # 3013 for Irf2
  genes.shTF2 <- rownames(res.shTF2)[res.shTF2$padj < 0.05] # 2024 for Irf2
  
  # get merge from 2 shRNA, combined fisher p-values
  my.merge.TF <- merge(data.frame(res.shTF1), data.frame(res.shTF2), by = "row.names", suffixes = c(".sh1",".sh2"))
  my.merge.TF$Sign <- sign(my.merge.TF$log2FoldChange.sh1) * sign(my.merge.TF$log2FoldChange.sh2)
  my.merge.TF$p_combined <- my.merge.TF$padj.sh1 * my.merge.TF$padj.sh2
  
  # Significant in both with the same direction
  # 0.01*0.01
  # use 1e-04 as threshold
  genes.TF  <- my.merge.TF$Row.names[bitAnd(my.merge.TF$p_combined < 1e-04, my.merge.TF$Sign==1) >0]
  my.num.TF <- length(genes.TF) # 1263 genes for Irf2
  
  # heatmap drawing - only if there is at least 2 gene
  my.heatmap.out <- paste0(my.outprefix,"sh", my.TF, "_CONSISTENT_regulated_Heatmap_FDR1e-4.pdf")
  
  pdf(my.heatmap.out, onefile = F)
  my.heatmap.title <- paste(my.TF, " significant (FDR<1e-04), ", my.num.TF, " genes",sep="")
  pheatmap(norm.cts[genes.TF,],
           cluster_cols = F,
           cluster_rows = T,
           colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
           show_rownames = F, scale="row",
           main = my.heatmap.title, 
           cellwidth = 15,
           border = NA)
  dev.off()
  
  write.table(my.merge.TF, file = paste0(my.outprefix,"_sh", my.TF, "_DEseq2_results_COMBINED_all_genes_statistics.txt"), col.names = T, row.names = F, quote = F, sep = "\t")
  write.table(my.merge.TF[my.merge.TF$Row.names %in% genes.TF,], file = paste0(my.outprefix,"_sh", my.TF, "_DEseq2_results_COMBINED_FDR1e-4_genes_statistics.txt"), col.names = T, row.names = F, quote = F, sep = "\t")
  
  
}
