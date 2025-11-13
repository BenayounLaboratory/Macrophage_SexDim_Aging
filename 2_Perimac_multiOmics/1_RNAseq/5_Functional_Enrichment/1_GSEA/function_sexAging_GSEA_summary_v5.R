options(stringsAsFactors = F)
library(ggplot2)
library(scales)
library(bitops)

# 
# my.m.gset    <- my.m.GO
# my.f.gset    <- my.f.GO
# my.gset.name <- "GO_ALL"
# my.thrs = 0.05

# my.m.gset    <- my.m.react
# my.f.gset    <- my.f.react
# my.gset.name <- "MSigDB_Reactome"
# my.thrs = 0.05

summarize_SexAge <- function(my.m.gset, my.f.gset, my.gset.name, my.thrs = 0.05){
  
  my.f.gset <- cbind(rownames(my.f.gset),my.f.gset)
  my.m.gset <- cbind(rownames(my.m.gset),my.m.gset)
  
  colnames(my.f.gset)[1] <- "Pathways"
  colnames(my.m.gset)[1] <- "Pathways"
  
  # merge to assess only pathways significant in both to identify divergent and convergent pathways
  my.merged.gset.both <- merge(my.f.gset, my.m.gset, by = "Pathways", suffixes = c(".F",".M"))
  
  # clean pathway names if needed (REACTOME)
  if (grepl("Reactome",my.gset.name, ignore.case = T)) {
    my.merged.gset.both$Pathways <- gsub("_"," ", my.merged.gset.both$Pathways)
    my.merged.gset.both$Pathways <- tolower(my.merged.gset.both$Pathways)
  }
  
  
  # filter for pathways significant with age in males OR females
  sig.filt <- bitOr(my.merged.gset.both$fdr.F < my.thrs, my.merged.gset.both$fdr.M < my.thrs) >0
  my.merged.gset.both <- my.merged.gset.both[sig.filt,]
  
  my.merged.gset.both$AGREEMENT <- ifelse(sign(my.merged.gset.both$nes.F) * sign(my.merged.gset.both$nes.M)>0,"CONVERGENT","DIVERGENT")
  
  pdf(paste0(Sys.Date(), "_GSEA_enrichment_", my.gset.name, "_GSEA_enrichment_PeriMac_Sex_Aging_FDR5_scatterplot.pdf"), width = 6, height = 6)
  plot(my.merged.gset.both$nes.F, 
       my.merged.gset.both$nes.M, 
       pch = 16,
       xlab = "GSEA NES with female aging",
       ylab = "GSEA NES with male aging",
       xlim = c(-4, 4),
       ylim = c(-4, 4),
       col = grDevices::adjustcolor( "black", alpha.f = 0.5),
       main = paste0(my.gset.name, " (FDR < ", 100* my.thrs, "% in either sex)") )
  abline(h = 0, col = "red", lty = "dashed")
  abline(v = 0, col = "red", lty = "dashed")
  dev.off()
  
  write.table(my.merged.gset.both, file = paste0(Sys.Date(), "_GSEA_enrichment_", my.gset.name, "_GSEA_enrichment_PeriMac_Sex_Aging_FDR", 100* my.thrs, "_Summary.txt"),
              sep = "\t", quote = F, row.names = F)
  
  
  #### Top divergent pathways
  path.disagree    <- my.merged.gset.both[my.merged.gset.both$AGREEMENT == "DIVERGENT"   ,]
  path.disagree$nes.F.abs <- abs(path.disagree$nes.F) # needed for sorting
  path.disagree$nes.M.abs <- abs(path.disagree$nes.M) # needed for sorting

  # lots of ties - need a sorting mechanism to report most interesting
  # sort on lowest p in one sex and highest in other, then on NES up or down
  top.div.m <- path.disagree[with(path.disagree, order(path.disagree$fdr.M, path.disagree$nes.M.abs, path.disagree$fdr.F,  decreasing = c(F,T,T))),] # inc pvalM, dec abs NES M, dec pval F
  top.div.f <- path.disagree[with(path.disagree, order(path.disagree$fdr.F, path.disagree$nes.F.abs, path.disagree$fdr.M,  decreasing = c(F,T,T))),] # inc pvalF, dec abs NES F, dec pval M
  
  
  # top 5 strongest enrichment/FDR up and down in F/M that disagree
  my.top.disagree  <- rbind(top.div.m[1:min(nrow(top.div.m),5), ] , # 5 or as many satisfying criteria
                            top.div.f[1:min(nrow(top.div.f),5), ] ) # 5 or as many satisfying criteria
  
  tmp.f <- my.top.disagree[,1:6]
  tmp.m <- my.top.disagree[,c(1,7:11)]
  colnames(tmp.f) <- c("Pathways","n","es","nes","pval.nes","fdr")
  colnames(tmp.m) <- c("Pathways","n","es","nes","pval.nes","fdr")
  
  # resort here so up and down are not mixed within each sex
  tmp.f <- tmp.f[with(tmp.f, order(tmp.f$fdr, tmp.f$nes, decreasing = T)),]
  tmp.m <- tmp.m[with(tmp.m, order(tmp.m$fdr, tmp.m$nes, decreasing = T)),]
  
  DIS.gsea.tab <- rbind(tmp.f, tmp.m)
  DIS.gsea.tab$Condition <- c(rep("F_aging",nrow(tmp.f)),rep("M_aging",nrow(tmp.m)))
  
  # create -log10 FDR for plotting
  DIS.gsea.tab$minlog10fdr  <- -log10(DIS.gsea.tab$fdr + 1e-30)
  
  # create and preserve wanted display order
  DIS.gsea.tab$Pathways <- factor(DIS.gsea.tab$Pathways, levels = unique(DIS.gsea.tab$Pathways))
  
  # Color scale
  my.max <- max(DIS.gsea.tab$nes)
  my.min <- min(DIS.gsea.tab$nes)
  
  my.values <- c(my.min,0.75*my.min,0.5*my.min,0.25*my.min,0,0.25*my.max,0.5*my.max,0.75*my.max,my.max)
  my.scaled <- rescale(my.values, to = c(0, 1))
  my.color.vector.age <- c("darkblue","dodgerblue4","dodgerblue3","dodgerblue1","white","lightcoral","brown1","firebrick2","firebrick4")
  
  my.plot <- ggplot(DIS.gsea.tab,aes(x=Condition,y=Pathways,colour=nes,size=minlog10fdr))+ theme_bw()+ geom_point(shape = 16)
  my.plot <- my.plot + ggtitle("PeriMac Aging") + labs(x = "-log10(pvalue)", y = "Sex DIVERGENT")
  my.plot <- my.plot + scale_colour_gradientn(colours = my.color.vector.age, na.value = "grey50", guide = "colourbar", values = my.scaled, limits = c(-3, 3))
  my.plot <- my.plot + scale_size_area(limits = c(0,30)) + scale_y_discrete(labels = wrap_format(40))
  my.plot
  
  pdf( paste(Sys.Date(),"GSEA_aging_BALLOON_plot",my.gset.name,"top_10_DIVERGENT_pathways_PeriMac_Sex_Aging.pdf", sep="_"), onefile=T, height = 4, width=5)
  print(my.plot)
  dev.off()
  
  
  #### Top convergent pathways
  path.agree      <- my.merged.gset.both[my.merged.gset.both$AGREEMENT == "CONVERGENT" ,]
  path.agree      <- path.agree[bitAnd(path.agree$fdr.F < my.thrs, path.agree$fdr.M < my.thrs) >0, ] # double significant ONLY
  
  # get average NES to select top
  path.agree$avNES <- apply(path.agree[,c("nes.F","nes.M")],1,mean)
  
  # top 5 up and down that agree; first NES, then F pval, then M pval
  top.conv.pos <- path.agree[with(path.agree, order(path.agree$avNES, path.agree$fdr.F, path.agree$fdr.M, decreasing = c(T,F,F))),]
  top.conv.neg <- path.agree[with(path.agree, order(path.agree$avNES, path.agree$fdr.F, path.agree$fdr.M, decreasing = c(F,F,F))),]
  
  my.top.agree  <- rbind(top.conv.pos[1:min(nrow(top.conv.pos),5),],  # largest value is top (positive)
                         top.conv.neg[1:min(nrow(top.conv.neg),5),])  # largest value is top (negative)
  
  tmp.f <- my.top.agree[,1:6]
  tmp.m <- my.top.agree[,c(1,7:11)]
  colnames(tmp.f) <- c("Pathways","n","es","nes","pval.nes","fdr")
  colnames(tmp.m) <- c("Pathways","n","es","nes","pval.nes","fdr")
  
  AGR.gsea.tab <- rbind(tmp.f, tmp.m)
  AGR.gsea.tab$Condition <- c(rep("F_aging",nrow(tmp.f)),rep("M_aging",nrow(tmp.m)))
  
  # create -log10 FDR for plotting
  AGR.gsea.tab$minlog10fdr  <- -log10(AGR.gsea.tab$fdr + 1e-30)
  
  # resort here so data looks cleanr
  AGR.gsea.tab <- AGR.gsea.tab[with(AGR.gsea.tab, order(AGR.gsea.tab$nes, AGR.gsea.tab$fdr, decreasing = c(F,T))),]

  # create and preserve wanted display order
  AGR.gsea.tab$Pathways <- factor(AGR.gsea.tab$Pathways, levels = unique(AGR.gsea.tab$Pathways))
  
  # Color scale
  my.max <- max(AGR.gsea.tab$nes)
  my.min <- min(AGR.gsea.tab$nes)
  
  my.values <- c(my.min,0.75*my.min,0.5*my.min,0.25*my.min,0,0.25*my.max,0.5*my.max,0.75*my.max,my.max)
  my.scaled <- rescale(my.values, to = c(0, 1))
  my.color.vector.age <- c("darkblue","dodgerblue4","dodgerblue3","dodgerblue1","white","lightcoral","brown1","firebrick2","firebrick4")
  
  my.plot <- ggplot(AGR.gsea.tab,aes(x=Condition,y=Pathways,colour=nes,size=minlog10fdr))+ theme_bw()+ geom_point(shape = 16)
  my.plot <- my.plot + ggtitle("PeriMac Aging (FDR <5%)") + labs(x = "-log10(pvalue)", y = "Sex CONVERGENT")
  my.plot <- my.plot + scale_colour_gradientn(colours = my.color.vector.age, na.value = "grey50", guide = "colourbar", values = my.scaled, limits = c(-3, 3))
  my.plot <- my.plot + scale_size_area(limits = c(0,30)) + scale_y_discrete(labels = wrap_format(40))
  my.plot
  
  pdf( paste(Sys.Date(),"GSEA_aging_BALLOON_plot",my.gset.name,"top_10_CONVERGENT_pathways_PeriMac_Sex_Aging.pdf", sep="_"), onefile=T, height = 4, width=5)
  print(my.plot)
  dev.off()
  
}

