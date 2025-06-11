options(stringsAsFactors = F)

# debug
# great.res    <- great.gobp
# my.gset.name <- "GOBP"
# my.thrs      <- 0.05


# great.res    <- great.gobp
# my.gset.name <- "GOBP"
# my.thrs      <- 0.05

# fix ordering

plot_enrich_ChIP <- function(great.res, my.gset.name, my.thrs = 0.05) {
  
  # top 10 enriched pathways; first FDR, then enrichment
  top.gsets <- great.res[with(great.res, order(great.res$HyperFdrQ, great.res$GeneFoldEnrich, decreasing = c(F,T))),]
  
  my.top.agree  <- top.gsets[1:min(nrow(top.gsets),10),]  # largest value is top (positive)
  my.top.agree$Condition <- rep("Irf2_Bound",nrow(my.top.agree))
  my.top.agree$Pathways  <- paste(my.top.agree$ID,my.top.agree$Desc, sep = "_")
    
  # create -log10 FDR for plotting
  my.top.agree$minlog10fdr  <- -log10(my.top.agree$HyperFdrQ + 1e-30)
  
  # create and preserve wanted display order
  my.max.char <- max(nchar(my.top.agree$Pathways))
  my.top.agree$Pathways <- factor(my.top.agree$Pathways, levels = unique(my.top.agree$Pathways))
  
  # Color scale
  my.max <- 40

  my.values <- c(0,0.25*my.max,0.5*my.max,0.75*my.max,my.max)
  my.scaled <- rescale(my.values, to = c(0, 1))
  my.color.vector.age <- c("mistyrose","lightcoral","brown1","firebrick2","firebrick4")
  
  my.plot <- ggplot(my.top.agree,aes(x=Condition,y=Pathways,colour=GeneFoldEnrich,size=minlog10fdr))+ theme_bw()+ geom_point(shape = 16)
  my.plot <- my.plot + ggtitle("Irf2 ChIPseq GREAT enrichments") + labs(x = "", y = "")
  my.plot <- my.plot + scale_colour_gradientn(colours = my.color.vector.age, na.value = "grey50", guide = "colourbar", values = my.scaled, limits = c(0,my.max ))
  my.plot <- my.plot + scale_size_area(limits = c(0,10)) + scale_y_discrete(labels = wrap_format(40))
  my.plot
  
  pdf( paste(Sys.Date(),"GREAT_BALLOON_plot",my.gset.name,"top_10_pathways_Irf2_BMDM_ChIP.pdf", sep="_"), onefile=T, height = 4, width = 5)
  print(my.plot)
  dev.off()
  
}