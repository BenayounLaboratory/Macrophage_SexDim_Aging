setwd('/Users/berenice/Desktop/Irf2_CUTRUN/GREAT/ChIP_CUTRUN_summary')
options(stringsAsFactors = FALSE)

# load libraries for plotting
library(ggplot2)
library(scales) 
library(bitops) 

options(java.parameters = "-Xmx16g" )
require(openxlsx)

# 2025-12-02
# Analyze/Plot Irf2 ChIP and CUT&RUN functional enrichment by GREAT

# 2025-12-07
# Use GREAT results on HOMER peaks to match CUTadnRUN


######################################################################
# 1. read and parse GREAT output
PeriMac.great.res <- read.csv('2025-12-01_Irf2_CUTRUN_greatExportAll.tsv', header = T, sep = "\t", skip = 3, nrows = 3500) # nrow -> to not read comment lines at bottom
BMDM.great.res    <- read.csv('2025-12-07_Irf2_ChIPseq_greatExportAll.tsv', header = T, sep = "\t", skip = 3, nrows = 3500) # nrow -> to not read comment lines at bottom

# Ontologies
unique(PeriMac.great.res$X..Ontology)
# [1] Ensembl Genes             GO Biological Process     GO Cellular Component     GO Molecular Function     Human Phenotype           Mouse Phenotype Single KO Mouse Phenotype  

PeriMac.great.sig <- PeriMac.great.res[PeriMac.great.res$HyperFdrQ < 0.05, c("X..Ontology","ID","Desc", "HyperRank", "HyperP", "HyperFdrQ", "GeneFoldEnrich" )] # 721 terms are significant
BMDM.great.sig    <- BMDM.great.res   [BMDM.great.res   $HyperFdrQ < 0.05, c("X..Ontology","ID","Desc", "HyperRank", "HyperP", "HyperFdrQ", "GeneFoldEnrich" )] # 486 terms are significant

my.GO.cat <- c("GO Biological Process",
               "GO Cellular Component",
               "GO Molecular Function")

# grab GO, all categories
PeriMac.goall <- PeriMac.great.sig[PeriMac.great.sig$X..Ontology %in% my.GO.cat, ] # 265
BMDM.goall    <- BMDM.great.sig[BMDM.great.sig$X..Ontology %in% my.GO.cat, ]       # 224
######################################################################

######################################################################
# 2. get enrichment plots

# get conserved categories
merged.goall  <- merge(PeriMac.goall, BMDM.goall[,-c(1,3)], by = "ID", suffixes = c(".p",".b")) # 127

merged.goall$av_enrich <- apply(merged.goall[,c("GeneFoldEnrich.p","GeneFoldEnrich.b")],1,mean)
merged.goall$av_fdr    <- apply(merged.goall[,c("HyperFdrQ.p"     ,"HyperFdrQ.b")],1,mean)

# sort on lowest p, then on enrichment
top.cat <- merged.goall[with(merged.goall, order(merged.goall$av_enrich,  merged.goall$av_fdr, decreasing = c(T,F))),] # inc fdr, dec enrichment

# get top 10 catagories
top.cat.10 <- top.cat[1:10,]

# generate ggplot output format
top.bmdm <- top.cat.10[,c(1,3,8:11)]
top.peri <- top.cat.10[,c(1,3:7)]

colnames(top.bmdm) <- c("ID", "Desc", "HyperRank", "HyperP", "HyperFdrQ","GeneFoldEnrich")
colnames(top.peri) <- c("ID", "Desc", "HyperRank", "HyperP", "HyperFdrQ","GeneFoldEnrich")

top.bmdm$Condition <- "BMDM"
top.peri$Condition <- "PeriMac"

top.plot <- rbind(top.bmdm,top.peri)

# create -log10 FDR for plotting
top.plot$minlog10fdr  <- -log10(top.plot$HyperFdrQ + 1e-30)

# create and preserve wanted display order
top.plot$Pathways <- factor(paste0(top.plot$ID,"_",top.plot$Desc), levels = unique(paste0(top.plot$ID,"_",top.plot$Desc)))

# Color scale
my.max <- 25

my.values <- c(0,0.25*my.max,0.5*my.max,0.75*my.max,my.max)
my.scaled <- rescale(my.values, to = c(0, 1))
my.color.vector <- c("mistyrose","lightcoral","brown1","brown2","firebrick2","firebrick3","firebrick4")

my.plot <- ggplot(top.plot,aes(x=Condition,y=Pathways,colour=GeneFoldEnrich,size=minlog10fdr))+ theme_bw()+ geom_point(shape = 16)
my.plot <- my.plot + ggtitle("Irf2 GREAT enrichments") + labs(x = "", y = "")
my.plot <- my.plot + scale_colour_gradientn(colours = my.color.vector, na.value = "grey50", guide = "colourbar", values = my.scaled, limits = c(0,my.max ))
my.plot <- my.plot + scale_size_area(limits = c(0,8)) + scale_y_discrete(labels = wrap_format(40))
my.plot

pdf( paste(Sys.Date(),"GREAT_BALLOON_plot_GOALL_top_10_pathways_Irf2_BMDM_PeriMac.pdf", sep="_"), onefile=T, height = 4, width = 5)
print(my.plot)
dev.off()
################################################################################################


######################################################################
# 3. export enrichment results

# export results
great.res <- list("BMDM_GO"       = BMDM.goall     ,
                  "PeriMac_GO"    = PeriMac.goall  ,
                  "Conserved_GO"  = top.cat        )

write.xlsx(great.res, rowNames = F, file = paste0(Sys.Date(),"_Irf2_ChIP_CUTRUN_GREAT_GOALL_FDR5_enrichment_Results.xlsx"))


################################################################################################

#######################
sink(file = paste(Sys.Date(),"_session_Info.txt", sep =""))
sessionInfo()
sink()

