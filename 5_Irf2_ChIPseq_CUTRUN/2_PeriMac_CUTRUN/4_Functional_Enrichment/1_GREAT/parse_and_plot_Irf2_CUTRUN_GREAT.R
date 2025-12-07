setwd('/Users/berenice/Desktop/Irf2_CUTRUN/GREAT/All_peaks')
options(stringsAsFactors = FALSE)

# load libraries for plotting
library(ggplot2)
library(scales) 
library(bitops) 

source('get_bubbles_GREAT.R')

options(java.parameters = "-Xmx16g" )
require(openxlsx)


# 2025-12-02
# run Irf2 CUT and RUN functional enrichment

######################################################################
# 1. read and parse GREAT output
my.great.results <- read.csv('2025-12-01_Irf2_CUTRUN_greatExportAll.tsv', header = T, sep = "\t", skip = 3, nrows = 3500) # nrow -> to not read comment lines at bottom
# 3500 terms

# Ontologies
unique(my.great.results$X..Ontology)
# [1] Ensembl Genes             GO Biological Process     GO Cellular Component     GO Molecular Function     Human Phenotype           Mouse Phenotype Single KO
# [7] Mouse Phenotype  

my.great.results.sig <- my.great.results[my.great.results$HyperFdrQ < 0.05, c("X..Ontology","ID","Desc", "HyperRank", "HyperP", "HyperFdrQ", "GeneFoldEnrich" )] # 721 terms are significant

great.gobp <- my.great.results.sig[my.great.results.sig$X..Ontology == "GO Biological Process"    ,]  
great.gocc <- my.great.results.sig[my.great.results.sig$X..Ontology == "GO Cellular Component"    ,]  
great.gomf <- my.great.results.sig[my.great.results.sig$X..Ontology == "GO Molecular Function"    ,]  
great.mpo  <- my.great.results.sig[my.great.results.sig$X..Ontology == "Mouse Phenotype Single KO",]  

great.goall <- rbind(great.gobp,
                     great.gocc,
                     great.gomf)



######################################################################
# 2. get enrichment plots
plot_enrich_ChIP(great.gobp , "GO_BP", my.thrs = 0.05)
plot_enrich_ChIP(great.gocc , "GO_CC", my.thrs = 0.05)
plot_enrich_ChIP(great.gomf , "GO_MF", my.thrs = 0.05)
plot_enrich_ChIP(great.mpo  , "Mouse_Phenotypes", my.thrs = 0.05)

plot_enrich_ChIP(great.goall , "GO_ALL", my.thrs = 0.05)

################################################################################################


######################################################################
# 3. export enrichment results

# export results
great.res <- list("GO_BP"         = great.gobp ,
                  "GO_CC"         = great.gocc ,
                  "GO_MF"         = great.gomf ,
                  "Mouse_pheno"   = great.mpo  )
write.xlsx(great.res, rowNames = F, file = paste0(Sys.Date(),"_Irf2_CUTRUN_GREAT_FDR5_enrichment_Results.xlsx"))


################################################################################################

#######################
sink(file = paste(Sys.Date(),"_session_Info.txt", sep =""))
sessionInfo()
sink()

