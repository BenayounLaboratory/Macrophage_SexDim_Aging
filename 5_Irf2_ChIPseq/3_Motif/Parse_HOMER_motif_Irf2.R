setwd('/Volumes/BB_HQ_1/Immune_sex_dimorphism_Aging/Macrophages/Clean_Macrophage_Data/ChIP-seq/HOMER_Motif')
options(stringsAsFactors = F)

# 2025-01-16
# parse/plot Irf2 ChIPseq motif enrichment data

############################################
# 1. Read in HOMER motif enrichment in sex-specifc changes with age ATAC peaks
Irf2.chip.motifs  <- read.csv("2025-01-16_HOMER_MOTIF_Irf2_BMDM_ConsensusPeaks/knownResults.txt" , sep = "\t", header = T)

# get clean colnames
colnames(Irf2.chip.motifs ) <- c("Motif.Name"            ,
                                "Consensus"             ,
                                "pvalue"                ,
                                "Log_pvalue"            ,
                                "q.value"               ,
                                "Num_FG_with_motif"     ,
                                "Percent_FG_with_motif" ,
                                "Num_BG_with_motif"     ,
                                "Percent_BG_with_motif" )


# extract percentages as numbers
Irf2.chip.motifs$Percent_FG_with_motif  <- as.numeric(sub("%", "", Irf2.chip.motifs$Percent_FG_with_motif )) # get percent number
Irf2.chip.motifs$Percent_BG_with_motif  <- as.numeric(sub("%", "", Irf2.chip.motifs$Percent_BG_with_motif )) # get percent number

Irf2.chip.motifs$EnrichFold             <- Irf2.chip.motifs$Percent_FG_with_motif/Irf2.chip.motifs$Percent_BG_with_motif
Irf2.chip.motifs$minusLog10pval         <- -Irf2.chip.motifs$Log_pvalue

# Filter to retain only significantly enriched
sig.Irf2.chip.motifs  <- Irf2.chip.motifs  [Irf2.chip.motifs  $q.value < 0.05,]  # 46

write.table(sig.Irf2.chip.motifs, file = paste0(Sys.Date(),"_HOMER_Motifs_Enriched_in_Irf2_chIPseq_FDR5.txt"), quote = F, sep = "\t", row.names = F)


############################################################################################
# 2. Make bubble chart summary
library(ggplot2) 
library(scales) 
theme_set(theme_bw())

# get filtered/merged datafame for ggplot
sort.order    <- sort(sig.Irf2.chip.motifs$minusLog10pval, decreasing = T, index.return = T)
filt.Irf2.srt <- sig.Irf2.chip.motifs[sort.order$ix[1:5],]

filt.Irf2.srt$Condition <- "Irf2_BMDM_ChIP"

#### plotting
my.max <- max(filt.Irf2.srt$EnrichFold)
my.values <- c(0,0.25*my.max,0.5*my.max,0.75*my.max,my.max)
my.scaled <- rescale(my.values, to = c(0, 1))
my.color.vector <- c("mistyrose","lightcoral","brown1","firebrick2","firebrick4")

# to preserve the wanted order
filt.Irf2.srt$Motif.Name  <- factor(filt.Irf2.srt$Motif.Name, levels = rev(unique(filt.Irf2.srt$Motif.Name)))

pdf(paste0(Sys.Date(),"HOMER_Motifs_Enriched_Top5_Irf2_BMDM_ChIPseq.pdf"), height = 5, width=8)
my.plot <- ggplot(filt.Irf2.srt,aes(x=Condition,y=Motif.Name,colour=EnrichFold,size=minusLog10pval))+ theme_bw()+ geom_point(shape = 16)
my.plot <- ggplot(filt.Irf2.srt,aes(x=Condition,y=Motif.Name,colour=EnrichFold,size=minusLog10pval))+ theme(text = element_text(size=16))+ geom_point(shape = 16)
my.plot <- my.plot + ggtitle("HOMER") + labs(x = "", y = "Enriched Motifs")
my.plot <- my.plot + scale_colour_gradientn(colours = my.color.vector,space = "Lab", na.value = "grey50", guide = "colourbar", values = my.scaled, limits = c(0,my.max))
my.plot <- my.plot + scale_y_discrete(labels = wrap_format(30)) + scale_size_area(limits = c(0,1000))
print(my.plot)
dev.off()  

##############################################################################################

#######################
sink(file = paste(Sys.Date(),"_session_Info.txt", sep =""))
sessionInfo()
sink()


