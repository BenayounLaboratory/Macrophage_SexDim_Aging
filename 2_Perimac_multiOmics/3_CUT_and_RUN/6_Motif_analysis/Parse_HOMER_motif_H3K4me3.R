setwd('/Volumes/BB_HQ_1/Immune_sex_dimorphism_Aging/Macrophages/Clean_Macrophage_Data/CUT_and_RUN/H3K4me3_Aging/H3K4me3_height/HOMER_motif')
options(stringsAsFactors = F)

# 2022-12-20
# parse/plot H3K4me3 motif enrichment data

############################################
# 1. Read in HOMER motif enrichment in sex-specifc changes with age H3K4 peaks
H3K4.motif.M.dwn  <- read.csv("2022-04-12_HOMER_MOTIF_H3K4me3_M_DWN_ONLY_FDR5/knownResults.txt" , sep = "\t", header = T)
H3K4.motif.M.up   <- read.csv("2022-04-12_HOMER_MOTIF_H3K4me3_M_UP_ONLY_FDR5/knownResults.txt"  , sep = "\t", header = T)
H3K4.motif.F.dwn  <- read.csv("2022-04-12_HOMER_MOTIF_H3K4me3_F_DWN_ONLY_FDR5/knownResults.txt" , sep = "\t", header = T)
H3K4.motif.F.up   <- read.csv("2022-04-12_HOMER_MOTIF_H3K4me3_F_UP_ONLY_FDR5/knownResults.txt"  , sep = "\t", header = T)

# get clean colnames
colnames(H3K4.motif.M.up ) <- c("Motif.Name"            ,
                                "Consensus"             ,
                                "pvalue"                ,
                                "Log_pvalue"            ,
                                "q.value"               ,
                                "Num_FG_with_motif"     ,
                                "Percent_FG_with_motif" ,
                                "Num_BG_with_motif"     ,
                                "Percent_BG_with_motif" )
colnames(H3K4.motif.M.dwn) <- colnames(H3K4.motif.M.up)
colnames(H3K4.motif.F.up ) <- colnames(H3K4.motif.M.up)
colnames(H3K4.motif.F.dwn) <- colnames(H3K4.motif.M.up)

# extract percentages as numbers
H3K4.motif.M.up$Percent_FG_with_motif  <- as.numeric(sub("%", "", H3K4.motif.M.up$Percent_FG_with_motif )) # get percent number
H3K4.motif.M.up$Percent_BG_with_motif  <- as.numeric(sub("%", "", H3K4.motif.M.up$Percent_BG_with_motif )) # get percent number
H3K4.motif.M.dwn$Percent_FG_with_motif <- as.numeric(sub("%", "", H3K4.motif.M.dwn$Percent_FG_with_motif)) # get percent number
H3K4.motif.M.dwn$Percent_BG_with_motif <- as.numeric(sub("%", "", H3K4.motif.M.dwn$Percent_BG_with_motif)) # get percent number
H3K4.motif.F.up$Percent_FG_with_motif  <- as.numeric(sub("%", "", H3K4.motif.F.up$Percent_FG_with_motif )) # get percent number
H3K4.motif.F.up$Percent_BG_with_motif  <- as.numeric(sub("%", "", H3K4.motif.F.up$Percent_BG_with_motif )) # get percent number
H3K4.motif.F.dwn$Percent_FG_with_motif <- as.numeric(sub("%", "", H3K4.motif.F.dwn$Percent_FG_with_motif)) # get percent number
H3K4.motif.F.dwn$Percent_BG_with_motif <- as.numeric(sub("%", "", H3K4.motif.F.dwn$Percent_BG_with_motif)) # get percent number

H3K4.motif.M.up$EnrichFold  <- H3K4.motif.M.up$Percent_FG_with_motif/H3K4.motif.M.up$Percent_BG_with_motif
H3K4.motif.M.dwn$EnrichFold <- H3K4.motif.M.dwn$Percent_FG_with_motif/H3K4.motif.M.dwn$Percent_BG_with_motif
H3K4.motif.F.up$EnrichFold  <- H3K4.motif.F.up$Percent_FG_with_motif/H3K4.motif.F.up$Percent_BG_with_motif
H3K4.motif.F.dwn$EnrichFold <- H3K4.motif.F.dwn$Percent_FG_with_motif/H3K4.motif.F.dwn$Percent_BG_with_motif

# Filter to retain only significantly enriched
sig.H3K4.motif.M.up  <- H3K4.motif.M.up  [H3K4.motif.M.up  $q.value < 0.05,]  # 0 --- no male motifs
sig.H3K4.motif.M.dwn <- H3K4.motif.M.dwn [H3K4.motif.M.dwn $q.value < 0.05,]  # 0 --- no male motifs
sig.H3K4.motif.F.up  <- H3K4.motif.F.up  [H3K4.motif.F.up  $q.value < 0.05,]  # 115
sig.H3K4.motif.F.dwn <- H3K4.motif.F.dwn [H3K4.motif.F.dwn $q.value < 0.05,]  # 64


# Put up and down together
sig.H3K4.motif.F                <- rbind(sig.H3K4.motif.F.up,sig.H3K4.motif.F.dwn)
sig.H3K4.motif.F$Direction      <- c(rep("UP",nrow(sig.H3K4.motif.F.up)),rep("DOWN",nrow(sig.H3K4.motif.F.dwn)))
sig.H3K4.motif.F$minusLog10pval <- -sig.H3K4.motif.F$Log_pvalue
sig.H3K4.motif.F$EnrichFold[-c(1:nrow(sig.H3K4.motif.F.up))] <-  -sig.H3K4.motif.F$EnrichFold[-c(1:nrow(sig.H3K4.motif.F.up))]   # make negative values for enrichment in regions of decreased accessibility

write.table(sig.H3K4.motif.F, file = paste0(Sys.Date(),"_HOMER_Motifs_Enriched_in_Female_ONLY_H3K4me3_DA_peaks_FDR5.txt"), quote = F, sep = "\t", row.names = F)

# grab motifs related to transcription factor with RNAseq expression changes
F.age.TFs <- c("Meis1"   ,
               "HPC7-Scl",   # Tal1
               "IRF2"    ,
               "TEAD1"   )


############################################################################################
# 2. Make bubble chart summary
library(ggplot2) 
library(scales) 
theme_set(theme_bw())

# get filtered/merged datafame for ggplot

sig.H3K4.motif.F.filt <- sig.H3K4.motif.F[grep(paste(F.age.TFs,collapse="|"), sig.H3K4.motif.F$Motif.Name),]

F.srt <- sort(sig.H3K4.motif.F.filt$EnrichFold,decreasing = T, index.return = T)

filt.F.srt <- sig.H3K4.motif.F.filt[F.srt$ix,]

filt.F.srt$Condition <- "F_ONLY_H3K4_DA"

#### Females
my.max <- max(filt.F.srt$EnrichFold)
my.min <- min(filt.F.srt$EnrichFold)
my.values <- c(my.min,0.75*my.min,0.5*my.min,0.25*my.min,0,0.25*my.max,0.5*my.max,0.75*my.max,my.max)
my.scaled <- rescale(my.values, to = c(0, 1))
my.color.vector <- c("darkblue","dodgerblue4","dodgerblue3","dodgerblue1","white","lightcoral","brown1","firebrick2","firebrick4")

# to preserve the wanted order
filt.F.srt$Motif.Name  <- factor(filt.F.srt$Motif.Name, levels = rev(unique(filt.F.srt$Motif.Name)))

pdf(paste0(Sys.Date(),"HOMER_Motifs_Enriched_H3K4me3_DA_FEMALE_ONLY_SelectedDETF_Aging_Perimac_FDR5.pdf"), height = 3, width=9)
my.plot <- ggplot(filt.F.srt,aes(x=Condition,y=Motif.Name,colour=EnrichFold,size=minusLog10pval))+ theme_bw()+ geom_point(shape = 16)
my.plot <- ggplot(filt.F.srt,aes(x=Condition,y=Motif.Name,colour=EnrichFold,size=minusLog10pval))+ theme(text = element_text(size=16))+ geom_point(shape = 16)
my.plot <- my.plot + ggtitle("HOMER") + labs(x = "PeriMac Aging", y = "Sex-Specific DA H3K4")
my.plot <- my.plot + scale_colour_gradientn(colours = my.color.vector,space = "Lab", na.value = "grey50", guide = "colourbar", values = my.scaled)
print(my.plot)
dev.off()  



