setwd('/Volumes/BB_Home_HQ/Immune_sex_dimorphism_Aging/Macrophages/Clean_Macrophage_Data/ATACseq/HOMER_Motif/')
options(stringsAsFactors = F)

# 2022-04-12
# parse/plot ATAC motif enrichment data


############################################
# 1. Read in HOMER motif enrichment in sex-specifc changes with age ATAC peaks
ATAC.motif.M.dwn  <- read.csv("2021-11-09_HOMER_MOTIF_ATAC_M_DWN_ONLY_FDR5/knownResults.txt" , sep = "\t", header = T)
ATAC.motif.M.up   <- read.csv("2021-11-09_HOMER_MOTIF_ATAC_M_UP_ONLY_FDR5/knownResults.txt"  , sep = "\t", header = T)
ATAC.motif.F.dwn  <- read.csv("2021-11-09_HOMER_MOTIF_ATAC_F_DWN_ONLY_FDR5/knownResults.txt" , sep = "\t", header = T)
ATAC.motif.F.up   <- read.csv("2021-11-09_HOMER_MOTIF_ATAC_F_UP_ONLY_FDR5/knownResults.txt"  , sep = "\t", header = T)

# get clean colnames
colnames(ATAC.motif.M.up ) <- c("Motif.Name"            ,
                                "Consensus"             ,
                                "pvalue"                ,
                                "Log_pvalue"            ,
                                "q.value"               ,
                                "Num_FG_with_motif"     ,
                                "Percent_FG_with_motif" ,
                                "Num_BG_with_motif"     ,
                                "Percent_BG_with_motif" )
colnames(ATAC.motif.M.dwn) <- colnames(ATAC.motif.M.up)
colnames(ATAC.motif.F.up ) <- colnames(ATAC.motif.M.up)
colnames(ATAC.motif.F.dwn) <- colnames(ATAC.motif.M.up)

# extract percentages as numbers
ATAC.motif.M.up$Percent_FG_with_motif  <- as.numeric(sub("%", "", ATAC.motif.M.up$Percent_FG_with_motif )) # get percent number
ATAC.motif.M.up$Percent_BG_with_motif  <- as.numeric(sub("%", "", ATAC.motif.M.up$Percent_BG_with_motif )) # get percent number
ATAC.motif.M.dwn$Percent_FG_with_motif <- as.numeric(sub("%", "", ATAC.motif.M.dwn$Percent_FG_with_motif)) # get percent number
ATAC.motif.M.dwn$Percent_BG_with_motif <- as.numeric(sub("%", "", ATAC.motif.M.dwn$Percent_BG_with_motif)) # get percent number
ATAC.motif.F.up$Percent_FG_with_motif  <- as.numeric(sub("%", "", ATAC.motif.F.up$Percent_FG_with_motif )) # get percent number
ATAC.motif.F.up$Percent_BG_with_motif  <- as.numeric(sub("%", "", ATAC.motif.F.up$Percent_BG_with_motif )) # get percent number
ATAC.motif.F.dwn$Percent_FG_with_motif <- as.numeric(sub("%", "", ATAC.motif.F.dwn$Percent_FG_with_motif)) # get percent number
ATAC.motif.F.dwn$Percent_BG_with_motif <- as.numeric(sub("%", "", ATAC.motif.F.dwn$Percent_BG_with_motif)) # get percent number

ATAC.motif.M.up$EnrichFold  <- ATAC.motif.M.up$Percent_FG_with_motif/ATAC.motif.M.up$Percent_BG_with_motif
ATAC.motif.M.dwn$EnrichFold <- ATAC.motif.M.dwn$Percent_FG_with_motif/ATAC.motif.M.dwn$Percent_BG_with_motif
ATAC.motif.F.up$EnrichFold  <- ATAC.motif.F.up$Percent_FG_with_motif/ATAC.motif.F.up$Percent_BG_with_motif
ATAC.motif.F.dwn$EnrichFold <- ATAC.motif.F.dwn$Percent_FG_with_motif/ATAC.motif.F.dwn$Percent_BG_with_motif

# Filter to retain only significantly enriched
sig.ATAC.motif.M.up  <- ATAC.motif.M.up  [ATAC.motif.M.up  $q.value < 0.05,]  # 335
sig.ATAC.motif.M.dwn <- ATAC.motif.M.dwn [ATAC.motif.M.dwn $q.value < 0.05,]  # 9
sig.ATAC.motif.F.up  <- ATAC.motif.F.up  [ATAC.motif.F.up  $q.value < 0.05,]  # 53
sig.ATAC.motif.F.dwn <- ATAC.motif.F.dwn [ATAC.motif.F.dwn $q.value < 0.05,]  # 322


# Put up and down together
sig.ATAC.motif.M                <- rbind(sig.ATAC.motif.M.up,sig.ATAC.motif.M.dwn)
sig.ATAC.motif.M$Direction      <- c(rep("UP",nrow(sig.ATAC.motif.M.up)),rep("DOWN",nrow(sig.ATAC.motif.M.dwn)))
sig.ATAC.motif.M$minusLog10pval <- -sig.ATAC.motif.M$Log_pvalue
sig.ATAC.motif.M$EnrichFold[-c(1:nrow(sig.ATAC.motif.M.up))] <-  -sig.ATAC.motif.M$EnrichFold[-c(1:nrow(sig.ATAC.motif.M.up))]   # make negative values for enrichment in regions of decreased accessibility

sig.ATAC.motif.F                <- rbind(sig.ATAC.motif.F.up,sig.ATAC.motif.F.dwn)
sig.ATAC.motif.F$Direction      <- c(rep("UP",nrow(sig.ATAC.motif.F.up)),rep("DOWN",nrow(sig.ATAC.motif.F.dwn)))
sig.ATAC.motif.F$minusLog10pval <- -sig.ATAC.motif.F$Log_pvalue
sig.ATAC.motif.F$EnrichFold[-c(1:nrow(sig.ATAC.motif.F.up))] <-  -sig.ATAC.motif.F$EnrichFold[-c(1:nrow(sig.ATAC.motif.F.up))]   # make negative values for enrichment in regions of decreased accessibility

write.table(sig.ATAC.motif.M, file = paste0(Sys.Date(),"_HOMER_Motifs_Enriched_in_Male_ONLY_ATAC_DA_peaks_FDR5.txt"), quote = F, sep = "\t", row.names = F)
write.table(sig.ATAC.motif.F, file = paste0(Sys.Date(),"_HOMER_Motifs_Enriched_in_Female_ONLY_ATAC_DA_peaks_FDR5.txt"), quote = F, sep = "\t", row.names = F)

# grab motifs related to transcription factor with RNAseq expression changes
M.age.TFs <- c("Stat4",
               "NFkB" , # Rela/Relb
               "STAT1",
               "IRF1" ,
               "Stat3")

F.age.TFs <- c("Mef2c"   ,
               "Meis1"   ,
               "HPC7-Scl",   # Tal1
               "IRF2"    ,
               "Nr2f2"   ,
               "TEAD1"   ,
               "Sox7"    )

# Jun in both up and down for F, excluded for simplicity


############################################################################################
# 2. Make bubble chart summary
library(ggplot2) 
library(scales) 
theme_set(theme_bw())

# get filtered/merged datafame for ggplot

sig.ATAC.motif.M.filt <- sig.ATAC.motif.M[grep(paste(M.age.TFs,collapse="|"), sig.ATAC.motif.M$Motif.Name),]
sig.ATAC.motif.F.filt <- sig.ATAC.motif.F[grep(paste(F.age.TFs,collapse="|"), sig.ATAC.motif.F$Motif.Name),]

M.srt <- sort(sig.ATAC.motif.M.filt$EnrichFold,decreasing = T, index.return = T)
F.srt <- sort(sig.ATAC.motif.F.filt$EnrichFold,decreasing = T, index.return = T)

filt.M.srt <- sig.ATAC.motif.M.filt[M.srt$ix,]
filt.F.srt <- sig.ATAC.motif.F.filt[F.srt$ix,]

filt.M.srt$Condition <- "M_ONLY_ATAC_DA"
filt.F.srt$Condition <- "F_ONLY_ATAC_DA"

#### Females
my.max <- max(filt.F.srt$EnrichFold)
my.min <- min(filt.F.srt$EnrichFold)
my.values <- c(my.min,0.75*my.min,0.5*my.min,0.25*my.min,0,0.25*my.max,0.5*my.max,0.75*my.max,my.max)
my.scaled <- rescale(my.values, to = c(0, 1))
my.color.vector <- c("darkblue","dodgerblue4","dodgerblue3","dodgerblue1","white","lightcoral","brown1","firebrick2","firebrick4")

# to preserve the wanted order
filt.F.srt$Motif.Name  <- factor(filt.F.srt$Motif.Name, levels = rev(unique(filt.F.srt$Motif.Name)))

pdf(paste0(Sys.Date(),"HOMER_Motifs_Enriched_ATAC_DA_FEMALE_ONLY_SelectedDETF_Aging_Perimac_FDR5.pdf"), height = 5, width=9)
my.plot <- ggplot(filt.F.srt,aes(x=Condition,y=Motif.Name,colour=EnrichFold,size=minusLog10pval))+ theme_bw()+ geom_point(shape = 16)
my.plot <- ggplot(filt.F.srt,aes(x=Condition,y=Motif.Name,colour=EnrichFold,size=minusLog10pval))+ theme(text = element_text(size=16))+ geom_point(shape = 16)
my.plot <- my.plot + ggtitle("HOMER") + labs(x = "PeriMac Aging", y = "Sex-Specific DA ATAC")
my.plot <- my.plot + scale_colour_gradientn(colours = my.color.vector,space = "Lab", na.value = "grey50", guide = "colourbar", values = my.scaled)
print(my.plot)
dev.off()  

#### Males
# need fake negative for color scale
filt.M.srt        <- rbind(filt.M.srt,filt.M.srt[1,])
filt.M.srt[7,1]   <- "NEG_CTL_FOR_PLOTTING"
filt.M.srt[7,3:5] <- 1
filt.M.srt[7,6:9] <- NA
filt.M.srt[7,10]  <- -0.5
filt.M.srt[7,11]  <- "DOWN"
filt.M.srt[7,12]  <- 5

my.max <- max(filt.M.srt$EnrichFold)
my.min <- min(filt.M.srt$EnrichFold)
my.values <- c(my.min,0.75*my.min,0.5*my.min,0.25*my.min,0,0.25*my.max,0.5*my.max,0.75*my.max,my.max)
my.scaled <- rescale(my.values, to = c(0, 1))
my.color.vector <- c("darkblue","dodgerblue4","dodgerblue3","dodgerblue1","white","lightcoral","brown1","firebrick2","firebrick4")

# to preserve the wanted order
filt.M.srt$Motif.Name  <- factor(filt.M.srt$Motif.Name, levels = rev(unique(filt.M.srt$Motif.Name)))

pdf(paste0(Sys.Date(),"HOMER_Motifs_Enriched_ATAC_DA_MALE_ONLY_SelectedDETF_Aging_Perimac_FDR5.pdf"), height = 4.5, width=9)
my.plot <- ggplot(filt.M.srt,aes(x=Condition,y=Motif.Name,colour=EnrichFold,size=minusLog10pval))+ theme_bw()+ geom_point(shape = 16)
my.plot <- ggplot(filt.M.srt,aes(x=Condition,y=Motif.Name,colour=EnrichFold,size=minusLog10pval))+ theme(text = element_text(size=16))+ geom_point(shape = 16)
my.plot <- my.plot + ggtitle("HOMER") + labs(x = "PeriMac Aging", y = "Sex-Specific DA ATAC")
my.plot <- my.plot + scale_colour_gradientn(colours = my.color.vector,space = "Lab", na.value = "grey50", guide = "colourbar", values = my.scaled)
print(my.plot)
dev.off()  


