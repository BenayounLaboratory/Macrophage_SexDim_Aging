setwd('/Volumes/BB_HQ_1//Immune_sex_dimorphism_Aging/Macrophages/Clean_Macrophage_Data/RNAseq/Functional_enrichment/TF_analysis/noSVA/')
options(stringsAsFactors = F)

# 2022-04-08
# check TF expression changing with age in a sex-specifc manner

# 2022-12-20
# also include the list of TF co-factors

#### load expression data
tissue.cts <- read.table('../../../DEseq2_analysis/DEseq2/2022-03-30_Peritoneal_Macrophages_DESeq2_model_with_AGE_SEX_log2_counts_matrix.txt', sep = "\t", header = T)

# get Sex-specific changes
fem.mal.age.comp <- read.table('../../../DEseq2_analysis/DEseq2/2022-03-30_Peritoneal_Macrophages_AGING_Separated_Sex_Merged_Table_ALL_genes_statistics.txt')
fem.spe          <- read.table('../../../DEseq2_analysis/DEseq2/2022-03-30_Peritoneal_Macrophages_AGING_Female.NOT.Male_FDR5_genes_statistics.txt')
mal.spe          <- read.table('../../../DEseq2_analysis/DEseq2/2022-03-30_Peritoneal_Macrophages_AGING_Male.NOT.Female_FDR5_genes_statistics.txt')

# read AnimalTFDB file
tfdb1 <- read.table('/Volumes/BB_HQ_1/Immune_sex_dimorphism_Aging/Macrophages/Macrophage_datasets_for_Comparison/AnimalTFDB3.0/2022-04-08_AnimalTFDB_3.0_Mus_musculus_TF.txt', sep = "\t", header = T)
tfdb2 <- read.table('/Volumes/BB_HQ_1/Immune_sex_dimorphism_Aging/Macrophages/Macrophage_datasets_for_Comparison/AnimalTFDB3.0/2022-04-15_AnimalTFDB_3.0_Mus_musculus_TF_cofactors.txt', sep = "\t", header = T)
# 1636 TFs, 970 co-factors

tfdb <- rbind(tfdb1[,-5],tfdb2)

########################################################################################################
# Extract TF genes that are changing with age in only one sex
# a. Females
fem.spe.tf <- intersect(rownames(fem.spe),tfdb$Symbol)
fem.spe.tf
# [1] "Eid1"    "Fam83g"  "Hoxb3"   "Hoxb5"   "Hoxb6"   "Hoxc4"   "Hoxc6"   "Hsf4"    "Irf2"    "Jdp2"    "Loxl2"   "Mef2c"   "Meis1"   "Nr1h3"   "Nr2f2"  
# [16] "Plag1"   "Pml"     "Rfx2"    "Sox7"    "Sp110"   "Tal1"    "Tbl1xr1" "Tead1"   "Tfec"    "Zfp697"   

pdf(paste0(Sys.Date(),"_Female_Specific_Aging_TF_genes_heatmap.pdf"))
pheatmap::pheatmap(tissue.cts[fem.spe.tf,], scale = "row", cluster_cols = F, 
                   colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
                   main = "Female specific TF age-regulation", cellheight = 12, cellwidth = 12)
dev.off()

write.table(fem.spe.tf, file = paste0(Sys.Date(),"_Female_Aging_Specific_TF_and_cofactors.txt"), sep = "\t", quote = F, col.names = F, row.names = F)

##### COMPARISON with Homer Known Motif Enrichment Results 
## Irf2 motif is enriched in ATAC up female specific (rank 1)                 IRF2(IRF)/Erythroblas-IRF2-ChIP-Seq(GSE36985)/Homer
## Jdp2-Jun dim; cJun motif is enriched in ATAC up female specific (rank 11)  c-Jun-CRE(bZIP)/K562-cJun-ChIP-Seq(GSE31477)/Homer
## Jdp2-Jun dim; JunB motif is enriched in ATAC up female specific (rank 23)  JunB(bZIP)/DendriticCells-Junb-ChIP-Seq(GSE36099)/Homer

## Meis1 motif is enriched in ATAC down female specific (rank 3)              Meis1(Homeobox)/MastCells-Meis1-ChIP-Seq(GSE48085)/Homer
## Mef2C motif is enriched in ATAC down female specific (rank 10)             Mef2c(MADS)/GM12878-Mef2c-ChIP-Seq(GSE32465)/Homer
## Tal1/SCL motif is enriched in ATAC down female specific (rank 53)          SCL(bHLH)/HPC7-Scl-ChIP-Seq(GSE13511)/Homer
## Nr2f2/COUP-TFII motif is enriched in ATAC down female specific (rank 125) 	COUP-TFII(NR)/Artia-Nr2f2-ChIP-Seq(GSE46497)/Homer
## Tead1 motif is enriched in ATAC down female specific (rank 261)            TEAD1(TEAD)/HepG2-TEAD1-ChIP-Seq(Encode)/Homer
## Sox7 motif is enriched in ATAC down female specific (rank 287)             Sox7(HMG)/ESC-Sox7-ChIP-Seq(GSE133899)/Homer


### Enrichments in H3K4me3 CUT and RUN
# Irf2, Meis1, Tead1 are up
# Scl both in up and down


# b. Males
mal.spe.tf <- intersect(rownames(mal.spe),tfdb$Symbol)
mal.spe.tf
# [1] "Stat4"   "Stat1"   "Sumo3"   "Tle2"    "Plek"    "Irf1"    "Ncor1"   "Stat3"   "Zfp493"  "Gtf2h2"  "Fam208a" "Ripk3"   "Clu"     "Rb1"    
# [15] "Ddx17"   "Paxbp1"  "Eno1b"   "Rela"    "Prdx5"   "Btaf1"   "Scai"    "Scand1"  "Znfx1"   "Txn1"    "Padi2"   "Relb"    "Rbpms"   "Calr"   
# [29] "Gtf2a2"  "Csrnp1"  "Klf8"   

pdf(paste0(Sys.Date(),"_Male_Specific_Aging_TF_genes_heatmap.pdf"))
pheatmap::pheatmap(tissue.cts[mal.spe.tf,], scale = "row", cluster_cols = F, 
                   colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
                   main = "Male specific TF age-regulation", cellheight = 12, cellwidth = 12)
dev.off()

write.table(mal.spe.tf, file = paste0(Sys.Date(),"_Male_Aging_Specific_TF_and_cofactors.txt"), sep = "\t", quote = F, col.names = F, row.names = F)

##### COMPARISON with Homer Known Motif Enrichment Results 
## Stat4 motif is enriched in ATAC up male specific (rank 5)
## Stat3 motif is enriched in ATAC up male specific (rank 24 and 68)
## Stat1 motif is enriched in ATAC up male specific (rank 144)
## NFKB motif is enriched in ATAC up male specific  (rank 149) [Rela/Relb]
## Irf1 motif is enriched in ATAC up male specific  (rank 270)

### No enrichments in H3K4me3 CUT and RUN


########################################################################################################
#### DEBUG
# my.gene    <- "Rps6ka3"
# my.cts     <- tissue.cts
# my.sig.res <- fem.mal.age.comp
# my.ylim    <- c(8,14)
# y.pval     <- 13

get_sex_age_bplots <- function(my.gene, my.cts, my.sig.res, my.ylim = c(0,20), y.pval = 19) {
  # extract data frame of expression
  my.gene.exp <- data.frame("YF" = as.numeric(my.cts[my.gene,1:5]  ) ,
                            "OF" = as.numeric(my.cts[my.gene,6:10]  ) ,
                            "YM" = as.numeric(my.cts[my.gene,11:15] ) ,
                            "OM" = as.numeric(my.cts[my.gene,16:20]) )
  
  pdfname <- paste0(Sys.Date(),"_",my.gene,"_Gene_Expression_boxplot_DEseq2.pdf")
  pdf(pdfname, height = 5, width = 3)
  boxplot(my.gene.exp,
          col = c("deeppink", "deeppink4","deepskyblue", "deepskyblue4"),
          ylab = "DESeq2 VST-normalized log2 counts",
          ylim = my.ylim, outline = F, main = my.gene)
  beeswarm::beeswarm(my.gene.exp, add = T, pch = 16, cex = 0.75)
  text(1.5, y.pval, signif(my.sig.res[my.gene,]$padj_F,2))
  text(3.5, y.pval, signif(my.sig.res[my.gene,]$padj_M,2))
  dev.off()
  
}


# ### Female specific TF with changes in ATAC accessibility
# get_sex_age_bplots("Irf2"    , tissue.cts, fem.mal.age.comp,  c(10, 16), 15)
# get_sex_age_bplots("Mef2c"   , tissue.cts, fem.mal.age.comp,  c(6, 12), 11)
# get_sex_age_bplots("Tead1"   , tissue.cts, fem.mal.age.comp,  c(5, 10), 9)
# get_sex_age_bplots("Nr2f2"   , tissue.cts, fem.mal.age.comp,  c(5, 10), 9)
# get_sex_age_bplots("Tal1"   , tissue.cts, fem.mal.age.comp,   c(5, 10), 9)
# get_sex_age_bplots("Meis1"   , tissue.cts, fem.mal.age.comp,  c(5, 10), 9)
# get_sex_age_bplots("Jdp2"   , tissue.cts, fem.mal.age.comp,   c(8, 12), 11)
# get_sex_age_bplots("Sox7"   , tissue.cts, fem.mal.age.comp,   c(6, 12), 11)
# 
# # tissue.cts["Sox7" ,]
# 
#######################
sink(file = paste(Sys.Date(),"_session_Info.txt", sep =""))
sessionInfo()
sink()




