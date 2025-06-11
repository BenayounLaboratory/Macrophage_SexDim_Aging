setwd('/Volumes/BB_HQ_1/Immune_sex_dimorphism_Aging/Macrophages/Clean_Macrophage_Data/Intervention_Data/shRNA_analyses/Comparison')
options(stringsAsFactors = F)

library(bitops)
library(UpSetR)
# # devtools::install_github("yanlinlin82/ggvenn")
library(ggvenn)
library(ggplot2)

library(clusterProfiler)

# 2025-01-30
# Compare J774 and RAW sh IRF2 RNAseq
### Relax FDR to 1% since we are only taking those that overlap

################################################################################################################
##### 1. Read DESeq2 results
j774.sh.all <- read.table('../J774_cells/DESeq2/2025-01-30_J774_sh_Irf2_data_DESeq2_shIrf2_DEseq2_results_COMBINED_all_genes_statistics.txt' , header = T)
raw.sh.all  <- read.table('../RAW267_cells/DEseq2/2024-08-16_RAW_sh_Irf2_data_DESeq2_shIrf2_DEseq2_results_COMBINED_all_genes_statistics.txt', header = T)

j774.sh.fdr1 <- j774.sh.all [j774.sh.all$p_combined < 0.01,] ## 3437
raw.sh.fdr1  <- raw.sh.all  [raw.sh.all $p_combined < 0.01,] ## 5209
################################################################################################################

################################################################################################################
##### 2. prepare upset
listInput <- list("RAW264.7_UP"    = raw.sh.fdr1$Row.names[raw.sh.fdr1$log2FoldChange.sh1 >0],
                  "RAW264.7_DWN"   = raw.sh.fdr1$Row.names[raw.sh.fdr1$log2FoldChange.sh1 <0],
                  "J774A.1_UP"     = j774.sh.fdr1$Row.names[j774.sh.fdr1$log2FoldChange.sh1 >0],
                  "J774A.1_DWN"    = j774.sh.fdr1$Row.names[j774.sh.fdr1$log2FoldChange.sh1 <0])

inputData <- fromList(listInput)

pdf(paste0(Sys.Date(),"_J774_RAW_shIrf2_Upset_FDR1.pdf"), height = 4, width = 5)
upset(inputData, 
      order.by = "freq", 
      keep.order= T,
      nsets = length(listInput),       # Ensure all sets are displayed
      nintersects = NA,                # Show all possible intersections
      main.bar.color = "peru", 
      point.size = 3, 
      line.size = 1)
dev.off()
################################################################################################################

################################################################################################################
##### 3. Venn diagrams

# genes expressed robustly in either dataset
universe <- unique(j774.sh.all$Row.names, raw.sh.all$Row.names  )
length(universe) # 14281

### UP/UP
my.ups <- list("RAW264.7_UP"    = raw.sh.fdr1$Row.names[raw.sh.fdr1$log2FoldChange.sh1 >0],
               "J774A.1_UP"     = j774.sh.fdr1$Row.names[j774.sh.fdr1$log2FoldChange.sh1 >0])

overlap.up.test <- fisher.test(matrix(c(601,1910,958,14281-601-1910-958),2,2))
overlap.up.test$p.value ### 5.571926e-98

pdf(paste0(Sys.Date(),"_J774_RAW_shIrf2_UP_Venn_FDR1.pdf"))
ggvenn(my.ups, 
       fill_color = c("seagreen1", "seagreen"),
       stroke_size = 0.5, set_name_size = 4,
       auto_scale = T, show_stats = "c") + ggtitle(paste0("p = ",signif(overlap.up.test$p.value,3)))
dev.off()


consistent.up <- intersect(raw.sh.fdr1$Row.names[raw.sh.fdr1$log2FoldChange.sh1 >0],j774.sh.fdr1$Row.names[j774.sh.fdr1$log2FoldChange.sh1 >0]) # 180

write.table(consistent.up, file = paste0(Sys.Date(),"_J774_RAW_shIrf2_UP_consistently_FDR1.txt"), quote = F, col.names = F, row.names = F)

### DOWN/DOWN
my.dwns <- list("RAW264.7_DWN"   = raw.sh.fdr1$Row.names[raw.sh.fdr1$log2FoldChange.sh1 <0],
                "J774A.1_DWN"    = j774.sh.fdr1$Row.names[j774.sh.fdr1$log2FoldChange.sh1 <0])

overlap.dwn.test <- fisher.test(matrix(c(599,2099,1279,14281-599-2099-1279),2,2))
overlap.dwn.test$p.value ### 3.642008e-48

pdf(paste0(Sys.Date(),"_J774_RAW_shIrf2_DWN_Venn_FDR1.pdf"))
ggvenn(my.dwns, 
       fill_color = c("seagreen1", "seagreen"),
       stroke_size = 0.5, set_name_size = 4,
       auto_scale = T, show_stats = "c") + ggtitle(paste0("p = ",signif(overlap.dwn.test$p.value,3)))
dev.off()

consistent.dwn <- intersect(raw.sh.fdr1$Row.names[raw.sh.fdr1$log2FoldChange.sh1 <0],j774.sh.fdr1$Row.names[j774.sh.fdr1$log2FoldChange.sh1 <0]) # 131

write.table(consistent.dwn, file = paste0(Sys.Date(),"_J774_RAW_shIrf2_DWN_consistently_FDR1.txt"), quote = F, col.names = F, row.names = F)
################################################################################################################


################################################################################################################
##### 4. Functional enrichment

source('ORA_support_functions.R')

### Load genesets of interest
Sym.m2.reactome        <- read.gmt("/Volumes/BB_HQ_1/PATHWAY_ANNOT/MSigDB/m2.cp.reactome.v2022.1.Mm.symbols.gmt"  )
Sym.ENS.GO.ALL         <- read.gmt("/Volumes/BB_HQ_1/PATHWAY_ANNOT/ENSEMBL/2022-12-21_mouse_Ens108_GO_ALL.gmt"    ) 

### Run ORA
cons.up.react     <-  run_enrich(consistent.up   , universe, "shIRF2_consistent_UP ", Sym.m2.reactome, "MSigDB_Reactome")
cons.dwn.react    <-  run_enrich(consistent.dwn  , universe, "shIRF2_consistent_DWN", Sym.m2.reactome, "MSigDB_Reactome")

cons.up.go        <-  run_enrich(consistent.up   , universe, "shIRF2_consistent_UP ", Sym.ENS.GO.ALL, "ENSEMBL_GO")
cons.dwn.go       <-  run_enrich(consistent.dwn  , universe, "shIRF2_consistent_DWN", Sym.ENS.GO.ALL, "ENSEMBL_GO")


# plot top genesets
top_ORA_sum(cons.up.go, cons.dwn.go, "shIRF2_GO")
top_ORA_sum(cons.up.react, cons.dwn.react, "shIRF2_Reactome")

# export results
options(java.parameters = "-Xmx16g" )
require(openxlsx)

ora.res <- list("sh_UP_REACTOME"  = cons.up.react @result [cons.up.react @result$p.adjust < 0.05,],
                "sh_DWN_REACTOME" = cons.dwn.react@result [cons.dwn.react@result$p.adjust < 0.05,],
                "sh_UP_GO_ALL"    = cons.up.go  @result [cons.up.go  @result$p.adjust < 0.05,],
                "sh_DWN_GO_ALL"   = cons.dwn.go @result [cons.dwn.go @result$p.adjust < 0.05,])
write.xlsx(ora.res, rowNames = F, file = paste0(Sys.Date(),"_shIRF2_consistent_ClusterProfiler_ORA_REACTOME_GO_Results.xlsx"))


#######################
sink(file = paste(Sys.Date(),"shIrf2_Comparison_RNAseq_session_Info.txt", sep =""))
sessionInfo()
sink()



