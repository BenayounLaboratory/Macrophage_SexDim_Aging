setwd('/Volumes/BB_HQ_1/PATHWAY_ANNOT/X_inactivation_escape/')
options(stringsAsFactors = F)

library(readxl)

my.data <- data.frame(read_xlsx('2023-02-24_XCI_escapees.xlsx', col_names = F))


my.xci2023 <- as.vector(na.omit(unique(unlist(my.data[,-1]))))

write.table(my.xci2023, file = "2023-11-01_updated_known_mouse_XCI_genes.txt", sep = "\r", col.names = F, row.names = F, quote = F)
