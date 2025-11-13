README

### 2023-10-30

# get significant genes only
cat 2022-03-30_Peritoneal_Macrophages_AGING_Female_ALL_genes_statistics.txt | perl -lane 'next if $_=~m/^base/; next if $F[6]=~m/NA/; print $_ if ($F[6] < 0.05);' > 2022-03-30_Peritoneal_Macrophages_AGING_Female_FDR5_genes_statistics.txt
cat 2022-03-30_Peritoneal_Macrophages_AGING_Male_ALL_genes_statistics.txt | perl -lane 'next if $_=~m/^base/; next if $F[6]=~m/NA/; print $_ if ($F[6] < 0.05);'   > 2022-03-30_Peritoneal_Macrophages_AGING_Male_FDR5_genes_statistics.txt

##     355 2022-03-30_Peritoneal_Macrophages_AGING_Female_FDR5_genes_statistics.txt
##     358 2022-03-30_Peritoneal_Macrophages_AGING_Male_FDR5_genes_statistics.txt
  
# Format for NetworkAnalyst
cat 2022-03-30_Peritoneal_Macrophages_AGING_Male_FDR5_genes_statistics.txt   | cut -f 1,3 > PeriMac_Aging_FDR5_MaleONLY_for_NetworkAnalyst.txt
cat 2022-03-30_Peritoneal_Macrophages_AGING_Female_FDR5_genes_statistics.txt | cut -f 1,3 > PeriMac_Aging_FDR5_FemaleONLY_for_NetworkAnalyst.txt
     
# NetworkAnalyst 3.0
   - FDR5 UP and DOWN sex genes
   - Official Gene Symbol
   - LogFC
   - List Enrichment network
   - Generic PPI, IMEx (InnateDB)
   - Background: White
   - Layout: Fruchterman-Reingold
   - View: Expression: color Blind 
		
### Male Aging, largest subnetwork (#1)
	Networks	Nodes	Edges	Seeds
	subnetwork1	832	987	98
	
### Female Aging, largest subnetwork (#1)
	Networks	Nodes	Edges	Seeds
	subnetwork1	868	955	88
	