slc = read_csv("data_sources/slc_list.csv")
nhr = read_csv("data_sources/nhr_list.csv")
cyp = read_csv("data_sources/cyp_list.csv")
cured = read_csv("data_sources/curated_list.csv")

prot_list = list(slc,nhr,cyp,cured)

genes_for_merging = list(c("ABCC5a", "LOC590074"))
merged_genes = AddModuleScore(object = urchin_1, features = genes_for_merging, name = "genes_for_merging")
FeaturePlot(object = merged_genes, features = "genes_for_merging1")

existing = map(.x = prot_list, .f = function(x)(extract_gene_locs(urchin_1,x)))