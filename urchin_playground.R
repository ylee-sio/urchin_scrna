slc = read_csv("slc_list.csv")
nhr = read_csv("nhr_list.csv")
cyps = read_csv("cyps_list.csv")
cured = read_csv("curated_list.csv")

genes_for_merging = list(c("ABCC5a", "LOC590074"))
merged_genes = AddModuleScore(object = urchin_1, features = genes_for_merging, name = "genes_for_merging")
FeaturePlot(object = merged_genes, features = "genes_for_merging1")

