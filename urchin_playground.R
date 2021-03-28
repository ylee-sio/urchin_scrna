slc = read_csv("data_sources/slc_list.csv")
nhr = read_csv("data_sources/nhr_list.csv")
cyp = read_csv("data_sources/cyp_list.csv")
cured = read_csv("data_sources/curated_list.csv")

prot_list = list(slc,nhr,cyp,cured)

genes_for_merging = list(c("ABCC5a", "LOC590074"))
merged_genes = AddModuleScore(object = urchin_1, features = genes_for_merging, name = "genes_for_merging")
FeaturePlot(object = merged_genes, features = "genes_for_merging1")

existing = map(.x = prot_list, .f = function(x)(extract_gene_locs(urchin_1,x)))

#shows all cells of this identity
germline = Cells(subset(urchin_1, idents = "germline"))

#all cells in the dataset
urchin_1_total_cells = Cells(urchin_1)
t1 = map(devolist, .f = function(x)(length(Cells(x)))) %>% unlist()
t2 = dev_stage_names_list
t3 = tibble(dev_stage = t2, num_cells = t1, percentage_of_total = t1/length(urchin_1_total_cells))


dev_stage_cell_pop_stats = function(dev_stage_scrna){

  dev_stage_available_idents = levels(Idents(dev_stage_scrna))
  dev_stage_cell_ident_num_cell = map(.x = dev_stage_available_idents, .f = function(x)(length(Cells(subset(dev_stage_scrna, idents = x))))) %>% unlist()
  dev_stage_cell_ident_percentage = map(.x = dev_stage_available_idents, .f = function(x)(length(Cells(subset(dev_stage_scrna, idents = x)))/length(Cells(dev_stage_scrna)))) %>% unlist()

  stats_df = tibble(dev_stage_available_idents, dev_stage_cell_ident_num_cell, dev_stage_cell_ident_percentage)

  return(stats_df)
}