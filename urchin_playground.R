slc = read_csv("slc_list.csv")


source("~/scrna_r_tools.R", chdir = TRUE)
features_abc = read_csv("abc_list.csv") %>%
  select(GeneID, protein, subfamily, member)
features_slc = read_csv("slc_list.csv") %>%
  select(GeneID, protein, subfamily, member)
features_cyps = read_csv("cyps_list.csv") %>%
  select(GeneID, protein, subfamily, member) 
features_nhr = read_csv("nhr_list.csv") %>%
  select(GeneID, Name) %>% 
  drop_na()
features_final = read_csv("final_curation.csv")

abcb_0 = process_protein_df(urchin_1, features_abc, "abc", "ABC") %>% 
  subset(subfamily=="b") 
abcb = abcb_0[-8,]
abcc_0 = process_protein_df(urchin_1, features_abc, "abc", "ABC") %>% 
  filter(subfamily=="c") 
abcc = abcc_0[-c(6,17,19),]
slc_0 = process_protein_df(urchin_1, features_slc, "slc", "SLC")
slc = slc_0[-c(2,3,4,5,6,10,11,12,13,14,19),]
cyps = process_protein_df(urchin_1, features_cyps, "cyp", "CYP")
nhrs = process_protein_df(urchin_1, features_nhr, "nr", "nhr") %>% 
  transmute(GeneID = GeneID, named_locus = Name) %>% 
  drop_na()

c5_test1 = rownames(urchin_1) %in% c("LOC591982",
                          "LOC752233",
                          "LOC755624",
                          "LOC579421",
                          "LOC590074")
c5_test2 = rownames(urchin_1)[c5_test1]

c5_test3 = tibble(GeneID = c(c5_test2, "ABCC5a", "cell_typeA_score1","ABCC1", "ABCC9a"), named_locus = c("ABCC4a", "ABCC4d","ABCC4c","ABCC4b","ABCC5A_just_added", "ABCC5a","ABCC5a_combined", "ABCC1", "ABCC9a"))

# temp_efficient_plot(abcb, "abcb_x")
# temp_efficient_plot(abcc, "abcc_x")
# temp_efficient_plot(slc, "slc")
# temp_efficient_plot(cyps, "cyp")
# temp_efficient_plot(nhrs, "nhr")

labelled_dot_plotter(mb, slc)
labelled_dot_plotter(eg, abcc[7:10,])
a = labelled_feature_plotter(mb, slc[1:3,]) 
b = DimPlot(mb)
View(abcc_0)


al = map2(.x = devolist, .y = unique_orig_idents_replacements, .f = function(x,y) (DimPlot(x) + ggtitle(y)))
pdf("all_dev_stage_tsne.pdf", onefile = T, width=12, height = 8)
al
dev.off()

bl = map2(.x = devolist, .y = unique_orig_idents_replacements, .f = function(x,y) (labelled_dot_plotter(x,c5_test3) + ggtitle(y)))
pdf("abcc_dotplot_cleaned_test1.pdf", onefile = T, width = 16, height = 9)
bl
dev.off()

cl = map2(.x = devolist, .y = unique_orig_idents_replacements, .f = function(x,y) (labelled_feature_plotter(x,abcc[11:20,]) + plot_annotation(title = y)))
pdf("abcc_featureplot_2.pdf", onefile = T, width = 21, height = 12)
cl
dev.off()

genes_for_merging = list(c("ABCC5a", "LOC590074"))
merged_genes = AddModuleScore(object = urchin_1, features = genes_for_merging, name = "genes_for_merging")
FeaturePlot(object = merged_genes, features = "genes_for_merging1")

slc[which(slc$GeneID %in% get_gene_locs(u3_gcmMO,slc)),]
