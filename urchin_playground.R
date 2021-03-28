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


dev_stage_cell_pop_stats = function(dev_stage_scrna,dev_stage){

  dev_stage_available_idents = levels(Idents(dev_stage_scrna))
  dev_stage_cell_ident_num_cell = map(.x = dev_stage_available_idents, .f = function(x)(length(Cells(subset(dev_stage_scrna, idents = x))))) %>% unlist()
  dev_stage_cell_ident_percentage = map(.x = dev_stage_available_idents, .f = function(x)(length(Cells(subset(dev_stage_scrna, idents = x)))/length(Cells(dev_stage_scrna)))) %>% unlist()

  stats_df = tibble(dev_stage,dev_stage_available_idents, dev_stage_cell_ident_num_cell, dev_stage_cell_ident_percentage)

  return(stats_df)
}

summaries_list = map2(.x = devolist,.y = dev_stage_names_abbreviated_list, .f = function(x,y)(dev_stage_cell_pop_stats(x,y)))









ggplot(summaries_list[[8]], aes(x="", y=dev_stage_cell_ident_percentage, fill=dev_stage_available_idents)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  
  theme_void()



fig <- plot_ly()
fig %>% add_pie(

	data = summaries_list[[3]],
	labels = ~dev_stage_available_idents,
	values = ~dev_stage_cell_ident_percentage,
	
	)
fig

Animals <- c("giraffes", "orangutans", "monkeys")
SF_Zoo <- c(20, 14, 23)
LA_Zoo <- c(12, 18, 29)
data <- data.frame(Animals, SF_Zoo, LA_Zoo)

fig <- plot_ly(summaries_list_export, x = ~dev_stage_available_idents, y = ~dev_stage_cell_ident_percentage, type = 'bar', name = 'SF Zoo')
fig <- fig %>% add_trace(y = ~LA_Zoo, name = 'LA Zoo')
fig <- fig %>% layout(yaxis = list(title = 'Count'), barmode = 'stack')

fig

ggplot(summaries_list_export, aes(x = dev_stage, y = dev_stage_cell_ident_percentage, fill = dev_stage_available_idents)) + 
    geom_bar(stat = "identity", width = 1) + 
    theme(legend.position = "none") +
    scale_x_discrete(NULL, expand = c(0,0)) +
    scale_y_continuous(NULL, expand = c(0,0)) + 
    coord_polar(theta = "y") +
    facet_wrap(~dev_stage)

