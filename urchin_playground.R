slc_conMO = read_csv("data_sources/slc_list.csv") %>% 
	extract_gene_locs(u3_conMO)
slc_gcmMO = read_csv("data_sources/slc_list.csv") %>% 
	extract_gene_locs(u3_gcmMO)
nhr_conMO = read_csv("data_sources/nhr_list.csv") %>% 
	extract_gene_locs(u3_conMO)
nhr_gcmMO = read_csv("data_sources/nhr_list.csv") %>% 
	extract_gene_locs(u3_gcmMO)
cyp_conMO = read_csv("data_sources/cyp_list.csv") %>% 
	extract_gene_locs(u3_conMO)
cyp_gcmMO = read_csv("data_sources/cyp_list.csv") %>% 
	extract_gene_locs(u3_gcmMO)
cured_conMO = read_csv("data_sources/curated_list.csv") %>% 
	extract_gene_locs(u3_conMO)
cured_gcmMO = read_csv("data_sources/curated_list.csv") %>% 
	extract_gene_locs(u3_gcmMO)

conMO_list = list(slc_conMO,nhr_conMO, cyp_conMO, cured_conMO)
gcmMO_list = list(slc_gcmMO, nhr_gcmMO, cyp_gcmMO, cured_gcmMO)

protein_family_list = list("slc", "nhr", "cyp", "curated")


lab_DP(u3_conMO, conMO_list[[1]], paste0("48 HPF, Control Morpholino, Viewing expression for protein family: ",protein_family_list[[1]]))

t1 = map2(.x = conMO_list, 
	.y = protein_family_list,
 	.f = function(x,y)
 		(
 			lab_DP(
 				u3_conMO, 
 				x, 
 				paste0("48 HPF, Control Morpholino, Viewing expression for protein family: ",y)
 				)
 			)
 		)

t2 = map2(.x = gcmMO_list, 
	.y = protein_family_list,
 	.f = function(x,y)
 		(
 			lab_DP(
 				u3_gcmMO, 
 				x, 
 				paste0("48 HPF, GCM Morpholino, Viewing expression for protein family: ",y)
 				)
 			)
 		)
# prot_list = list(slc,nhr,cyp,cured)

# genes_for_merging = list(c("ABCC5a", "LOC590074"))
# merged_genes = AddModuleScore(object = urchin_1, features = genes_for_merging, name = "genes_for_merging")
# FeaturePlot(object = merged_genes, features = "genes_for_merging1")

# existing = map(.x = prot_list, .f = function(x)(extract_gene_locs(urchin_1,x)))

# #shows all cells of this identity
# germline = Cells(subset(urchin_1, idents = "germline"))

# #all cells in the dataset
# urchin_1_total_cells = Cells(urchin_1)
# t1 = map(devolist, .f = function(x)(length(Cells(x)))) %>% unlist()
# t2 = dev_stage_names_list
# t3 = tibble(dev_stage = t2, num_cells = t1, percentage_of_total = t1/length(urchin_1_total_cells))

# dev_stage_cell_pop_stats = function(dev_stage_scrna,dev_stage){

#   dev_stage_available_idents = levels(Idents(dev_stage_scrna))
#   dev_stage_cell_ident_num_cell = map(.x = dev_stage_available_idents, .f = function(x)(length(Cells(subset(dev_stage_scrna, idents = x))))) %>% unlist()
#   dev_stage_cell_ident_percentage = map(.x = dev_stage_available_idents, .f = function(x)(length(Cells(subset(dev_stage_scrna, idents = x)))/length(Cells(dev_stage_scrna)))) %>% unlist()

#   stats_df = tibble(dev_stage,dev_stage_available_idents, dev_stage_cell_ident_num_cell, dev_stage_cell_ident_percentage)

#   return(stats_df)
# }

#summaries_list = map2(.x = devolist,.y = dev_stage_names_abbreviated_list, .f = function(x,y)(dev_stage_cell_pop_stats(x,y)))









# ggplot(summaries_list[[8]], aes(x="", y=dev_stage_cell_ident_percentage, fill=dev_stage_available_idents)) +
#   geom_bar(stat="identity", width=1, color="white") +
#   coord_polar("y", start=0) +
  
#   theme_void()



# fig <- plot_ly()
# fig %>% add_pie(

# 	data = summaries_list[[3]],
# 	labels = ~dev_stage_available_idents,
# 	values = ~dev_stage_cell_ident_percentage,
	
# 	)
# fig

# Animals <- c("giraffes", "orangutans", "monkeys")
# SF_Zoo <- c(20, 14, 23)
# LA_Zoo <- c(12, 18, 29)
# data <- data.frame(Animals, SF_Zoo, LA_Zoo)

# fig <- plot_ly(summaries_list_export, x = ~dev_stage_available_idents, y = ~dev_stage_cell_ident_percentage, type = 'bar', name = 'SF Zoo')
# fig <- fig %>% add_trace(y = ~LA_Zoo, name = 'LA Zoo')
# fig <- fig %>% layout(yaxis = list(title = 'Count'), barmode = 'stack')

# fig

# ggplot(summaries_list_export, aes(x = dev_stage, y = dev_stage_cell_ident_percentage, fill = dev_stage_available_idents)) + 
#     geom_bar(stat = "identity", width = 1) + 
#     theme(legend.position = "none") +
#     scale_x_discrete(NULL, expand = c(0,0)) +
#     scale_y_continuous(NULL, expand = c(0,0)) + 
#     coord_polar(theta = "y") +
#     facet_wrap(~dev_stage)


t2 = unlist(json_data, use.names = F)
a = map(.x = str_extract_all(t2, "[:digit:]{1,}"), .f = function(x)(x[1])) %>% 
unlist() 
b = a[which(nchar(a) > 5)]
c = paste0("LOC",b) %>% 
as_tibble() %>% 
transmute(GeneID = value, Name = value)
d = extract_gene_locs(c, urchin_1)

l1.names = map(.x = json_data$children, .f = function(x)(x$name)) %>% unlist()
l1.l1.names = map(.x = json_data$children[[1]]$children, .f = function(x)(x$name)) %>% unlist()
l1.l2.names = map(.x = json_data$children[[2]]$children, .f = function(x)(x$name)) %>% unlist()
l1.l3.names = map(.x = json_data$children[[3]]$children, .f = function(x)(x$name)) %>% unlist()
l1.l4.names = map(.x = json_data$children[[4]]$children, .f = function(x)(x$name)) %>% unlist()
l1.l5.names = map(.x = json_data$children[[5]]$children, .f = function(x)(x$name)) %>% unlist()
l1.l6.names = map(.x = json_data$children[[6]]$children, .f = function(x)(x$name)) %>% unlist()

l1.children = map(.x = json_data$children[[1]]$children, .f = function(x)(x$children))
	l1.l1.children = map(.x = json_data$children[[1]]$children[[1]]$children, .f = function(x)(x$name))
		l1.l1.l1.l1.children = map(.x = json_data$children[[1]]$children[[1]]$children[[1]]$children, .f = function(x)(x$name)) %>% unlist()
		l1.l1.l1.l2.children = map(.x = json_data$children[[1]]$children[[1]]$children[[2]]$children, .f = function(x)(x$name)) %>% unlist()
		#####################
		l1.l1.l2.l1.children = map(.x = json_data$children[[1]]$children[[2]]$children[[1]]$children, .f = function(x)(x$name))
		l1.l1.l2.l2.children = map(.x = json_data$children[[1]]$children[[2]]$children[[2]]$children, .f = function(x)(x$name))
		l1.l1.l2.l3.children = map(.x = json_data$children[[1]]$children[[2]]$children[[3]]$children, .f = function(x)(x$name))
		l1.l1.l2.l4.children = map(.x = json_data$children[[1]]$children[[2]]$children[[4]]$children, .f = function(x)(x$name))
		l1.l1.l2.l5.children = map(.x = json_data$children[[1]]$children[[2]]$children[[5]]$children, .f = function(x)(x$name))
		######################
		l1.l1.l3.l1.children = map(.x = json_data$children[[1]]$children[[3]]$children[[1]]$children, .f = function(x)(x$name))
		l1.l1.l3.l2.children = map(.x = json_data$children[[1]]$children[[3]]$children[[2]]$children, .f = function(x)(x$name))
		l1.l1.l3.l3.children = map(.x = json_data$children[[1]]$children[[3]]$children[[3]]$children, .f = function(x)(x$name))
		l1.l1.l3.l4.children = map(.x = json_data$children[[1]]$children[[3]]$children[[4]]$children, .f = function(x)(x$name))
		######################
		l1.l1.l4.l1.children = map(.x = json_data$children[[1]]$children[[4]]$children[[1]]$children, .f = function(x)(x$name))
		l1.l1.l4.l2.children = map(.x = json_data$children[[1]]$children[[4]]$children[[2]]$children, .f = function(x)(x$name))
		l1.l1.l4.l3.children = map(.x = json_data$children[[1]]$children[[4]]$children[[3]]$children, .f = function(x)(x$name))
		l1.l1.l4.l4.children = map(.x = json_data$children[[1]]$children[[4]]$children[[4]]$children, .f = function(x)(x$name))
		l1.l1.l4.l5.children = map(.x = json_data$children[[1]]$children[[4]]$children[[5]]$children, .f = function(x)(x$name))
		######################
		l1.l1.l5.l1.children = map(.x = json_data$children[[1]]$children[[5]]$children[[1]]$children, .f = function(x)(x$name))
		l1.l1.l5.l2.children = map(.x = json_data$children[[1]]$children[[5]]$children[[2]]$children, .f = function(x)(x$name))
		l1.l1.l5.l3.children = map(.x = json_data$children[[1]]$children[[5]]$children[[3]]$children, .f = function(x)(x$name))
		######################
		l1.l1.l6.l1.children = map(.x = json_data$children[[1]]$children[[6]]$children[[1]]$children, .f = function(x)(x$name))
		l1.l1.l6.l2.children = map(.x = json_data$children[[1]]$children[[6]]$children[[2]]$children, .f = function(x)(x$name))
	


		a1 = tibble(l1 = l1.l1.names[[1]], l2 = l1.l1.children[[1]], l3 = l1.l1.l1.l1.children %>% unlist())
		a2 = tibble(l1 = l1.l1.names[[1]], l2 = l1.l1.children[[1]], l3 = l1.l1.l1.l2.children %>% unlist())
		b1 = tibble(l1 = l1.l1.names[[2]], l2 = l1.l1.children[[2]], l3 = l1.l1.l2.l1.children %>% unlist())
		b2 = tibble(l1 = l1.l1.names[[2]], l2 = l1.l1.children[[2]], l3 = l1.l1.l2.l2.children %>% unlist())
		b3 = tibble(l1 = l1.l1.names[[2]], l2 = l1.l1.children[[2]], l3 = l1.l1.l2.l3.children %>% unlist())
		b4 = tibble(l1 = l1.l1.names[[2]], l2 = l1.l1.children[[2]], l3 = l1.l1.l2.l4.children %>% unlist())
		b5 = tibble(l1 = l1.l1.names[[2]], l2 = l1.l1.children[[2]], l3 = l1.l1.l2.l5.children %>% unlist())
















	########################################################################################################################	
	l1.l2.children = map(.x = json_data$children[[2]]$children, .f = function(x)(x$children))
	########################################################################################################################
	l1.l3.children = map(.x = json_data$children[[3]]$children, .f = function(x)(x$children))
		l1.l3.1.children = map(.x = json_data$children[[3]]$children[[1]]$children, .f = function(x)(x$name))
		l1.l3.2.children = map(.x = json_data$children[[3]]$children[[2]]$children, .f = function(x)(x$name))
		l1.l3.3.children = map(.x = json_data$children[[3]]$children[[3]]$children, .f = function(x)(x$name))
		l1.l3.4.children = map(.x = json_data$children[[3]]$children[[4]]$children, .f = function(x)(x$name))
		l1.l3.5.children = map(.x = json_data$children[[3]]$children[[5]]$children, .f = function(x)(x$name))
		l1.l3.6.children = map(.x = json_data$children[[3]]$children[[6]]$children, .f = function(x)(x$name))
		l1.l3.7.children = map(.x = json_data$children[[3]]$children[[7]]$children, .f = function(x)(x$name))
		l1.l3.8.children = map(.x = json_data$children[[3]]$children[[8]]$children, .f = function(x)(x$name))
		l1.l3.9.children = map(.x = json_data$children[[3]]$children[[9]]$children, .f = function(x)(x$name))
		l1.l3.10.children = map(.x = json_data$children[[3]]$children[[10]]$children, .f = function(x)(x$name))
		l1.l3.11.children = map(.x = json_data$children[[3]]$children[[11]]$children, .f = function(x)(x$name))
		l1.l3.12.children = map(.x = json_data$children[[3]]$children[[12]]$children, .f = function(x)(x$name))
		l1.l3.13.children = map(.x = json_data$children[[3]]$children[[13]]$children, .f = function(x)(x$name))
		l1.l3.14.children = map(.x = json_data$children[[3]]$children[[14]]$children, .f = function(x)(x$name))
		l1.l3.15.children = map(.x = json_data$children[[3]]$children[[15]]$children, .f = function(x)(x$name))
		l1.l3.16.children = map(.x = json_data$children[[3]]$children[[16]]$children, .f = function(x)(x$name))
		l1.l3.17.children = map(.x = json_data$children[[3]]$children[[17]]$children, .f = function(x)(x$name))
		l1.l3.18.children = map(.x = json_data$children[[3]]$children[[18]]$children, .f = function(x)(x$name))
		l1.l3.19.children = map(.x = json_data$children[[3]]$children[[19]]$children, .f = function(x)(x$name))
		l1.l3.20.children = map(.x = json_data$children[[3]]$children[[20]]$children, .f = function(x)(x$name))
		l1.l3.21.children = map(.x = json_data$children[[3]]$children[[21]]$children, .f = function(x)(x$name))
		l1.l3.22.children = map(.x = json_data$children[[3]]$children[[22]]$children, .f = function(x)(x$name))
		l1.l3.23.children = map(.x = json_data$children[[3]]$children[[23]]$children, .f = function(x)(x$name))
		l1.l3.24.children = map(.x = json_data$children[[3]]$children[[24]]$children, .f = function(x)(x$name))
		l1.l3.25.children = map(.x = json_data$children[[3]]$children[[25]]$children, .f = function(x)(x$name))
		l1.l3.26.children = map(.x = json_data$children[[3]]$children[[26]]$children, .f = function(x)(x$name))
		l1.l3.27.children = map(.x = json_data$children[[3]]$children[[27]]$children, .f = function(x)(x$name))
		l1.l3.28.children = map(.x = json_data$children[[3]]$children[[28]]$children, .f = function(x)(x$name))
		l1.l3.29.children = map(.x = json_data$children[[3]]$children[[29]]$children, .f = function(x)(x$name))
		l1.l3.30.children = map(.x = json_data$children[[3]]$children[[30]]$children, .f = function(x)(x$name))
		l1.l3.31.children = map(.x = json_data$children[[3]]$children[[31]]$children, .f = function(x)(x$name))
		l1.l3.32.children = map(.x = json_data$children[[3]]$children[[32]]$children, .f = function(x)(x$name))
		l1.l3.33.children = map(.x = json_data$children[[3]]$children[[33]]$children, .f = function(x)(x$name))
		l1.l3.34.children = map(.x = json_data$children[[3]]$children[[34]]$children, .f = function(x)(x$name))
		l1.l3.35.children = map(.x = json_data$children[[3]]$children[[35]]$children, .f = function(x)(x$name))
		l1.l3.36.children = map(.x = json_data$children[[3]]$children[[36]]$children, .f = function(x)(x$name))
		l1.l3.37.children = map(.x = json_data$children[[3]]$children[[37]]$children, .f = function(x)(x$name))
		l1.l3.38.children = map(.x = json_data$children[[3]]$children[[38]]$children, .f = function(x)(x$name))
		l1.l3.39.children = map(.x = json_data$children[[3]]$children[[39]]$children, .f = function(x)(x$name))
		l1.l3.40.children = map(.x = json_data$children[[3]]$children[[40]]$children, .f = function(x)(x$name))
		l1.l3.41.children = map(.x = json_data$children[[3]]$children[[41]]$children, .f = function(x)(x$name))
		l1.l3.42.children = map(.x = json_data$children[[3]]$children[[42]]$children, .f = function(x)(x$name))
		l1.l3.43.children = map(.x = json_data$children[[3]]$children[[43]]$children, .f = function(x)(x$name))
		l1.l3.44.children = map(.x = json_data$children[[3]]$children[[44]]$children, .f = function(x)(x$name))
		l1.l3.45.children = map(.x = json_data$children[[3]]$children[[45]]$children, .f = function(x)(x$name))
		l1.l3.46.children = map(.x = json_data$children[[3]]$children[[46]]$children, .f = function(x)(x$name))
		l1.l3.47.children = map(.x = json_data$children[[3]]$children[[47]]$children, .f = function(x)(x$name))
		l1.l3.48.children = map(.x = json_data$children[[3]]$children[[48]]$children, .f = function(x)(x$name))
		l1.l3.49.children = map(.x = json_data$children[[3]]$children[[49]]$children, .f = function(x)(x$name))
		l1.l3.50.children = map(.x = json_data$children[[3]]$children[[50]]$children, .f = function(x)(x$name))
		l1.l3.51.children = map(.x = json_data$children[[3]]$children[[51]]$children, .f = function(x)(x$name))
		l1.l3.52.children = map(.x = json_data$children[[3]]$children[[52]]$children, .f = function(x)(x$name))
		l1.l3.53.children = map(.x = json_data$children[[3]]$children[[53]]$children, .f = function(x)(x$name))
		l1.l3.54.children = map(.x = json_data$children[[3]]$children[[54]]$children, .f = function(x)(x$name))
		l1.l3.55.children = map(.x = json_data$children[[3]]$children[[55]]$children, .f = function(x)(x$name))
		l1.l3.56.children = map(.x = json_data$children[[3]]$children[[56]]$children, .f = function(x)(x$name))
		l1.l3.57.children = map(.x = json_data$children[[3]]$children[[57]]$children, .f = function(x)(x$name))
		l1.l3.58.children = map(.x = json_data$children[[3]]$children[[58]]$children, .f = function(x)(x$name))
		l1.l3.59.children = map(.x = json_data$children[[3]]$children[[59]]$children, .f = function(x)(x$name))
		########################################################################################################################
	l1.l4.children = map(.x = json_data$children[[4]]$children, .f = function(x)(x$children))
		l1.l4.l1.children = map(.x = json_data$children[[4]]$children[[1]]$children, .f = function(x)(x$name))
			l1.l4.l1.l1.children = map(.x = json_data$children[[4]]$children[[1]]$children[[1]]$children, .f = function(x)(x$name))
			l1.l4.l1.l2.children = map(.x = json_data$children[[4]]$children[[1]]$children[[2]]$children, .f = function(x)(x$name))
			l1.l4.l1.l3.children = map(.x = json_data$children[[4]]$children[[1]]$children[[3]]$children, .f = function(x)(x$name))
			l1.l4.l1.l4.children = map(.x = json_data$children[[4]]$children[[1]]$children[[4]]$children, .f = function(x)(x$name))
			l1.l4.l1.l5.children = map(.x = json_data$children[[4]]$children[[1]]$children[[5]]$children, .f = function(x)(x$name))
			######################
		l1.l4.l2.children = map(.x = json_data$children[[4]]$children[[2]]$children, .f = function(x)(x$name))
			l1.l4.l2.l1.children = map(.x = json_data$children[[4]]$children[[2]]$children[[1]]$children, .f = function(x)(x$name))
			l1.l4.l2.l2.children = map(.x = json_data$children[[4]]$children[[2]]$children[[2]]$children, .f = function(x)(x$name))
			l1.l4.l2.l3.children = map(.x = json_data$children[[4]]$children[[2]]$children[[3]]$children, .f = function(x)(x$name))
			l1.l4.l2.l4.children = map(.x = json_data$children[[4]]$children[[2]]$children[[4]]$children, .f = function(x)(x$name))
			l1.l4.l2.l5.children = map(.x = json_data$children[[4]]$children[[2]]$children[[5]]$children, .f = function(x)(x$name))	
			l1.l4.l2.l6.children = map(.x = json_data$children[[4]]$children[[2]]$children[[6]]$children, .f = function(x)(x$name))
			l1.l4.l2.l7.children = map(.x = json_data$children[[4]]$children[[2]]$children[[7]]$children, .f = function(x)(x$name))
			l1.l4.l2.l8.children = map(.x = json_data$children[[4]]$children[[2]]$children[[8]]$children, .f = function(x)(x$name))
			l1.l4.l2.l9.children = map(.x = json_data$children[[4]]$children[[2]]$children[[9]]$children, .f = function(x)(x$name))
			l1.l4.l2.l10.children = map(.x = json_data$children[[4]]$children[[2]]$children[[10]]$children, .f = function(x)(x$name))
			l1.l4.l2.l11.children = map(.x = json_data$children[[4]]$children[[2]]$children[[11]]$children, .f = function(x)(x$name))
			######################
		l1.l4.l3.children = map(.x = json_data$children[[4]]$children[[3]]$children, .f = function(x)(x$name))
			l1.l4.l3.l1.children = map(.x = json_data$children[[4]]$children[[3]]$children[[1]]$children, .f = function(x)(x$name))
			l1.l4.l3.l2.children = map(.x = json_data$children[[4]]$children[[3]]$children[[2]]$children, .f = function(x)(x$name))
			######################
		l1.l4.l4.children = map(.x = json_data$children[[4]]$children[[4]]$children, .f = function(x)(x$name))
		######################
		l1.l4.l5.children = map(.x = json_data$children[[4]]$children[[5]]$children, .f = function(x)(x$name))
			l1.l4.l5.l1.children = map(.x = json_data$children[[4]]$children[[5]]$children[[1]]$children, .f = function(x)(x$name))
			l1.l4.l5.l2.children = map(.x = json_data$children[[4]]$children[[5]]$children[[2]]$children, .f = function(x)(x$name))
		######################
		l1.l4.l6.children = map(.x = json_data$children[[4]]$children[[6]]$children, .f = function(x)(x$name))
			l1.l4.l6.l1.children = map(.x = json_data$children[[4]]$children[[6]]$children[[1]]$children, .f = function(x)(x$name))
			l1.l4.l6.l2.children = map(.x = json_data$children[[4]]$children[[6]]$children[[2]]$children, .f = function(x)(x$name))
			l1.l4.l6.l3.children = map(.x = json_data$children[[4]]$children[[6]]$children[[3]]$children, .f = function(x)(x$name))
			l1.l4.l6.l4.children = map(.x = json_data$children[[4]]$children[[6]]$children[[4]]$children, .f = function(x)(x$name))
			l1.l4.l6.l5.children = map(.x = json_data$children[[4]]$children[[6]]$children[[5]]$children, .f = function(x)(x$name))
		######################
		l1.l4.l7.children = map(.x = json_data$children[[4]]$children[[7]]$children, .f = function(x)(x$name))
			l1.l4.l7.l1.children = map(.x = json_data$children[[4]]$children[[7]]$children[[1]]$children, .f = function(x)(x$name))
			l1.l4.l7.l2.children = map(.x = json_data$children[[4]]$children[[7]]$children[[2]]$children, .f = function(x)(x$name))
			l1.l4.l7.l3.children = map(.x = json_data$children[[4]]$children[[7]]$children[[3]]$children, .f = function(x)(x$name))
		######################
		l1.l4.l8.children = map(.x = json_data$children[[4]]$children[[8]]$children, .f = function(x)(x$name))
			l1.l4.l8.l1.children = map(.x = json_data$children[[4]]$children[[8]]$children[[1]]$children, .f = function(x)(x$name))
			l1.l4.l8.l2.children = map(.x = json_data$children[[4]]$children[[8]]$children[[2]]$children, .f = function(x)(x$name))
			l1.l4.l8.l3.children = map(.x = json_data$children[[4]]$children[[8]]$children[[3]]$children, .f = function(x)(x$name))
			l1.l4.l8.l4.children = map(.x = json_data$children[[4]]$children[[8]]$children[[4]]$children, .f = function(x)(x$name))
			l1.l4.l8.l5.children = map(.x = json_data$children[[4]]$children[[8]]$children[[5]]$children, .f = function(x)(x$name))
			l1.l4.l8.l6.children = map(.x = json_data$children[[4]]$children[[8]]$children[[6]]$children, .f = function(x)(x$name))
		l1.l4.l9.children = map(.x = json_data$children[[4]]$children[[9]]$children, .f = function(x)(x$name))
			l1.l4.l9.l1.children = map(.x = json_data$children[[4]]$children[[9]]$children[[1]]$children, .f = function(x)(x$name))
			l1.l4.l9.l2.children = map(.x = json_data$children[[4]]$children[[9]]$children[[2]]$children, .f = function(x)(x$name))
			l1.l4.l9.l3.children = map(.x = json_data$children[[4]]$children[[9]]$children[[3]]$children, .f = function(x)(x$name))
			l1.l4.l9.l4.children = map(.x = json_data$children[[4]]$children[[9]]$children[[4]]$children, .f = function(x)(x$name))








l1.l4.children = map(.x = json_data$children[[4]]$children, .f = function(x)(x$children))
l1.l5.children = map(.x = json_data$children[[5]]$children, .f = function(x)(x$children))
l1.l6.children = map(.x = json_data$children[[6]]$children, .f = function(x)(x$children))








a1 = tibble(l1 = l1.l1.names[[1]], l2 = l1.l1.children[[1]], l3 = l1.l1.l1.l1.children)
a2 = tibble(l1 = l1.l1.names[[1]], l2 = l1.l1.children[[1]], l3 = l1.l1.l1.l2.children)













































































