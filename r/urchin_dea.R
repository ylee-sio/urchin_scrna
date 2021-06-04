ae_urchin_3_cell8 = FindAllMarkers(urchin_3_cell8,only.pos=T)
ae_urchin_3_cell64 = FindAllMarkers(urchin_3_cell64,only.pos=T)
ae_urchin_3_morula = FindAllMarkers(urchin_3_morula,only.pos=T)
ae_urchin_3_eb = FindAllMarkers(urchin_3_eb,only.pos=T)
ae_urchin_3_hb = FindAllMarkers(urchin_3_hb,only.pos=T)
ae_urchin_3_mb = FindAllMarkers(urchin_3_mb,only.pos=T)
ae_urchin_3_eg = FindAllMarkers(urchin_3_eg,only.pos=T)
s = FindAllMarkers(urchin_3_lg,only.pos=T)

#find in each, any slc or abc that might be present
#save data frames of present slcs and abcs from each differential expression df
#make dotplots using dfs
#somehow label/title each dotplots so that viewer is aware of the parameters used to make differential expression df

q1 = rownames(ae_urchin_3_cell8)[ae_urchin_3_cell8 %>% 
									rownames() %in% 
									slcabc$GeneID %>% 
									which()
									]

q1 = ae_urchin_3_eg[ae_urchin_3_eg$gene %in% slcabc$GeneID %>% which(),] %>% 
mutate(GeneID=gene)

q2 = left_join(q1, slcabc, by = "GeneID")
q3 = q2[c("GeneID","Name")]


a = lab_DP(urchin_3_eg, unique(q3), "eg", T)
# a = DotPlot(urchin_0, features = slc$GeneID,) + coord_flip()
b = ggplot_build(a)$plot$data
b = mutate(b, GeneID=rownames(b))
c = subset(b, avg.exp.scaled>=1.5)
c = arrange(c,desc(avg.exp.scaled))
cluster_name = "nsm_pigment_cells"
c = subset(c,id == cluster_name)
d = c$features.plot %>% unique() %>% droplevels()
d = tibble(GeneID=d,temp_name=d)
e = left_join(d,slcabc,by="GeneID")
f = unique(e$GeneID)
f %in% slcabc$GeneID %>% sum()
g = which(slcabc$GeneID %in% f)
slc2 = slcabc[g,]
slc3 = mutate(slc2, Name = paste0(GeneID,"-",Name))
qq = lab_DP(urchin_3_eg,slc3, paste0("EG-",cluster_name),T)
qq

find_cluster_diff = function(
	scrna_df, 
	avg_expression_df, 
	feature_df, 
	target_cluster, 
	avg_exp_scaled_thresh,
	stage_name
	)
{
	q1 = avg_expression_df[avg_expression_df$gene %in% feature_df$GeneID %>% which(),] %>% 
	mutate(GeneID=gene)
	
	q2 = left_join(q1, feature_df, by = "GeneID")
	q3 = q2[c("GeneID","Name")]
	
	
	a = lab_DP(scrna_df, unique(q3), "eg", T)
	# a = DotPlot(urchin_0, features = slc$GeneID,) + coord_flip()
	b = ggplot_build(a)$plot$data
	b = mutate(b, GeneID=rownames(b))
	c = subset(b, avg.exp.scaled>=avg_exp_scaled_thresh)
	c = arrange(c,desc(avg.exp.scaled))
	c = subset(c,id == target_cluster)
	d = c$features.plot %>% unique() %>% droplevels()
	d = tibble(GeneID=d,temp_name=d)
	e = left_join(d,feature_df,by="GeneID")
	f = unique(e$GeneID)
	f %in% feature_df$GeneID %>% sum()
	g = which(feature_df$GeneID %in% f)
	slc2 = feature_df[g,]
	slc3 = mutate(slc2, Name = paste0(GeneID,"-",Name))
	qq = lab_DP(scrna_df,slc3, paste0(stage_name,"-",target_cluster),T)

	return(qq)
	qq


}