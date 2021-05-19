slc = read_csv("data_sources/list/slc_final.csv")
abc = read_csv("data_sources/list/confirmed_abc.csv")
slc_pigment=read_csv("data_sources/cell_type_enriched/slc_lg_nsm_pigment_cells.csv")

additional_markers = tibble(GeneID = c("gcm","Six1"), Name = c("gcm","Six1"))
slcabc = bind_rows(slc, abc,additional_markers)


pigment = subset(urchin_2, features = slcabc$GeneID, idents="nsm_pigment_cells")
# pigment2 <- CreateSeuratObject(counts = pigment2)
pigment <- NormalizeData(object = pigment)
pigment <- FindVariableFeatures(object = pigment)
pigment <- ScaleData(object = pigment)
pigment <- RunPCA(object = pigment)
pigment <- FindNeighbors(object = pigment, dims = 1:30)
pigment <- FindClusters(object = pigment, resolution = 0.4)
pigment <- RunTSNE(object = pigment)
pigment_split_stages = SplitObject(pigment,split.by="orig.ident")

# DimPlot(object = pigment, reduction = "tsne")


lab_multi_DP(pigment_split_stages,dev_stage_names_list[-2],slc_pigment,"pigment_cell_subclusters",F,T)


# pdf("test_pigments.pdf",width=24, height = 16)
# lab_DP(pigment,slc_pigment,"Pigment cell subclusters and expression",T)
# dev.off()

# lab_DP_2 = function(scrna_df, feature_df, plot_title, flip_coord){

#   dp = DotPlot(scrna_df, features = feature_df$GeneID,split.by="orig.ident") + 
#     scale_x_discrete(breaks=c(feature_df$GeneID),
#                      labels=c(feature_df$Name)) +
#     RotatedAxis() +
#     ggtitle(as.character(plot_title))

#     if (flip_coord == T){
#       dp = DotPlot(scrna_df, features = feature_df$GeneID) + 
#       scale_x_discrete(breaks=c(feature_df$GeneID),
#                        labels=c(feature_df$Name)) +
#       RotatedAxis() +
#       ggtitle(as.character(plot_title)) +
#       coord_flip()
#   }
  
#   return(dp)
# }

plg = pigment_split_stages[7]
plg = plg$`Late Gastrula`
plg_count = plg@meta.data %>% count(seurat_clusters)

pigment_meta_cell_barplot = map(.x = pigment_split_stages, .f=function(x)(x@meta.data %>% 
																			count(seurat_clusters) %>% 
																			mutate(percentage=paste0(round(n/sum(n),digits=4)*100,"%")
																				)
																			)
)



meta_cell_barplot = function(meta_cell_df_list,title){
 plot = ggplot(meta_cell_df_list,mapping = aes(x=seurat_clusters,y=n)) +
		geom_bar(stat="identity") +
		geom_text(aes(label = percentage), vjust = 1.5, colour = "white") +
		xlab("Cluster ID") +
		ylab("Number of cells from a subcluster out of total in pigment cluster") +
		theme_classic() +
		ggtitle(title)

return(plot)

}

p1 = map2(.x=pigment_meta_cell_barplot, .y=names(pigment_split_stages), .f=function(x,y)(meta_cell_barplot(x,y)))

pdf("pigment_cell_subcluster_metadata.pdf",width=8,height=8)
p1
dev.off()