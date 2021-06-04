library(tidyverse)
library(Seurat)
library(plotly)
# library(cowplot)
# library(grid)
# library(tidyr)
# library(patchwork)
# library(gridExtra)
# library(parallel)

# extract_gene_locs finds which genes, from a user supplied list (features_list) of genes, are present in the 2020 Wessel single cell RNA seq dataset (scrna_df).
# features_list takes a list two column dataframe with the names GeneID and Name. The GeneID is the locus ID of a gene of interest.
# scrna_df is a SeuratObject.
# returns a two column tibble with columns 'GeneID' and 'Name'.
extract_gene_locs = function(features_list, scrna_df){
  
  feature_loc_list = rownames(scrna_df)[rownames(scrna_df) %in% features_list$GeneID]
  
  extracted_gene_locs = features_list[which(features_list$GeneID %in% feature_loc_list),] %>% 
  select(GeneID, Name) %>%
  unique()

  dup_check = duplicated(extracted_gene_locs$GeneID) %>% 
              sum()
  
  if(dup_check >= 1){

    extracted_gene_locs = extracted_gene_locs[-which(duplicated(extracted_gene_locs$GeneID)),]
    
    return(extracted_gene_locs)
    
    } else {
      
      return(extracted_gene_locs)
    
    }
}

# intramutate finds, within a dataframe (df), items in a list (search_replace_list$orig_search_list) which one seeks to replace items from another list (search_replace_list$replacement_list).
# returns a tibble with replaced elements
intramutate = function(df, orig_search_list, replacement_list){
  
  for (i in 1:length(replacement_list)){
    
    df = case_when(
      df == orig_search_list[i] ~ replacement_list[i],
      df != orig_search_list[i] ~ df
    )
  }
  
  return(df)
}

# produces a table summarizing metadata describing a SeuratObject
# returns a tibble 
dev_stage_cell_pop_stats = function(dev_stage_scrna){

  dev_stage_available_idents = levels(Idents(dev_stage_scrna))

  dev_stage_cell_ident_num_cell = map(.x = dev_stage_available_idents, 
                                      .f = function(x)(length(Cells(subset(dev_stage_scrna, idents = x))))) %>% 
                                  unlist()

  dev_stage_cell_ident_percentage = map(.x = dev_stage_available_idents,
                                        .f = function(x)(length(Cells(subset(dev_stage_scrna, idents = x)))/length(Cells(dev_stage_scrna)))) %>% 
                                    unlist()

  stats_df = tibble(
    dev_stage_available_idents, 
    dev_stage_cell_ident_num_cell, 
    dev_stage_cell_ident_percentage
          )

  return(stats_df)
}

lab_FP = function(scrna_df, feature_df){
  
  plot_list = FeaturePlot(scrna_df, features = feature_df$GeneID)
  
  for (i in 1:nrow(feature_df)){

    plot_list[[i]] = plot_list[[i]] + ggtitle(feature_df$Name[i])

  }

  return(plot_list)

}

lab_DP = function(scrna_df, feature_df, plot_title, flip_coord){

  dp = DotPlot(scrna_df, features = feature_df$GeneID) + 
    scale_x_discrete(breaks=c(feature_df$GeneID),
                     labels=c(feature_df$Name)) +
    RotatedAxis() +
    ggtitle(as.character(plot_title))

    if (flip_coord == T){
      dp = DotPlot(scrna_df, features = feature_df$GeneID) + 
      scale_x_discrete(breaks=c(feature_df$GeneID),
                       labels=c(feature_df$Name)) +
      RotatedAxis() +
      ggtitle(as.character(plot_title)) +
      coord_flip()
  }
  
  return(dp)

}

lab_multi_DP = function(scrna_df_list, title_list, lab_DP_feature_df, pdf_title, interactive = F, flip_coord = T){
  
  DP_list = map2(

    .x = scrna_df_list, 
    .y = title_list, 
    .f = function(x,y) (lab_DP(x, lab_DP_feature_df, y, flip_coord)

      )
    )
  
  if(interactive == F){

    pdf(paste0(pdf_title, ".pdf"), onefile=T, width=12, height=8)
    
    map(DP_list, print) %>% 
    invisible()
    
    dev.off()
    return(DP_list)  
  
  } else {

    dir.create(pdf_title)
    
    DP_list_interactive = map(DP_list, ggplotly) %>% invisible()
    
    map2(
      .x=DP_list_interactive, 
      .y = title_list, 
      .f = function(x,y)(htmlwidgets::saveWidget(as_widget(x), paste0(pdf_title,"/",y,".html"),selfcontained=TRUE))) %>% 
    invisible()
  }
}

find_markers_export = function(devolist, file_names_list){

  map2(.x = devolist,
     .y = file_names_list,
     .f = function(x,y)
      (
        FindAllMarkers(x,only.pos=T) %>% 
        write_csv(paste0("tmp.out/",y,".csv"))

        )
      )

}

cluster_wise_dea = function(scrna_df, avg_expression_df, feature_df, target_cluster, avg_exp_scaled_thresh, stage_name, export_plotted_loci = NULL, export_name = NULL){
  
  q1 = avg_expression_df[avg_expression_df$gene %in% feature_df$GeneID %>% which(),] %>% mutate(GeneID=gene)
  q2 = left_join(q1, feature_df, by = "GeneID")
  q3 = q2[c("GeneID","Name")]
  
  a = lab_DP(scrna_df, unique(q3), "eg", T)

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

  check = feature_df[g,]
  plot_df = mutate(check, Name = paste0(GeneID,"-",Name))
  qq = lab_DP(scrna_df,plot_df, paste0(stage_name,"-",target_cluster),T)

  return(qq)
  qq

}

# calculate_subclusters = function(){

#       # additional_markers = tibble(GeneID = c("gcm","Six1"), Name = c("gcm","Six1"))
      
      
#       # pigment = subset(urchin_2, features = slcabc$GeneID, idents="nsm_pigment_cells")
#       # # pigment2 <- CreateSeuratObject(counts = pigment2)
#       # pigment <- NormalizeData(object = pigment)
#       # pigment <- FindVariableFeatures(object = pigment)
#       # pigment <- ScaleData(object = pigment)
#       # pigment <- RunPCA(object = pigment)
#       # pigment <- FindNeighbors(object = pigment, dims = 1:30)
#       # pigment <- FindClusters(object = pigment, resolution = 0.4)
#       # pigment <- RunTSNE(object = pigment)
#       # pigment_split_stages = SplitObject(pigment,split.by="orig.ident")
      
#       # lab_multi_DP(pigment_split_stages,dev_stage_names_list[-2],slc_pigment,"pigment_cell_subclusters",F,T)
      
#       # plg = pigment_split_stages[7]
#       # plg = plg$`Late Gastrula`
#       # plg_count = plg@meta.data %>% count(seurat_clusters)
      
#       # pigment_meta_cell_barplot = map(.x = pigment_split_stages, .f=function(x)(x@meta.data %>% 
#       #                                       count(seurat_clusters) %>% 
#       #                                       mutate(percentage=paste0(round(n/sum(n),digits=4)*100,"%")
#       #                                         )
#       #                                       )
#       # )
      
      
      
#       # meta_cell_barplot = function(meta_cell_df_list,title){
#       #  plot = ggplot(meta_cell_df_list,mapping = aes(x=seurat_clusters,y=n)) +
#       #     geom_bar(stat="identity") +
#       #     geom_text(aes(label = percentage), vjust = 1.5, colour = "white") +
#       #     xlab("Cluster ID") +
#       #     ylab("Number of cells from a subcluster out of total in pigment cluster") +
#       #     theme_classic() +
#       #     ggtitle(title)
      
#       #   return(plot)
      
#       # }

#       # p1 = map2(.x=pigment_meta_cell_barplot, .y=names(pigment_split_stages), .f=function(x,y)(meta_cell_barplot(x,y)))
      
#       # pdf("pigment_cell_subcluster_metadata.pdf",width=8,height=8)
#       # p1
#       # dev.off()

# }