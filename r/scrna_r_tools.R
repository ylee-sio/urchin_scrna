library(Seurat)
library(patchwork)
library(tidyverse)
library(gridExtra)
library(parallel)
library(plotly)
library(cowplot)
library(grid)
library(tidyr)

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
intramutate = function(df, search_replace_list){
  
  for (i in 1:length(search_replace_list$replacement_list)){
    
    df = case_when(
      df == search_replace_list$orig_search_list[i] ~ search_replace_list$replacement_list[i],
      df != search_replace_list$orig_search_list[i] ~ df
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

# resub = function(cbmc, neighbors_dim, tsne_dim, cluster_res){

# cbmc <- NormalizeData(cbmc)
# cbmc <- FindVariableFeatures(cbmc)
# cbmc <- ScaleData(cbmc)
# cbmc <- RunPCA(cbmc, verbose = FALSE)
# cbmc <- FindNeighbors(cbmc, dims = neighbors_dim)
# cbmc <- FindClusters(cbmc, resolution = as.numeric(cluster_res), verbose = FALSE)
# cbmc <- RunTSNE(cbmc, dims = tsne_dim)

# return(cbmc)

# }

# lg2 <- NormalizeData(lg2)
# lg2 <- FindVariableFeatures(lg2)
# lg2 <- ScaleData(lg2)
# lg2 <- FindVariableFeatures(lg2)
# lg2 <- RunPCA(lg2, verbose = FALSE)
# lg2 <- FindNeighbors(lg2, dims = 1:30)
# lg2 <- FindClusters(lg2, resolution = as.numeric(4.0), verbose = FALSE)
# lg2 <- RunTSNE(lg2, dims = 1:30)


# t9 = read_csv("data_sources/list/slc_final.csv")
# t0 = read_csv("data_sources/list/slc_final.csv")
# t0 = bind_rows(t9,t0[which(!(t0$GeneID %in% t9$GeneID)),])

# # t1 = read_csv("~/projects/purp_scrna/data_sources/kegg_filtered/MPL.csv")
# # t2 = t1[str_which(t1$Name, "nhr"),]
# # t3 = str_extract_all(t2$Name, boundary("word")) %>% 
# # str_extract_all("SLC[:alpha:]*[:digit:]*[:alpha:]*[:digit:]*", simplify = T)
# # t4 = t3[,1]
# # t5 = mutate(t2, Name = t4)

# # pdf("test.pdf", height = 48, width = 24)
# # DotPlot(lg, features = t5$GeneID) + 
# #     scale_x_discrete(breaks=c(t5$GeneID),
# #                      labels=c(t5$Name)) +
# # coord_flip() + 
# # RotatedAxis()
# # dev.off()

# t1 = FindAllMarkers(lg2,features=t0$GeneID[1:440],logfc.threshold=.5, only.pos=T)
# t1 = group_by(t1, cluster) %>% group_map(.f=function(x,...)(subset(x,p_val==0)),.keep=T)


a = DotPlot(lg, features = slc$GeneID,) + coord_flip()
b = ggplot_build(a)$plot$data
b = mutate(b, GeneID=rownames(b))
c = subset(b, avg.exp.scaled<=.5 & pct.exp <= 1)
d = merge(c, slc, by = "GeneID")
# e = d %>% na.omit()
e = d
f = unique(e$GeneID)
f %in% slc$GeneID %>% sum()
g = which(f %in% slc$GeneID)
slc2 = slc[-g,]
slc3 = mutate(slc2, Name = paste0(GeneID,"-",Name))
qq = lab_DP(lg,slc3, "meep",T)

pdf("lg_extraction_list.pdf", width = 8, height = 48)
qq
dev.off()



