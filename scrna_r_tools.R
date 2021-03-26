library(Seurat)
library(patchwork)
library(tidyverse)
library(gridExtra)
library(parallel)
library(plotly)
library(cowplot)
library(grid)
library(tidyr)

extract_gene_locs = function(scrna_df, features_list){
  
  feature_loc_list = rownames(scrna_df)[rownames(scrna_df) %in% features_list$GeneID]
  extracted_gene_locs = features_list[which(features_list$GeneID %in% feature_loc_list),] %>% unique()

  return(extracted_gene_locs)
}
intramutate = function(df, orig_search_list, replacement_list){
  
  for (i in 1:length(replacement_list)){
    
    df = case_when(
      df == orig_search_list[i] ~ replacement_list[i],
      df != orig_search_list[i] ~ df
    )
  }
  
  return(df)
}
lab_FP = function(scrna_df, feature_df){
  
  plot_list = FeaturePlot(scrna_df, features = feature_df$GeneID)
  
  for (i in 1:nrow(feature_df)){
    plot_list[[i]] = plot_list[[i]] + ggtitle(feature_df$Name[i])
  }
  return(plot_list)
}
lab_DP = function(scrna_df, feature_df, plot_title){
  
  dp = DotPlot(scrna_df, features = feature_df$GeneID) + 
    scale_x_discrete(breaks=c(feature_df$GeneID),
                     labels=c(feature_df$Name)) +
    RotatedAxis() +
    ggtitle(as.character(plot_title))
  
  return(dp)
}
lab_multi_DP = function(scrna_df_list, title_list, lab_DP_feature_df, pdf_title, interactive = F){
  
  DP_list = map2(.x = scrna_df_list, .y = title_list, .f = function(x,y) (lab_DP(x, lab_DP_feature_df, y)))
  
  if(interactive == F){
    pdf(paste0(pdf_title, ".pdf"), onefile=T, width=16, height=9)
    map(DP_list, print)
    dev.off()
    return(DP_list)  
  
  } else {
    DP_list_interactive = map(DP_list, ggplotly)
    map2(.x=DP_list_interactive, 
         .y = title_list, 
         .f = function(x,y)(htmlwidgets::saveWidget(as_widget(x), paste0(y,".html")))
         )
  }
}
