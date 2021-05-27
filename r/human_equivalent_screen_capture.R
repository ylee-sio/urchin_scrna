
df3 = left_join(df1,df2,by=Gene)
df3_tissue = df3[c(1,ncol(df3),str_which(names(df3),"Tissue"))]
df3_single = df3[c(1,ncol(df3),str_which(names(df3),"Single"))]

df6 = df3_single %>% gather(test,value,-Gene,-type) 
# map(.x = unique(df6$Gene), .f = function(x)(unique(df6 %>% subset(Gene == x) %>% arrange(desc(value)))))


df7 = map(.x = unique(df6$Gene), .f = function(x)(unique(df6 %>% subset(Gene == x) %>% arrange(desc(value)))))
df8 = map(df7,.f = function(x)((x$value)-mean(x$value)/sd(x$value)))
df9 = bind_rows(df7)
df10 = unlist(df8)
df11 = bind_cols(df9,zscore =df10) %>% na.omit()
df12 = df11 %>% group_by(Gene) %>% slice_max(order_by = zscore, n = 10) 
write_csv(df12, "top_10_humantissues_perSMT.csv")


temp_df11_zscore1 = case_when(df11$zscore<0~0, df11$zscore>=0~df11$zscore)
df11 = add_column(df11, temp_df11_zscore1) %>% mutate(zscore=temp_df11_zscore1)

df13 = df11 %>% group_by(Gene) %>% group_map(~ quantile(.x$zscore, probs = c(0.5,0.6, 0.7, 0.8, 0.9,.95, .96, .97, .98, .99), na.rm=T)) 
df14 = df13 %>% bind_rows() %>% add_column(Gene = unique(df12$Gene))
df15 = map2(.x=df14$Gene,.y=df14$`95%`,.f=function(x,y)(subset(df11, Gene == x & zscore >= y)))
df16 = map2(.x=df14$Gene,.y=df14$`90%`,.f=function(x,y)(subset(df11, Gene == x & zscore >= y))) %>% bind_rows()


df16 = mutate(df16, found_in_urchin_cell_type=type, urchin_cell_type = type ,found_in_human_cell_type = str_sub(test,14,-6),human_cell_type = test)

a = ggplot(df16) + geom_point(aes(x=Gene,y=zscore,found_in_human_cell_type = found_in_human_cell_type,found_in_urchin_cell_type=found_in_urchin_cell_type, col=human_cell_type)) + theme(axis.text.x=element_text(angle=90, hjust=1))
b = ggplotly(a, tooltip=c("Gene","found_in_urchin_cell_type","found_in_human_cell_type"))
htmlwidgets::saveWidget(as_widget(b), "top90p_of_z_byhumancell.html") 