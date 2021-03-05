##data= data, y=cell line
diffmet<-function(data,y){
  
  f<-data %>% 
    filter(!is.na(SUPER_PATHWAY)) %>% 
    select(contains(y),BIOCHEMICAL) %>% 
    pivot_longer(cols=-c("BIOCHEMICAL"),names_to = "muestra") %>% 
    mutate(condicion=case_when(grepl("PARENTAL", muestra) ~ "Parental",
                               grepl("MESENCHYMAL",muestra) ~ "Mesenchymal" )) %>% 
    select(-muestra) %>% 
    group_by(BIOCHEMICAL,condicion) %>% 
    summarise(meanMet = mean(value)) %>% #group_by(BIOCHEMICAL) %>% 
    summarise(foldchange = log2(first(meanMet)/last(meanMet) ) ) ###meanMet[1]/meanMet[2]
  
  
  p<-data %>% 
    filter(!is.na(SUPER_PATHWAY)) %>% 
    select(contains(y),BIOCHEMICAL) %>% 
    group_by(BIOCHEMICAL) %>% 
    pivot_longer(cols=-c("BIOCHEMICAL"),names_to = "muestra") %>% 
    mutate(condicion=case_when(grepl("PARENTAL", muestra) ~ "Parental",
                               grepl("MESENCHYMAL",muestra) ~ "Mesenchymal" )) %>% 
    select(-muestra) %>% 
    t_test(value ~ condicion, data = .,paired=T ) %>% 
    adjust_pvalue(method = "fdr") %>%
    add_significance("p.adj") #%>% group_by(p.adj.signif) %>% count()
  
  
  z<-left_join(f,p) %>% mutate(threshold=case_when(
    foldchange > 1 & p.adj<.05 ~ 1, 
    foldchange < -1 & p.adj<.05 ~ 2,
    TRUE ~ 0))%>% 
    mutate(label=case_when(
      threshold>0 ~ BIOCHEMICAL, 
      TRUE ~ "")) %>%
    select(BIOCHEMICAL,foldchange,p.adj,p.adj.signif,threshold,label) %>% 
    mutate(across(threshold,as.factor))
  
  return(z)
}



