---
title: "Metabolomic analysis"
author: "Meztli Matadamas"
date: "20/10/2020"
output:
  html_document: default
  pdf_document: default
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)

library(Rtsne)
library(RColorBrewer) 

library(tidyverse)
library(rstatix)
library(cowplot)
```

# Epithelial to Mesenchymal transition
Epithelial-to-mesenchymal transition (EMT) relates to many molecular and cellular alterations that occur when epithelial cells undergo a switch in differentiation generating mesenchymal-like cells with newly acquired migratory and invasive properties. In cancer cells, EMT leads to drug resistance and metastasis. Moreover, differences in genetic backgrounds, even between patients with the same type of cancer, also determine resistance to some treatments. Metabolic rewiring is essential to induce EMT, hence it is important to identify key metabolic elements for this process, which can be later used to treat cancer cells with different genetic backgrounds. 

This protocol analyze Metabolomic data from [Sun Y., *et. al.* 2014](https://doi.org/10.1186/2049-3002-2-20) to infer reactions changing after EMT induction with TGF-$\beta$.  

The dataset contains data from three non-small cell lung cancer cell lines: **A549, HCC827 and NCI-H358**

These cell lines have different genetic backgrounds, hence we aimed to identify the metabolic similarities among them.
A549 and NCI-H358 have mutations in KRAS, while HCC827 has mutations in EGFR. To find key EMT features, we compared the three cell lines data and found their similarities. 

# Reading the data

```{r import,message=FALSE}
EMT <- vroom::vroom("input/Metabolomica_antes-despues_EMT.txt")

EMT <- janitor::clean_names(dat = EMT, 
                            case = "screaming_snake")
EMT<- EMT %>% filter(!is.na(BIOCHEMICAL))
EMT
```

# PCA and tSNE generation

In order to explore the data we performed two dimentional reduction algorithms: PCA and tSNE. 

```{r pca, echo=T, message=FALSE}
PCAMatrix<-EMT %>% 
  select(contains(c("A549","HCC","NCI"))) %>% drop_na()

pca = prcomp(t(PCAMatrix))
nametable<-as.data.frame(cbind(colnames(PCAMatrix),gsub("\\_[1-9]","",colnames(PCAMatrix))))
nametable$V2<-as.factor(nametable$V2)
nametable$V2<-factor(nametable$V2,levels = levels(nametable$V2)[c(2,1,4,3,6,5)])

library(ggfortify)
aplot<-autoplot(pca,data=nametable,colour='V2',frame=F,size=3)+
  scale_color_brewer(palette = "Paired")+
  theme(text = element_text(size=12),legend.position="bottom",legend.title = element_blank())

set.seed(80)
results=Rtsne(t(PCAMatrix),dims=2,perplexity=5,verbose=F,max_iter=10000)
# plotting
bplot<-ggplot(as.data.frame(results$Y),aes(x=results$Y[,1],y=results$Y[,2],col=nametable$V2))+
  xlab("tSNE1") + ylab("tSNE2")+
  geom_point(size=3)+ #Size and alpha just for fun
  scale_color_brewer(palette = "Paired")+ #your colors here
  theme(text = element_text(size=12),legend.position="bottom",legend.title = element_blank())

prow<-cowplot::plot_grid(
  aplot+ theme(legend.position="none"),
  bplot+ theme(legend.position="none"),
  labels = "AUTO",label_size = 20,label_fontfamily = "serif")

legend_b <- get_legend(
  aplot 
)

cowplot::plot_grid(prow, legend_b, ncol = 1, rel_heights = c(1, .2))

```

As we observed PCA do not show clear groups among the data. tSNE shows that cell lines are different between each other however, each cell line has a clear difference between parental and mesenchymal samples. Parental corresponds to cell lines prior to EMT induction, while Mesenchymal corresponds to cell lines after EMT induction.

# Differential analysis

We calculated the log-fold change and identified the metabolites that significantly changed their concentration before and after EMT induction, for each cell line.
As we need to do it for each cell line contained in the table we developed a function to calculate *p*-value and foldchange, returning a tribble with the data to graph the results.

```{r}
#only known metabolites have an assigned pathway
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
```


## A549

```{r A549,warning=F,message=F}
A549<-diffmet(EMT,"A549")

ggplot(A549,aes(x=foldchange, y=-log10(p.adj), colour=threshold)) +
  geom_point(size=3) +
  ggtitle(label = "A549")+
  xlab("log2 fold change") + 
  ylab("-log10 p-value") +
  theme(text = element_text(size=25),legend.position="none",
        plot.title = element_text(hjust = 0.5))+
  scale_color_manual(values=c("black", "red", "blue"))+
  geom_hline(yintercept = -log(0.05,10), linetype = "dashed",color="red") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed",color="red") 

#get list of metabolites differentially present
subset(A549, p.adj<.05 & abs(foldchange)>1.5)[,1]
```

## HCC827

```{r HCC827, warning=FALSE,message=F}
HCC827<-diffmet(EMT,"HCC")

ggplot(HCC827,aes(x=foldchange, y=-log10(p.adj), colour=threshold)) +
  geom_point(size=3) +
  ggtitle(label = "HCC827")+
  xlab("log2 fold change") + 
  ylab("-log10 p-value") +
  theme(text = element_text(size=25),legend.position="none",
        plot.title = element_text(hjust = 0.5))+
  scale_color_manual(values=c("black", "red", "blue"))+
  geom_hline(yintercept = -log(0.05,10), linetype = "dashed",color="red") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed",color="red") 

#get list of metabolites differentially present
subset(HCC827, p.adj<.05 & abs(foldchange)>1.5)[,1]
```

## NCI-H358

```{r H358, warning=FALSE,message=F}
H358<-diffmet(EMT,"NCI")

ggplot(H358,aes(x=foldchange, y=-log10(p.adj), colour=threshold)) +
  geom_point(size=3) +
  ggtitle(label = "NCI-H358")+
  xlab("log2 fold change") + 
  ylab("-log10 p-value") +
  theme(text = element_text(size=25),legend.position="none",
        plot.title = element_text(hjust = 0.5))+
  scale_color_manual(values=c("black", "red", "blue"))+
  geom_hline(yintercept = -log(0.05,10), linetype = "dashed",color="red") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed",color="red") 

#get list of metabolites differentially present
subset(H358, p.adj<.05 & abs(foldchange)>1.5)[,1]

```

Let´s calculate how many metabolites change in common between cell lines. We represented this on a Venn diagram.

```{r Venn,echo=F,fig.cap="Right side of the venn diagram shows the metabolites that change in common in all the cell lines"}
#eBayes data
A549$BIOCHEMICAL[which(!A549$threshold==0)]->a
HCC827$BIOCHEMICAL[which(!HCC827$threshold==0)]->b
H358$BIOCHEMICAL[which(!H358$threshold==0)]->c

Matrix<-matrix(0,nrow=dim(A549)[1],ncol=3)
#sort(unique(c(a,b,c)))->todos
rownames(Matrix)<-A549$BIOCHEMICAL
colnames(Matrix)<-c("A549","HCC827","NCI358")
Matrix[which(rownames(Matrix) %in% a),1]<-1
Matrix[which(rownames(Matrix) %in% b),2]<-1
Matrix[which(rownames(Matrix) %in% c),3]<-1

limma::vennCounts(Matrix)->conteo

as.data.frame(Matrix)->Todo
Todo[which(Todo$A549==1&Todo$HCC827==1&Todo$NCI358==1),]

library(ggforce)
library(grid)

df.venn <- data.frame(x = c(-1, 1, 0),
                      y = c(1, 1, -.7),
                      labels = colnames(Todo))

vdc <- conteo
class(vdc) <- 'matrix'
df.vdc <- as.data.frame(vdc) %>%
  mutate(x = c(2,0, 1.2, 0.7, -1.2, -0.7, 0, 0),
         y = c(-2,-1, 1.2, 0, 1.2, 0, 1.25, 0.4))

V<-ggplot(df.venn) +
  geom_circle(aes(x0 = x, y0 = y, r = 1.5, fill = labels,color=labels), alpha = .5,size=1) +
  coord_fixed() +
  theme_void() +
  theme(legend.position = 'none') +
  labs(fill = NULL) +
  annotate("text", x = df.vdc$x, y = df.vdc$y, label = df.vdc$Counts)+
  geom_label(label=c("A549","HCC827","NCI-H358"),x=c(-1,1,0),y=c(2.8,2.8,-2.5),label.size = NA,alpha=0)+
  ylim(-2.5,3)

metcom<- rownames(Todo[which(Todo$A549==1&Todo$HCC827==1&Todo$NCI358==1),])
metcom<-data.frame(met=metcom,
                   x=rep(0,length(metcom)),y=length(metcom):1)
metcom$met<-gsub("(^[[:alpha:]])", "\\U\\1", metcom$met, perl=TRUE)

commonmet<-ggplot(metcom,aes(x,y))+
  geom_text(aes(label=met))+theme_void()+ylim(0,9)

cowplot::plot_grid(
  V,commonmet,ncol = 2,
  labels = c("",""),label_fontfamily = "serif"
)
```



Venn Diagram for Enrichment pathway results from [https://www.metaboanalyst.ca/](Metaboanalyst 4.0):

```{r,fig.cap="On the right side common enriched pathways among all the cell lines"}

#setwd("../Dropbox/inmegen/EMT/MSEA2/")

A549<-read.csv("input/EnrichmentA549.csv")
HCC<-read.csv("input/EnrichmentHCC827.csv")
NCI<-read.csv("input/EnrichmentNCIH358.csv")



as.character(A549[which(A549$Holm.p<.01),1])->a
as.character(HCC[which(HCC$Holm.p<.01),1])->b
as.character(NCI[which(NCI$Holm.p<.01),1])->c

length(sort(unique(c(a,b,c))))

Matrix<-matrix(0,nrow=length(A549$X),ncol=3)
rownames(Matrix)<-A549$X
colnames(Matrix)<-c("A549","HCC827","NCI358")

Matrix[which(rownames(Matrix) %in% a),1]<-1
Matrix[which(rownames(Matrix) %in% b),2]<-1
Matrix[which(rownames(Matrix) %in% c),3]<-1

library("limma")
vennCounts(Matrix)->conteo

as.data.frame(Matrix)->Todo
Todo[which(Todo$A549==1&Todo$HCC827==1&Todo$NCI358==1),]

library(ggforce)
library(grid)
library(tidyverse)

df.venn <- data.frame(x = c(-1, 1, 0),
                      y = c(1, 1, -.7),
                      labels = colnames(Todo))
vdc <- conteo
class(vdc) <- 'matrix'
##positions on Venn Diagram
df.vdc <- as.data.frame(vdc) %>%
  mutate(x = c(2,0, 1.2, 0.7, -1.2, -0.7, 0, 0),
         y = c(-2,-1, 1.2, 0, 1.2, 0, 1.25, 0.4))

VE<-ggplot(df.venn) +
  geom_circle(aes(x0 = x, y0 = y, r = 1.5, fill = labels,color=labels), alpha = .5) +
  coord_fixed() +
  theme_void() +
  theme(legend.position = 'none') +
  labs(fill = NULL) +
  annotate("text", x = df.vdc$x, y = df.vdc$y, label = df.vdc$Counts)+
  geom_label(label=c("A549","HCC827","NCI-H358"),x=c(-1,1,0),y=c(2.8,2.8,-2.5),label.size = NA,alpha=0)+
  ylim(-2.5,3)

#ggsave("output/VennEnrich.png",width = 6,height=6,units = "cm")

path<-rownames(Todo[which(Todo$A549==1&Todo$HCC827==1&Todo$NCI358==1),])

pathcom<-data.frame(path=path,
                   x=length(path),y=length(path):1)
pathcom$path<-gsub("(^[[:alpha:]])", "\\U\\1", pathcom$path, perl=TRUE)

commonpath<-ggplot(pathcom,aes(x,y))+
  geom_text(aes(label=path))+theme_void()


cowplot::plot_grid(VE,commonpath,rel_widths = c(4,6))

ggsave("output/Enrich.png",width = 18,height = 8,units = "cm")

```




```{r warning=FALSE}
feb4<-vroom::vroom("input/Feb4.txt")
feb19<-vroom::vroom("input/Feb19.txt")

##transform
feb4<-feb4 %>% 
  pivot_longer(cols = -Time) %>% 
  mutate(name=as_factor(name),condicion=case_when(grepl("SIN", name) ~ "SIN",
                                                  grepl("CON", name) ~ "CON"),
         cellline=case_when(grepl("A549",name)~"A549",grepl("dKO",name)~"dKO")) %>% 
  arrange(desc(condicion))%>% mutate(condicion=as_factor(condicion))


feb19<-feb19 %>% pivot_longer(cols = -Time)%>% 
  mutate(name=as_factor(name),condicion=case_when(grepl("SIN", name) ~ "SIN",
                                                  grepl("CON", name) ~ "CON"),
         cellline=case_when(grepl("A549",name)~"A549",grepl("dKO",name)~"dKO")) %>% 
  arrange(desc(condicion))%>% mutate(condicion=as_factor(condicion))

###graph
#colors=c("#E0A33D","#FF4335",  "#51AFDF","#2B45F5")
colors=c("red2","red2",  "gold2","gold2")

wound4<-ggplot(feb4,aes(Time,value,color=name))+
  geom_line(aes(linetype = condicion)) +
  geom_point(aes(shape=cellline),size=2) +
  scale_color_manual(values=colors,guide=F)+
  theme_bw()+ggtitle("8 days")+ylab("Wound closure %")+theme(legend.position="bottom")+
  ylim(0,45)+xlab("Time (h)")


wound21<-ggplot(feb19,aes(Time,value,color=name))+
  geom_line(aes(linetype = condicion)) +
  geom_point(aes(shape=cellline),size=2) +
  scale_color_manual(values=colors,guide=F)+
  theme_bw()+ggtitle("21 days")+ylab("Wound closure %")+theme(legend.position="bottom")+
  ylim(0,45)+xlab("Time (h)")

prow<-cowplot::plot_grid(
  wound4+ theme(legend.position="none"),
  wound21+ theme(legend.position="none"),
  labels = "AUTO",label_size = 20,label_fontfamily = "serif")

legend_b <- cowplot::get_legend(
  wound21 
)

cowplot::plot_grid(prow, legend_b, ncol = 1, rel_heights = c(1, .1))


ggsave("output/HeridaGraphs.png",units = "cm",width = 15,height = 8)
```

