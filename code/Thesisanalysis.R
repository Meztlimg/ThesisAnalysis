
#R version 3.6.3
#

library(Rtsne) # implementation of t-SNE algorithm
library(RColorBrewer) # library to access easily multiple colors for plotting

library(tidyverse)
library(rstatix) #Pipe-Friendly Framework for Basic Statistical Tests
library(ggfortify)


#import data
EMT <- vroom::vroom("input/Metabolomica_antes-despues_EMT.txt")

EMT <- janitor::clean_names(dat = EMT, 
                            case = "screaming_snake")
EMT<- EMT %>% filter(!is.na(BIOCHEMICAL))
##############################
#PCA
PCAMatrix<-EMT %>% 
  select(contains(c("A549","HCC","NCI"))) %>% drop_na()

pca = prcomp(t(PCAMatrix))
nametable<-as.data.frame(cbind(colnames(PCAMatrix),gsub("\\_[1-9]","",colnames(PCAMatrix))))
nametable$V2<-as.factor(nametable$V2)
nametable$V2<-factor(nametable$V2,levels = levels(nametable$V2)[c(2,1,4,3,6,5)])

aplot<-autoplot(pca,data=nametable,colour='V2',frame=F,size=3)+
  scale_color_brewer(palette = "Paired")+
  theme(text = element_text(size=12),legend.position="bottom",legend.title = element_blank())

#ggsave("ThesisPCA.png")


##############################
#t-SNE
print('running t-SNE...')
# 2.2.1. low resolution of t-SNE
print('running 2D t-SNE with small perplexity...')
set.seed(80)
results=Rtsne(t(PCAMatrix),dims=2,perplexity=5,verbose=TRUE,max_iter=10000)
# plotting
bplot<-ggplot(as.data.frame(results$Y),aes(x=results$Y[,1],y=results$Y[,2],col=nametable$V2))+
  xlab("tSNE1") + ylab("tSNE2")+
  geom_point(size=3)+ #Size and alpha just for fun
  scale_color_brewer(palette = "Paired")+ #your colors here
  theme(text = element_text(size=12),legend.position="bottom",legend.title = element_blank())

#ggsave("ThesisTSNE.png")

prow<-cowplot::plot_grid(
  aplot+ theme(legend.position="none"),
  bplot+ theme(legend.position="none"),
  labels = "AUTO",label_size = 20,label_fontfamily = "serif")

legend_b <- cowplot::get_legend(
  aplot 
)

cowplot::plot_grid(prow, legend_b, ncol = 1, rel_heights = c(1, .2))

ggsave("output/Clustering.png",height = 10,width = 19,units = "cm")

###############################
#Differential analysis

source("code/ThesisFunctions.R")

#A549 ######################

A549<-diffmet(EMT,"A549")

A<-ggplot(A549,aes(x=foldchange, y=-log10(p.adj), colour=threshold)) +
  geom_point() +
  ggtitle(label = "A549")+
  xlab("log2 fold change") + 
  ylab("-log10 p-value") +
  theme(legend.position="none",
        plot.title = element_text(hjust = 0.5))+
  scale_color_manual(values=c("black", "red", "blue"))+
  geom_hline(yintercept = -log(0.05,10), linetype = "dashed",color="red") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed",color="red") 

#get list of metabolites differentially present
subset(A549, p.adj<.05 & abs(foldchange)>1.5)[,1]


#HCC827 ######################

HCC827<-diffmet(EMT,"HCC")

H<-ggplot(HCC827,aes(x=foldchange, y=-log10(p.adj), colour=threshold)) +
  geom_point() +
  ggtitle(label = "HCC827")+
  xlab("log2 fold change") + 
  ylab("-log10 p-value") +
  theme(legend.position="none",
        plot.title = element_text(hjust = 0.5))+
  scale_color_manual(values=c("black", "red", "blue"))+
  geom_hline(yintercept = -log(0.05,10), linetype = "dashed",color="red") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed",color="red") 

#get list of metabolites differentially present
subset(HCC827, p.adj<.05 & abs(foldchange)>1.5)[,1]

#NCI-H358 ####################
H358<-diffmet(EMT,"NCI")

N<-ggplot(H358,aes(x=foldchange, y=-log10(p.adj), colour=threshold)) +
  geom_point() +
  ggtitle(label = "NCI-H358")+
  xlab("log2 fold change") + 
  ylab("-log10 p-value") +
  theme(legend.position="none",
        plot.title = element_text(hjust = 0.5))+
  scale_color_manual(values=c("black", "red", "blue"))+
  geom_hline(yintercept = -log(0.05,10), linetype = "dashed",color="red") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed",color="red") 

#get list of metabolites differentially present
subset(H358, p.adj<.05 & abs(foldchange)>1.5)[,1]




# Venn Diagram ######################

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
#png("MetabComun2020conocidos.png",res = 300,units = "cm",width = 15,height = 15)
limma::vennDiagram(conteo,circle.col = c("red","blue","green"))
#dev.off()

library(ggforce)
library(grid)

as.data.frame(Matrix)->Todo
Todo[which(Todo$A549==1&Todo$HCC827==1&Todo$NCI358==1),]

df.venn <- data.frame(x = c(-1, 1, 0),
                      y = c(1, 1, -.7),
                      labels = colnames(Todo))
ggplot(df.venn, aes(x0 = x, y0 = y, r = 1.5, fill = labels)) +
  geom_circle(alpha = .3) +
  coord_fixed() +
  theme_void()

vdc <- conteo
class(vdc) <- 'matrix'
##positions on Venn Diagram
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
  A,H,N,NULL,V,commonmet,ncol = 3,
  labels = c("A","","","","B"),label_fontfamily = "serif",label_size = 16,
  axis = "b"
)

ggsave("output/Metabolites.png",width = 20,height = 11,units = "cm")


### Enrichment analysis