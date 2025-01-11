###AF-DEG
rm(list = ls())  
setwd("C:\\Users\\Xuhongjie\\Desktop\\RAMDD")
load("re_mergeexp.Rdata")

exp<-re_mergeexp
library(stringr)
group_list=ifelse(str_detect(pd$title,"NC"),"Normal","AF")

group_list = factor(group_list,
                    levels = c("Normal","AF"))

library(limma)
design=model.matrix(~group_list)
fit=lmFit(exp,design)
fit=eBayes(fit)
deg=topTable(fit,coef=2,number = Inf)

library(dplyr)
deg <- mutate(deg,symbol=rownames(deg))
head(deg)

logFC_t=1
P.Value_t = 0.05
k1 = (deg$P.Value< P.Value_t)&(deg$logFC < -logFC_t)
k2 = (deg$P.Value< P.Value_t)&(deg$logFC > logFC_t)
change = ifelse(k1,"down",ifelse(k2,"up","stable"))
deg <- mutate(deg,change)

library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
s2e <- bitr(deg$symbol, 
            fromType = "SYMBOL",
            toType = "ENTREZID",
            OrgDb = org.Hs.eg.db)

deg <- inner_join(deg,s2e,by=c("symbol"="SYMBOL"))

save(exp,group_list,deg,logFC_t,P.Value_t,file = "mergestep4output.Rdata")

cg = deg$symbol[deg$change !="stable"]
length(cg)
n=exp[cg,]
degexp<-n

save(exp,group_list,deg,logFC_t,P.Value_t,degexp,file = "mergestep4output.Rdata")

#AF-pheatmap plot and  volcano plot
rm(list = ls()) 
load(file = "mergestep4output.Rdata")

library(dplyr)
library(ggplot2)
dat  = deg

p <- ggplot(data = dat, 
            aes(x = logFC, 
                y = -log10(P.Value))) +
  geom_point(alpha=0.4, size=3.5, 
             aes(color=change)) +
  ylab("-log10(Pvalue)")+
  scale_color_manual(values=c("grey", "red","blue"))+
  geom_vline(xintercept=c(-logFC_t,logFC_t),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(P.Value_t),lty=4,col="black",lwd=0.8) +
  theme_bw()
p


if(F){

  for_label <- dat%>% 
    filter(symbol %in% c("GPX4","TFRC")) 
}
if(T){

  for_label <- dat %>% head(10)
}
if(F) {

  x1 = dat %>% 
    filter(change == "up") %>% 
    head(3)
  x2 = dat %>% 
    filter(change == "down") %>% 
    head(3)
  for_label = rbind(x1,x2)
}

volcano_plot <- p +
  geom_point(size = 3, shape = 1, data = for_label) +
  ggrepel::geom_label_repel(
    aes(label = symbol),
    data = for_label,
    color="black"
  )
volcano_plot
ggsave(plot = volcano_plot,filename = paste0(gse,"volcano.png"))


if(T){
 
  cg = deg$symbol[deg$change !="stable"]
  length(cg)
}else{
  
  x=deg$logFC[deg$change !="stable"] 
  names(x)=deg$probe_id[deg$change !="stable"] 
  cg=c(names(head(sort(x),30)),names(tail(sort(x),30)))
  length(cg)
}
n=exp[cg,]
dim(n)

group_list

library(pheatmap)
annotation_col=data.frame(group=group_list)
rownames(annotation_col)=colnames(n) 
library(ggplotify)
heatmap_plot <- as.ggplot(pheatmap(n,show_colnames =F,
                                   show_rownames = F,
                                   scale = "row",
                                   cluster_cols = F, 
                                   annotation_col=annotation_col)) 
heatmap_plot

###SLE-DEG
rm(list = ls())  
setwd("C:\\Users\\Xuhongjie\\Desktop\\SLE\\GSE50772")

mergeexp<-GSE50772exp
mergeexp<- mergeexp[!duplicated(mergeexp$symbol),]

rownames(mergeexp) = mergeexp[,1]

exp<-mergeexp[,c(2:82)]


library(stringr)
group_list=ifelse(str_detect(pd$title,"NC"),"Normal","SLE")

group_list = factor(group_list,
                    levels = c("Normal","SLE"))

library(limma)
design=model.matrix(~group_list)
fit=lmFit(exp,design)
fit=eBayes(fit)
deg=topTable(fit,coef=2,number = Inf)


library(dplyr)
deg <- mutate(deg,symbol=rownames(deg))
head(deg)

logFC_t=1
P.Value_t = 0.05
k1 = (deg$P.Value< P.Value_t)&(deg$logFC < -logFC_t)
k2 = (deg$P.Value< P.Value_t)&(deg$logFC > logFC_t)
change = ifelse(k1,"down",ifelse(k2,"up","stable"))
deg <- mutate(deg,change)

library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
s2e <- bitr(deg$symbol, 
            fromType = "SYMBOL",
            toType = "ENTREZID",
            OrgDb = org.Hs.eg.db)

deg <- inner_join(deg,s2e,by=c("symbol"="SYMBOL"))

cg = deg$symbol[deg$change !="stable"]
length(cg)
n=exp[cg,]
degexp<-n

save(degexp,exp,group_list,deg,logFC_t,P.Value_t,file = "GSE50772step4output.Rdata")

###SLE-pheatmap plot and  volcano plot
rm(list = ls()) 
load(file = "GSE50772step4output.Rdata")
#1.火山图----
library(dplyr)
library(ggplot2)
dat  = deg

p <- ggplot(data = dat, 
            aes(x = logFC, 
                y = -log10(P.Value))) +
  geom_point(alpha=0.4, size=3.5, 
             aes(color=change)) +
  ylab("-log10(Pvalue)")+
  scale_color_manual(values=c("blue", "grey","red"))+
  geom_vline(xintercept=c(-logFC_t,logFC_t),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(P.Value_t),lty=4,col="black",lwd=0.8) +
  theme_bw()
p


if(F){
  
  for_label <- dat%>% 
    filter(symbol %in% c("GPX4","TFRC")) 
}
if(T){
  
  for_label <- dat %>% head(10)
}
if(F) {
  
  x1 = dat %>% 
    filter(change == "up") %>% 
    head(3)
  x2 = dat %>% 
    filter(change == "down") %>% 
    head(3)
  for_label = rbind(x1,x2)
}

volcano_plot <- p +
  geom_point(size = 3, shape = 1, data = for_label) +
  ggrepel::geom_label_repel(
    aes(label = symbol),
    data = for_label,
    color="black"
  )
volcano_plot
ggsave(plot = volcano_plot,filename = paste0(gse,"volcano.png"))



load(file = 'GSE50772step4output.Rdata')
if(T){
  cg = deg$symbol[deg$change !="stable"]
  length(cg)
}else{
  x=deg$logFC[deg$change !="stable"] 
  names(x)=deg$symbol[deg$change !="stable"] 
  cg=c(names(head(sort(x),30)),names(tail(sort(x),30)))
  length(cg)
}

n=exp[cg,]
dim(n)

library(pheatmap)
annotation_col=data.frame(group=group_list)
rownames(annotation_col)=colnames(n) 
library(ggplotify)
heatmap_plot <- as.ggplot(pheatmap(n,show_colnames =F,
                                   show_rownames = F,
                                   scale = "row",
                                  cluster_cols = F, 
                                   annotation_col=annotation_col)) 
heatmap_plot
ggsave(heatmap_plot,filename = paste0(gse,"heatmap.png"))
load("pca_plot.Rdata")
library(patchwork)
(pca_plot + volcano_plot +heatmap_plot)+ plot_annotation(tag_levels = "A")

