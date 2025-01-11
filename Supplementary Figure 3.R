###ssGSEA analysis of SLE dataset
rm(list = ls()) 

library(dplyr)
load("mergestep4output.Rdata")
group_list
NC<-exp[,c(1,2,3,4,5,6,14,15,16,33,34,39,40,41,42,43)]
NC<-mutate(NC,symbol=rownames(NC))
AF<-exp[,-c(1,2,3,4,5,6,14,15,16,33,34,39,40,41,42,43)]
AF<-mutate(AF,symbol=rownames(AF))
exp_NC_AF <- merge(NC,AF,by="symbol")
rownames(exp_NC_AF) = exp_NC_AF[,1]
exp_NC_AF<-exp_NC_AF[,c(2:49)]
exp<-exp_NC_AF

BiocManager::install("GSVA")
library(GSVA)



geneset <- mmc3

geneset2 <- split(geneset$Metagene,geneset$`Cell type`)
View(geneset2)
geneset2[1:2]

res <- gsva(
  as.matrix(exp), 
  geneset2, 
  method = "ssgsea", 
  kcdf = "Gaussian",
  mx.diff = F,
  verbose   = F 
)

resm <- res
for (i in colnames(res)) {
  resm[,i] <- (res[,i] -min(res[,i]))/(max(res[,i] )-min(res[,i] ))
}


View(res)
View(resm)
save(res,resm,file = c("KIRC_immune_cell_result.Rdata")) 

library(pheatmap)

pheatmap(res, show_colnames = F)


library(stringr)
group_list=c(rep("Normal",times=16),rep("AF",times=32))
group_list = factor(group_list,
                    levels = c("Normal","AF"))
group_list
annotation <- data.frame(group_list)
rownames(annotation) <- colnames(res)
table(annotation)
head(annotation)

pheatmap(res,
         show_colnames = F,
         show_rownames = T,
         scale = "row",
         annotation_col = annotation,
         cluster_cols = F,
         fontsize = 10)

library('tidyr')
library('ggplot2')
library('ggpubr')
ssGSEA=res
mydata <- t(as.data.frame(ssGSEA))%>%as.data.frame()

mydata$group <- group_list
table(mydata$group) 
ggdata <- gather(mydata,key = 'Cell',value = 'Value',-29) 

ggplot(ggdata,aes(x=Cell,y=Value,fill=group))+
  geom_boxplot(width=0.7,size=0.3,outlier.color = NA)+
  theme_bw()+
  theme(panel.grid = element_blank())+
  stat_compare_means(symnum.args = list(cutpoints = c(0,0.001, 0.01, 0.05, 1), 
                                        symbols = c("***", "**", "*", "ns")),
                     label = "p.signif")+
  theme(axis.text.x = element_text(angle = 60,hjust = 1))+
  theme(legend.position = 'top')+
  xlab('')+ylab('Infiltration Abundance')+
  labs(fill='Group')


genes <- c( "ITGB2", "NFKBIA", "TMEM45A")


exp_genes <- exp[genes,]

rb <- rbind(res,exp_genes)
rb=t(rb)
rownames(rb)

library(ggplot2)
library(tinyarray)
library(GSVA)
library(dplyr)
library(Hmisc)
library(pheatmap)
library(ggpubr)
rb=as.matrix(rb)
m = rcorr(rb)$r[1:nrow(res),(ncol(rb)-length(genes)+1):ncol(rb)]

p = rcorr(rb)$P[1:nrow(res),(ncol(rb)-length(genes)+1):ncol(rb)]


library(corrplot)
split=m
splitp=p
splitp<-na.omit(splitp)

mark <- matrix(case_when(as.vector(splitp) < 0.001~"***",
                         as.vector(splitp)< 0.01~"**",
                         as.vector(splitp)< 0.05~"*",
                         T~""),
               nrow = nrow(splitp))
View(mark)


col2<-colorRampPalette(c("blue","white", "red"))(95)

library(pheatmap)
pheatmap(t(split),
         display_numbers = t(mark),
         number_color = "black",
         fontsize_number = 13,
         color = col2,
         border_color = "white",
         treeheight_col = 0,
         treeheight_row = 0,
         fontsize_col = 8,
         fontsize_row = 8,
         angle_col = "90")


###lollipop chart 
genes <- c("ITGB2", "NFKBIA","TMEM45A")

exp_genes <- exp[genes,]

rb <- rbind(res,exp_genes)
rb=t(rb)
rownames(rb)

library(ggplot2)
library(GSVA)
library(dplyr)
library(Hmisc)
library(pheatmap)
library(ggpubr)

m = rcorr(rb)$r[1:nrow(res),(ncol(rb)-length(genes)+1):ncol(rb)]

p = rcorr(rb)$P[1:nrow(res),(ncol(rb)-length(genes)+1):ncol(rb)]

library(corrplot)
split=m
splitp=p
splitp<-na.omit(splitp)
write.table(splitp,file="splitp.xls",sep='\t',quote = F)
write.table(split,file="split.xls",sep='\t',quote = F)


library(ggplot2)
library(forcats)
A<-TMEM45A
A=A[order(A$cor, decreasing=F),]
A$cell<- as.factor(A$cell)
A$cell <- fct_inorder(A$cell)

p <- ggplot(A,aes(x=cell,y=cor))+
  geom_point(data=A,aes(size=abs(cor),color=p.value))+
  scale_size(rang = c(0,6))+
  scale_color_viridis_c()+  
  geom_segment(aes(x=cell,xend=cell,y=0,yend=cor),
               size=1,linetype="solid")+
  labs(x="", y='Correlation(r)',title = 'TMEM45A')+
  coord_flip()+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5),
        panel.border = element_blank(),
        axis.text =element_text(size = 10, color = "black"),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust=0.5))
p

p1 <- p +geom_point(data=A,aes(size=cor,color=p.value))+guides(size=guide_legend(title="Cor"))
p1

p2<-p1 + annotate(geom = "text", x = unique(A$cell),
                  label = A$p.value,
                  y =1.02, hjust = 0)
p2


library(ggplot2)
library(forcats)
A<-ITGB2
A=A[order(A$cor, decreasing=F),]
A$cell<- as.factor(A$cell)
A$cell <- fct_inorder(A$cell)

p <- ggplot(A,aes(x=cell,y=cor))+
  geom_point(data=A,aes(size=abs(cor),color=p.value))+
  scale_size(rang = c(0,6))+
  scale_color_viridis_c()+  
  geom_segment(aes(x=cell,xend=cell,y=0,yend=cor),
               size=1,linetype="solid")+
  labs(x="", y='Correlation(r)',title = 'ITGB2')+
  coord_flip()+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5),
        panel.border = element_blank(),
        axis.text =element_text(size = 10, color = "black"),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust=0.5))
p

p1 <- p +geom_point(data=A,aes(size=cor,color=p.value))+guides(size=guide_legend(title="Cor"))
p1

p2<-p1 + annotate(geom = "text", x = unique(A$cell),
                  label = A$p.value,
                  y =1.02, hjust = 0)
p2

library(ggplot2)
library(forcats)
A<-NFKBIA
A=A[order(A$cor, decreasing=F),]
A$cell<- as.factor(A$cell)
A$cell <- fct_inorder(A$cell)

p <- ggplot(A,aes(x=cell,y=cor))+
  geom_point(data=A,aes(size=abs(cor),color=p.value))+
  scale_size(rang = c(0,6))+
  scale_color_viridis_c()+  
  geom_segment(aes(x=cell,xend=cell,y=0,yend=cor),
               size=1,linetype="solid")+
  labs(x="", y='Correlation(r)',title = 'NFKBIA')+
  coord_flip()+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5),
        panel.border = element_blank(),
        axis.text =element_text(size = 10, color = "black"),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust=0.5))
p

p1 <- p +geom_point(data=A,aes(size=cor,color=p.value))+guides(size=guide_legend(title="Cor"))
p1

p2<-p1 + annotate(geom = "text", x = unique(A$cell),
                  label = A$p.value,
                  y =1.02, hjust = 0)
p2
