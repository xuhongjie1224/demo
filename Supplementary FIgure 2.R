rm(list = ls()) 

mergeexp <- merge(GSE14975exp,GSE31281exp,by="symbol")
mergeexp <- merge(mergeexp,GSE41177exp,by="symbol")
mergeexp <- merge(mergeexp,GSE79768exp,by="symbol")

save(mergeexp , file = "mergestep3output.Rdata")

write.table(mergeexp,file="mergeexp.txt",sep='\t',quote = F,row.names = F)

mergeexp<-read.table("mergeexp.txt",header=T,sep="\t")
write.table(mergeexp,file="mergeexp.xls",sep='\t',quote = F)


rm(list = ls())

load("mergestep3output.Rdata")
mergeexp<- mergeexp[!duplicated(mergeexp$symbol),]
rownames(mergeexp) = mergeexp[,1]

exp<-mergeexp[,c(2:49)]


BiocManager::install("sva")
library(sva)
library(stringr)
Group = ifelse(str_detect(pd$title,"NC"),"NC","AF")
Group = factor(Group,levels = c("NC","AF"))


batch <- c(rep("A",13),rep("B",19),rep("C",6),rep("D",10))
mod = model.matrix(~Group)
outTab <- data.frame(ComBat(exp, batch,mod, par.prior=TRUE))  

re_mergeexp<-outTab
save(re_mergeexp,file = "re_mergeexp.Rdata")


rm(list = ls())  
load(file = "re_mergeexp.Rdata")

exp<-re_mergeexp

library(stringr)

group_list=c(rep("GSE79768",times=13),rep("GSE41177",times=19),
             rep("GSE31821",times=6),rep("GSE14975",times=10))

group_list = factor(group_list,
                    levels = c("GSE79768","GSE41177","GSE31821","GSE14975"))

dat=as.data.frame(t(exp))
dat[1:4,1:4]
library(FactoMineR)
library(factoextra) 

dat.pca <- PCA(dat, graph = FALSE)
pca_plot <- fviz_pca_ind(dat.pca,
                         geom.ind = "point", # show points only (nbut not "text")
                         col.ind = group_list, # color by groups
                         #palette = c("#00AFBB", "#E7B800"),
                         addEllipses = TRUE, # Concentration ellipses
                         legend.title = "Groups"
)
pca_plot


rm(list = ls())  

load("mergestep3output.Rdata")
mergeexp<- mergeexp[!duplicated(mergeexp$symbol),]
rownames(mergeexp) = mergeexp[,1]

exp<-mergeexp[,c(2:49)]




library(stringr)

group_list=c(rep("GSE79768",times=13),rep("GSE41177",times=19),
             rep("GSE31821",times=6),rep("GSE14975",times=10))

group_list = factor(group_list,
                    levels = c("GSE79768","GSE41177","GSE31821","GSE14975"))

dat=as.data.frame(t(exp))
dat[1:4,1:4]
library(FactoMineR)
library(factoextra) 

dat.pca <- PCA(dat, graph = FALSE)
pca_plot <- fviz_pca_ind(dat.pca,
                         geom.ind = "point", # show points only (nbut not "text")
                         col.ind = group_list, # color by groups
                         #palette = c("#00AFBB", "#E7B800"),
                         addEllipses = TRUE, # Concentration ellipses
                         legend.title = "Groups"
)
pca_plot
ggsave(plot = pca_plot,filename = paste0(gse,"PCA.png"))
save(pca_plot,file = "pca_plot.Rdata")
