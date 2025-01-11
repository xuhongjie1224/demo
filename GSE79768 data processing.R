###GSE79768的原始数据处理
rm(list = ls()) 


BiocManager::install("ggplot2")
library(affyPLM)
library(affy)
library(RColorBrewer)
library(impute)
library(limma)
library(pheatmap)
library(ggplot2)
   
setwd("C:\\Users\\Xuhongjie\\Desktop\\AFMMD1\\AF\\GSE79768\\NC")
library(affyPLM)
library(affy)
Data<-ReadAffy()
sampleNames(Data)
N=length(Data)
eset.rma<-rma(Data)
Young_exprs<-exprs(eset.rma)
probeid<-rownames(Young_exprs)
Young_exprs<-cbind(probeid,Young_exprs)
write.table(Young_exprs,file="NC_exprs.txt",sep='\t',quote = F,row.names = F)

setwd("C:\\Users\\Xuhongjie\\Desktop\\AFMMD1\\AF\\GSE79768\\AF")
library(affyPLM)
library(affy)
Data<-ReadAffy()
sampleNames(Data)
N=length(Data)
eset.rma<-rma(Data)
Old_exprs<-exprs(eset.rma)
probeid<-rownames(Old_exprs)
Old_exprs<-cbind(probeid,Old_exprs)
write.table(Old_exprs,file="AF_exprs.txt",sep='\t',quote = F,row.names = F)


setwd("C:\\Users\\Xuhongjie\\Desktop\\AFMMD1\\AF\\GSE79768\\cel")
NC_exprs<-read.table("NC_exprs.txt",header=T,sep="\t")
AF_exprs<-read.table("AF_exprs.txt",header=T,sep="\t")
GSE79768.probe_exprs<-merge(NC_exprs,AF_exprs,by="probeid")
write.table(GSE79768.probe_exprs,file="GSE79768.probeid_exprs.txt",sep='\t',quote = F,row.names = F)

RAND.probeid_exprs<-read.table("GSE79768.probeid_exprs.txt",header=T,sep="\t")
write.table(RAND.probeid_exprs,file="GSE79768.probeid_exprs.xls",sep='\t',quote = F)



exp<- GSE79768.probe_exprs
X=as.matrix(exp[2:14])
rownames(X)<-exp$probeid
exp<-X
exp[1:4,1:4] 
range(exp)


boxplot(exp)

save(GSE79768.probe_exprs,file = "GSE79768step1output.Rdata")

rm(list = ls())  
load(file = "GSE79768step1output.Rdata")

probe_exprs<-GSE79768.probe_exprs



library(stringr)

GPL97pd <- read_excel("cel/GPL97pd.xlsx")
pd<-GPL97pd_ 
pd$title

library(stringr)
group_list=ifelse(str_detect(pd$title,"NC"),"Normal","AF")

group_list = factor(group_list,
                    levels = c("Normal","AF"))

gpl 

if(!require(hgu133plus2.db))BiocManager::install("hgu133plus2.db")
library(hgu133plus2.db)
ls("package:hgu133plus2.db")
ids <- toTable(hgu133plus2SYMBOL)
head(ids)
ids=as.matrix(ids[1:2])

colnames(probe_exprs)[1] <-"probe_id" 
exp <- merge(ids,probe_exprs,by="probe_id")
head(exp)

GSE79768exp<-exp

GSE79768exp<-GSE79768exp[,c(2:15)]
GSE79768ids<-ids
GSE79768group_list<-group_list
save(GSE79768exp,GSE79768group_list,GSE79768ids,
     file = "GSE79768step2output.Rdata")

rm(list = ls())  
setwd("C:\\Users\\Xuhongjie\\Desktop\\RAMDD")

write.table(mergeexp,file="mergeRA.txt",sep='\t',quote = F,row.names = F)

write.table(mergeexp,file="mergeRA.xls",sep='\t',quote = F)

mergeexp<-GSE79768exp
mergeexp<- mergeexp[!duplicated(mergeexp$symbol),]

rownames(mergeexp) = mergeexp[,1]

exp<-mergeexp[,c(2:14)]

pd <- read_excel("pd.xlsx")
group_list<-GSE79768group_list

library(stringr)
group_list=ifelse(str_detect(pd$title,"ND"),"Normal","RA")

group_list = factor(group_list,
                    levels = c("Normal","RA"))


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

save(exp,group_list,deg,logFC_t,P.Value_t,file = "GSE79768step4output.Rdata")

