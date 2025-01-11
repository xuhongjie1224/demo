rm(list = ls()) 

BiocManager::install("ggplot2")
library(affyPLM)
library(affy)
library(RColorBrewer)
library(impute)
library(limma)
library(pheatmap)
library(ggplot2)


setwd("C:\\Users\\Xuhongjie\\Desktop\\SLE\\GSE50772\\NC")
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


setwd("C:\\Users\\Xuhongjie\\Desktop\\SLE\\GSE50772\\SLE")
library(affyPLM)
library(affy)
Data<-ReadAffy()
sampleNames(Data)
N=length(Data)
eset.rma<-rma(Data)
Old_exprs<-exprs(eset.rma)
probeid<-rownames(Old_exprs)
Old_exprs<-cbind(probeid,Old_exprs)
write.table(Old_exprs,file="SLE_exprs.txt",sep='\t',quote = F,row.names = F)

setwd("C:\\Users\\Xuhongjie\\Desktop\\SLE\\GSE50772\\cel")
NC_exprs<-read.table("NC_exprs.txt",header=T,sep="\t")
SLE_exprs<-read.table("SLE_exprs.txt",header=T,sep="\t")
GSE50772.probe_exprs<-merge(NC_exprs,SLE_exprs,by="probeid")
write.table(GSE50772.probe_exprs,file="GSE50772.probeid_exprs.txt",sep='\t',quote = F,row.names = F)

RAND.probeid_exprs<-read.table("GSE50772.probeid_exprs.txt",header=T,sep="\t")
write.table(RAND.probeid_exprs,file="GSE50772.probeid_exprs.xls",sep='\t',quote = F)

exp<- GSE50772.probe_exprs
X=as.matrix(exp[2:82])
rownames(X)<-exp$probeid
exp<-X
exp[1:4,1:4] 
range(exp)


boxplot(exp)

save(GSE50772.probe_exprs,file = "GSE50772step1output.Rdata")

rm(list = ls())  
load(file = "GSE50772step1output.Rdata")

probe_exprs<-GSE50772.probe_exprs



library(stringr)

GPL97pd <- read_excel("cel/GPL97pd.xlsx")
pd<-GPL97pd_ 
pd$title


library(stringr)
group_list=ifelse(str_detect(pd$title,"NC"),"Normal","SLE")

group_list = factor(group_list,
                    levels = c("Normal","SLE"))

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

GSE50772exp<-exp

GSE50772exp<-GSE50772exp[,c(2:83)]
GSE50772ids<-ids
GSE50772group_list<-group_list
save(GSE50772exp,GSE50772group_list,GSE50772ids,exp,
     file = "GSE50772step2output.Rdata")



