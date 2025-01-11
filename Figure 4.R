### 69 core genes
rm(list = ls()) 
load("GSE50772step4output.Rdata")
SLE_deg<-rownames(degexp)

load("mergestep4output.Rdata")
AF_deg<-rownames(degexp)

k1 = AF_deg%in% SLE_deg;table(k1)
DEG<- AF_deg[k1]

write.table(AF_deg,
            file="AF_deg.txt",
            row.names = F,
            col.names = F,
            quote = F)

write.table(SLE_deg,
            file="SLE_deg.txt",
            row.names = F,
            col.names = F,
            quote = F)

DEG_PPI<-c(DEG,PPI$X1)
DEG_PPI<-DEG_PPI[!duplicated(DEG_PPI)]
write.table(DEG_PPI,
            file="DEG_PPI.txt",
            row.names = F,
            col.names = F,
            quote = F)

save(DEG,DEG_PPI,
     file = "DEG_PPI.Rdata")

###GO of 69 core genes
rm(list = ls()) 
load("mergestep4output.Rdata")
load("DEG_PPI.Rdata")

library(dplyr)
AF_exp <- mutate(exp,symbol=rownames(exp))
k1 = AF_exp$symbol%in% DEG_PPI;table(k1)
AF_lasso= AF_exp[k1,]

library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
s2e <- bitr(AF_lasso$symbol, 
            fromType = "SYMBOL",
            toType = "ENTREZID",
            OrgDb = org.Hs.eg.db)
AF_fj <- inner_join(AF_lasso,s2e,by=c("symbol"="SYMBOL"))

gene_all = AF_fj[,'ENTREZID']


if(T){
  ego_CC <- enrichGO(gene = gene_all,
                     OrgDb= org.Hs.eg.db,
                     ont = "CC",
                     pAdjustMethod = "BH",
                     minGSSize = 1,
                     pvalueCutoff = 0.9,
                     qvalueCutoff = 0.9,
                     readable = TRUE)
  ego_BP <- enrichGO(gene = gene_all,
                     OrgDb= org.Hs.eg.db,
                     ont = "BP",
                     pAdjustMethod = "BH",
                     minGSSize = 1,
                     pvalueCutoff = 0.9,
                     qvalueCutoff = 0.9,
                     readable = TRUE)
  ego_MF <- enrichGO(gene = gene_all,
                     OrgDb= org.Hs.eg.db,
                     ont = "MF",
                     pAdjustMethod = "BH",
                     minGSSize = 1,
                     pvalueCutoff = 0.9,
                     qvalueCutoff = 0.9,
                     readable = TRUE)
  save(ego_CC,ego_BP,ego_MF,file = "ego_AF_fj.Rdata")
}
load(file = "ego_AF_fj.Rdata")

dotplot(ego_BP,showCategory=20,label_format=100)


dotplot(ego_MF,showCategory=20,label_format=100)

dotplot(ego_CC,showCategory=20,label_format=100)


##KEGG pathway analysis of 69 core genes

library(R.utils)

R.utils::setOption("clusterProfiler.download.method",'auto')
if(T){
  
  AF_fj<- enrichKEGG(gene         = gene_all,
                        organism     = 'hsa',
                        pvalueCutoff = 0.05)
  save(AF_fj,file = "AF_fj_kegg.Rdata")
}
load("AF_fj_kegg.Rdata")
dotplot(AF_fj,showCategory=20,label_format=100)
