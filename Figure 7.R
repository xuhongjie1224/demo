###GSEA analysis of TMEM45A. 
rm(list = ls()) 
load("mergestep4output.Rdata")
exprSet <- exp
test <- exprSet[1:10,1:10]

batch_cor <- function(gene){
  y <- as.numeric(exprSet[gene,])
  rownames <- rownames(exprSet)
  do.call(rbind,future_lapply(rownames, function(x){
    dd  <- cor.test(as.numeric(exprSet[x,]),y,type="spearman")
    data.frame(gene=gene,mRNAs=x,cor=dd$estimate,p.value=dd$p.value )
  }))
}

library(future.apply)
plan(multiprocess)
system.time(dd <- batch_cor("TMEM45A"))

gene <- dd$mRNAs

library(clusterProfiler)
gene = bitr(gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)

gene_df <- data.frame(logFC=dd$cor,
                      SYMBOL = dd$mRNAs)
gene_df <- merge(gene_df,gene,by="SYMBOL")


geneList <- gene_df$logFC

names(geneList) = gene_df$ENTREZID

geneList = sort(geneList, decreasing = TRUE)

library(clusterProfiler)

hallmarks <- read.gmt("c2.cp.kegg.v7.4.entrez.gmt")

y <- GSEA(geneList,TERM2GENE =hallmarks)

library(ggplot2)
dotplot(y,showCategory=30,label_format=100,split=".sign")+facet_grid(~.sign)


###GSEA analysis of single pathway
yd <- data.frame(y)
library(enrichplot)
gseaplot2(y,"KEGG N GLYCAN BIOSYNTHESIS",color = "red",pvalue_table = T)
gseaplot2(y,"KEGG RIBOSOME",color = "red",pvalue_table = T)
gseaplot2(y,"KEGG SPHINGOLIPID METABOLISM",color = "red",pvalue_table = T)
gseaplot2(y,"KEGG FC GAMMA R MEDIATED PHAGOCYTOSIS",color = "red",pvalue_table = T)
gseaplot2(y,"KEGG_TGF_BETA_SIGNALING_PATHWAY",color = "red",pvalue_table = T)
