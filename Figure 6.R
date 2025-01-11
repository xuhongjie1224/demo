###Lasso analysis of SLE
rm(list = ls()) 
load("DEG_PPI.Rdata")
load("GSE50772step4output.Rdata")

library(dplyr)
SLE_exp <- mutate(exp,symbol=rownames(exp))
k1 = SLE_exp$symbol%in% DEG_PPI;table(k1)
SLE_lasso= SLE_exp[k1,]
SLE_lasso<-SLE_lasso[,-c(82)]

meta <- read.csv("ClinicalTraits.csv")

x=t(SLE_lasso)
y=meta$SLE
library(glmnet)
model_lasso <- glmnet(x, y,nlambda=1000, alpha=1)
print(model_lasso)
View(print(model_lasso))
NROW(model_lasso$lambda)
model_lasso$lambda[782]
plot(model_lasso)
plot(model_lasso,xvar = "lambda")

set.seed(1234)

cv_fit <- cv.glmnet(x=x, y=y, nlambda = 1000,alpha = 1)
plot(cv_fit)

model_lasso_min <- glmnet(x=x, y=y, alpha = 1, lambda=cv_fit$lambda.min)
head(model_lasso_min$beta,20)
choose_gene_min=rownames(model_lasso_min$beta)[as.numeric(model_lasso_min$beta)!=0]
length(choose_gene_min)
SLE_choose_gene_min<-choose_gene_min

save(SLE_choose_gene_min,
     file = "SLE_lasso_result.Rdata")


###Lasso analysis of AF
load("mergestep4output.Rdata")

library(dplyr)
AF_exp <- mutate(exp,symbol=rownames(exp))
k1 = AF_exp$symbol%in% DEG_PPI;table(k1)
AF_lasso= AF_exp[k1,]
AF_lasso<-AF_lasso[,-c(49)]

meta <- read.csv("AF_ClinicalTraits.csv")

x=t(AF_lasso)
y=meta$AF
library(glmnet)
model_lasso <- glmnet(x, y,nlambda=1000, alpha=1)
print(model_lasso)
View(print(model_lasso))
NROW(model_lasso$lambda)
model_lasso$lambda[1000]
plot(model_lasso)
plot(model_lasso,xvar = "lambda")

set.seed(1234)
cv_fit <- cv.glmnet(x=x, y=y, nlambda = 1000,alpha = 1)
plot(cv_fit)

model_lasso_min <- glmnet(x=x, y=y, alpha = 1, lambda=cv_fit$lambda.min)
head(model_lasso_min$beta,20)
choose_gene_min=rownames(model_lasso_min$beta)[as.numeric(model_lasso_min$beta)!=0]
length(choose_gene_min)
AF_choose_gene_min<-choose_gene_min

save(AF_choose_gene_min,
     file = "AF_lasso_result.Rdata")

k2 = AF_choose_gene_min %in% SLE_choose_gene_min;table(k2)
SLE_AF_lasso= AF_choose_gene_min[k2]
SLE_AF_lasso

write.table(AF_choose_gene_min,
            file="AF_choose_gene_min.txt",
            row.names = F,
            col.names = F,
            quote = F)

write.table(SLE_choose_gene_min,
            file="SLE_choose_gene_min.txt",
            row.names = F,
            col.names = F,
            quote = F)

write.table(SLE_AF_lasso,
            file="SLE_AF_lasso.txt",
            row.names = F,
            col.names = F,
            quote = F)

###SYM-RFE analysis of SLE
rm(list = ls()) 
load("DEG_PPI.Rdata")
load("GSE50772step4output.Rdata")
set.seed(5)

library(mlbench)
library(caret)

library(dplyr)
SLE_exp <- mutate(exp,symbol=rownames(exp))
k1 = SLE_exp$symbol%in% DEG_PPI;table(k1)
SLE_lasso= SLE_exp[k1,]
SLE_lasso<- SLE_lasso[,-c(82)]

meta <- read.csv("ClinicalTraits.csv")

x=t(SLE_lasso)
y=meta$SLE

control <- rfeControl(functions = caretFuncs, method = "cv", number = 10)

results <- rfe(x, 
               y, 
               sizes = c(1:69), 
               rfeControl = control,
               method = "svmRadial")

print(results)

predictors(results)

plot(results, type=c("g", "o"))
SLE_SVM<-predictors(results)

save(SLE_SVM,
     file = "SLE_SVM.Rdata") 

###SYM-RFE analysis of AF
rm(list = ls()) 
load("DEG_PPI.Rdata")
load("mergestep4output.Rdata")
set.seed(765544)

library(mlbench)
library(caret)

library(dplyr)
AF_exp <- mutate(exp,symbol=rownames(exp))
k1 = AF_exp$symbol%in% DEG_PPI;table(k1)
AF_lasso= AF_exp[k1,]
AF_lasso<- AF_lasso[,-c(49)]

meta <- read.csv("AF_ClinicalTraits.csv")

x=t(AF_lasso)
y=meta$AF

control <- rfeControl(functions = caretFuncs, method = "cv", number = 10)

results <- rfe(x, 
               y, 
               sizes = c(1:69), 
               rfeControl = control,
               method = "svmRadial")

print(results)

predictors(results)

plot(results, type=c("g", "o"))
AF_SVM<-predictors(results)

save(AF_SVM,
     file = "AF_SVM.Rdata") 


###intersection of results from Lasso and SYM-RFE analysis. 
rm(list = ls())
k2 = SLE_SVM %in% AF_SVM;table(k2)
AF_SLE_SYM= SLE_SVM[k2]

k3 = AF_SLE_SYM %in%  AF_SLE_Lasso;table(k2)
AF_lasso_SYM= AF_SLE_SYM[k2]

write.table(AF_lasso_SYM,file="SLE_SVM.txt",sep='\t',quote = F,row.names = F)
