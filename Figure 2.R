###AF-WGCNA
rm(list = ls())  
library(WGCNA)

options(stringsAsFactors=F)


load(file = "mergestep4output.Rdata")
femData<-degexp
library(dplyr)
femData<- mutate(femData,symbol=rownames(femData)) 

dim(femData)
names(femData)

datExpr0 <- as.data.frame(t(femData[, -c(49)]))
names(datExpr0) <- femData$symbol
rownames(datExpr0) <- names(femData)[-c(49)]

gsg <- goodSamplesGenes(datExpr0, verbose=3)

gsg$allOK

if (!gsg$allOK){
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse=", ")))
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse=", ")))
  datExpr0 <- datExpr0[gsg$goodSamples, gsg$goodGenes]
}


sampleTree <- hclust(dist(datExpr0), method="average")
sizeGrWindow(12, 9)
par(cex=0.6)
par(mar=c(0,4,2,0))
plot(sampleTree, main="Sample clustering to detect outliers", sub="", xlab="", cex.lab=1.5, 
     cex.axis=1.5, cex.main=2)
abline(h=25, col="red")


clust <- cutreeStatic(sampleTree, cutHeight=25, minSize=10)
table(clust)
keepSamples <- (clust==1)
datExpr <- datExpr0[keepSamples, ]
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)

datExpr <- datExpr0

traitData <- read.csv("ClinicalTraits.csv")
allTraits <- traitData[, c(1:3)]

dim(allTraits)
names(allTraits)


femaleSamples <- rownames(datExpr)
traitRows <- match(femaleSamples, allTraits$Sample)
datTraits <- allTraits[traitRows, -1]
rownames(datTraits) <- allTraits[traitRows, 1]


collectGarbage()


sampleTree2 <- hclust(dist(datExpr), method="average")
traitColors <- numbers2colors(datTraits, signed=F)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels=names(datTraits), 
                    main="Sample dendrogram and trait heatmap")


save(datExpr, datTraits, file="dataInput.RData")

allowWGCNAThreads()

enableWGCNAThreads(nThreads=16)

powers <- c(c(1:10), seq(from=12, to=20, by=2))

sft <- pickSoftThreshold(datExpr, powerVector=powers, verbose=5, networkType="unsigned")

sft$powerEstimate

sizeGrWindow(9, 5)
par(mfrow=c(1, 2))
cex1 <- 0.9
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3])*sft$fitIndices[, 2],
     xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit,signed R^2", type="n",
     main=paste("Scale independence"))
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3])*sft$fitIndices[, 2],
     labels=powers, cex=cex1, col="red")
abline(h=0.8,col="red")
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab="Soft Threshold (power)", ylab="Mean Connectivity", type="n",
     main=paste("Mean connectivity"))
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels=powers, cex=cex1, col="red")

sft$powerEstimate <- 14

net <- blockwiseModules(datExpr, corType="pearson",
                        power=sft$powerEstimate,
                        TOMType="unsigned", saveTOMs=TRUE, saveTOMFileBase="femaleMouseTOM",
                        deepSplit=4, minModuleSize=30,
                        reassignThreshold=0, mergeCutHeight=0.25,
                        numericLabels=T, pamRespectsDendro=F, nThreads=0,
                        verbose=3,maxBlockSize = 30000)

table(net$colors)

sizeGrWindow(50, 29)

mergedColors <- labels2colors(net$colors)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels=F, hang=0.03,
                    addGuide=T, guideHang=0.05)

moduleLabels <- net$colors
moduleColors <- labels2colors(net$colors)
MEs <- net$MEs
geneTree <- net$dendrograms[[1]]

save(net,MEs, moduleLabels, moduleColors, geneTree,
     file="networkConstruction-auto.RData")


nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)

MEs0 <- moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)

moduleTraitCor <- cor(MEs, datTraits, use="p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)


sizeGrWindow(10, 8)

textMatrix <-  paste(signif(moduleTraitCor, 2), "\n(",
                     signif(moduleTraitPvalue, 1), ")", sep="")
dim(textMatrix) <- dim(moduleTraitCor)
par(mar=c(4, 8.5, 3, 3))

labeledHeatmap(Matrix=moduleTraitCor,
               xLabels=names(datTraits),
               yLabels=names(MEs),
               ySymbols=names(MEs),
               colorLabels=F,
               colors=greenWhiteRed(50),
               textMatrix=textMatrix,
               setStdMargins=F,
               cex.text=1,
               zlim=c(-1,1),
               main=paste("Module-trait relationships"))


modNames <- substring(names(MEs), 3)

geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use="p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) <- paste("MM", modNames, sep="")
names(MMPvalue) <- paste("p.MM", modNames, sep="")

geneTraitSignificance <- as.data.frame(cor(datExpr, datTraits, use="p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) <- paste("GS.", names(datTraits), sep="")
names(GSPvalue) <- paste("p.GS.", names(datTraits), sep="")

geneInfo <- cbind(geneModuleMembership, MMPvalue, geneTraitSignificance, GSPvalue)
write.table(geneInfo, file="geneInfo.txt", sep="\t", quote=F)


TOM <- TOMsimilarityFromExpr(datExpr, power=sft$powerEstimate)
dissTOM <- 1 - TOM

plotTOM <- dissTOM^7
diag(plotTOM) <- NA
TOMplot(plotTOM, geneTree, moduleColors, main="Network heatmap plot, all genes")

nSelect <- 400
set.seed(10)
select <- sample(nGenes, size=nSelect)
selectTOM <- dissTOM[select, select]


selectTree <- hclust(as.dist(selectTOM), method="average")
selectColors <- moduleColors[select]

plotDiss <- selectTOM^7
diag(plotDiss) <- NA
TOMplot(plotDiss, selectTree, selectColors, main="Network heatmap plot, selected genes")

plotEigengeneNetworks(MEs, "", marDendro=c(0, 4, 1, 2),
                      marHeatmap=c(3, 4, 1, 2), cex.lab=0.8,
                      plotDendrograms=T, plotHeatmaps=T,
                      xLabelsAngle=90)


rm(list = ls())  
getwd()
workingDir = "."
library(WGCNA)
options(stringsAsFactors = FALSE)
load("networkConstruction-auto.RData")
load("dataInput.RData")

names(datExpr)
names(datExpr)[moduleColors=="turquoise"]
X<-(datExpr)[moduleColors=="turquoise"]
AF_turquoise<-colnames(X)

names(datExpr)[moduleColors=="grey"]
X<-(datExpr)[moduleColors=="grey"]
AF_grey<-colnames(X)
save(AF_grey,AF_turquoise, file = "AF_module.Rdata")



###SLE-WGCNA
rm(list = ls())  
library(WGCNA)

options(stringsAsFactors=F)

load(file = "mergestep4output.Rdata")
femData<-degexp
library(dplyr)
femData<- mutate(femData,symbol=rownames(femData)) 

dim(femData)
names(femData)

datExpr0 <- as.data.frame(t(femData[, -c(82)]))
names(datExpr0) <- femData$symbol
rownames(datExpr0) <- names(femData)[-c(82)]


gsg <- goodSamplesGenes(datExpr0, verbose=3)

gsg$allOK

if (!gsg$allOK){
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse=", ")))
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse=", ")))
  datExpr0 <- datExpr0[gsg$goodSamples, gsg$goodGenes]
}


sampleTree <- hclust(dist(datExpr0), method="average")
sizeGrWindow(12, 9)
par(cex=0.6)
par(mar=c(0,4,2,0))
plot(sampleTree, main="Sample clustering to detect outliers", sub="", xlab="", cex.lab=1.5, 
     cex.axis=1.5, cex.main=2)
abline(h=300, col="red")

clust <- cutreeStatic(sampleTree, cutHeight=15, minSize=10)
table(clust)
keepSamples <- (clust==1)
datExpr <- datExpr0[keepSamples, ]
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)

datExpr <- datExpr0

traitData <- read.csv("ClinicalTraits.csv")
allTraits <- traitData[, c(1:3)]

dim(allTraits)
names(allTraits)


femaleSamples <- rownames(datExpr)
traitRows <- match(femaleSamples, allTraits$Sample)
datTraits <- allTraits[traitRows, -1]
rownames(datTraits) <- allTraits[traitRows, 1]

collectGarbage()

sampleTree2 <- hclust(dist(datExpr), method="average")
traitColors <- numbers2colors(datTraits, signed=F)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels=names(datTraits), 
                    main="Sample dendrogram and trait heatmap")

save(datExpr, datTraits, file="dataInput.RData")

allowWGCNAThreads()

enableWGCNAThreads(nThreads=16)


powers <- c(c(1:10), seq(from=12, to=20, by=2))

sft <- pickSoftThreshold(datExpr, powerVector=powers, verbose=5, networkType="unsigned")

sft$powerEstimate

sizeGrWindow(9, 5)
par(mfrow=c(1, 2))
cex1 <- 0.9
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3])*sft$fitIndices[, 2],
     xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit,signed R^2", type="n",
     main=paste("Scale independence"))
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3])*sft$fitIndices[, 2],
     labels=powers, cex=cex1, col="red")
abline(h=0.7, col="red")
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab="Soft Threshold (power)", ylab="Mean Connectivity", type="n",
     main=paste("Mean connectivity"))
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels=powers, cex=cex1, col="red")

sft$powerEstimate <- 8

net <- blockwiseModules(datExpr, corType="pearson",
                        power=sft$powerEstimate,
                        TOMType="unsigned", saveTOMs=TRUE, saveTOMFileBase="femaleMouseTOM",
                        deepSplit=4, minModuleSize=30,
                        reassignThreshold=0, mergeCutHeight=0.25,
                        numericLabels=T, pamRespectsDendro=F, nThreads=0,
                        verbose=3,maxBlockSize = 30000)


table(net$colors)

sizeGrWindow(50, 29)

mergedColors <- labels2colors(net$colors)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels=F, hang=0.03,
                    addGuide=T, guideHang=0.05)

moduleLabels <- net$colors
moduleColors <- labels2colors(net$colors)
MEs <- net$MEs
geneTree <- net$dendrograms[[1]]

save(MEs, moduleLabels, moduleColors, geneTree,
     file="networkConstruction-auto.RData")

nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)


MEs0 <- moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)

moduleTraitCor <- cor(MEs, datTraits, use="p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

sizeGrWindow(10, 8)

textMatrix <-  paste(signif(moduleTraitCor, 2), "\n(",
                     signif(moduleTraitPvalue, 1), ")", sep="")
dim(textMatrix) <- dim(moduleTraitCor)
par(mar=c(4, 8.5, 3, 3))
labeledHeatmap(Matrix=moduleTraitCor,
               xLabels=names(datTraits),
               yLabels=names(MEs),
               ySymbols=names(MEs),
               colorLabels=F,
               colors=greenWhiteRed(50),
               textMatrix=textMatrix,
               setStdMargins=F,
               cex.text=1,
               zlim=c(-1,1),
               main=paste("Module-trait relationships"))


modNames <- substring(names(MEs), 3)

geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use="p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) <- paste("MM", modNames, sep="")
names(MMPvalue) <- paste("p.MM", modNames, sep="")

geneTraitSignificance <- as.data.frame(cor(datExpr, datTraits, use="p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) <- paste("GS.", names(datTraits), sep="")
names(GSPvalue) <- paste("p.GS.", names(datTraits), sep="")

geneInfo <- cbind(geneModuleMembership, MMPvalue, geneTraitSignificance, GSPvalue)
write.table(geneInfo, file="geneInfo.txt", sep="\t", quote=F)

TOM <- TOMsimilarityFromExpr(datExpr, power=sft$powerEstimate)
dissTOM <- 1 - TOM

plotTOM <- dissTOM^7
diag(plotTOM) <- NA
TOMplot(plotTOM, geneTree, moduleColors, main="Network heatmap plot, all genes")



nSelect <- 400
set.seed(10)
select <- sample(nGenes, size=nSelect)
selectTOM <- dissTOM[select, select]

selectTree <- hclust(as.dist(selectTOM), method="average")
selectColors <- moduleColors[select]
plotDiss <- selectTOM^7
diag(plotDiss) <- NA
TOMplot(plotDiss, selectTree, selectColors, main="Network heatmap plot, selected genes")

plotEigengeneNetworks(MEs, "", marDendro=c(0, 4, 1, 2),
                      marHeatmap=c(3, 4, 1, 2), cex.lab=0.8,
                      plotDendrograms=T, plotHeatmaps=T,
                      xLabelsAngle=90)

rm(list = ls())  
getwd()
workingDir = "."
library(WGCNA)
options(stringsAsFactors = FALSE)
load("networkConstruction-auto.RData")
load("dataInput.RData")

names(datExpr)
names(datExpr)[moduleColors=="brown"]
X<-(datExpr)[moduleColors=="brown"]
SLE_brown<-colnames(X)

names(datExpr)[moduleColors=="blue"]
X<-(datExpr)[moduleColors=="blue"]
SLE_blue<-colnames(X)

names(datExpr)[moduleColors=="turquoise"]
X<-(datExpr)[moduleColors=="turquoise"]
SLE_turquoise<-colnames(X)

names(datExpr)[moduleColors=="grey"]
X<-(datExpr)[moduleColors=="grey"]
SLE_grey<-colnames(X)

save(SLE_blue, SLE_brown,SLE_grey,SLE_turquoise, file = "SLE_module.Rdata")

