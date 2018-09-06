setwd("/geschwindlabshares/CrossDisorder_transcriptome_comparison/ArgMouse/test3/")
library(WGCNA)
options(stringsAsFactors = FALSE)

bsize = 5000
powers = c(seq(1,9,by=1),seq(10,30,by=2))
enableWGCNAThreads()
allowWGCNAThreads()

load("datExprForBlckMod_CBPure.RData")
datExpr = t(tdatExpr)

pdf("softThresh_CBPure.pdf")  ## *do on the server

sft = pickSoftThreshold(data= tdatExpr, networkType = "signed", corFnc="bicor",verbose=5,powerVector=powers,blockSize = bsize)

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Thresh Power", ylab="Scale free R^2",type="n")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels = powers, cex = 0.7, col="red",  xlab="Soft Thresh Power", ylab="Scale free R^2",main=names(datExpr))
abline(h=0.8, col="black")
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab = "Soft threshold power", ylab = "Mean connectivity", type = "n")
text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, cex = 0.7, col="black")

dev.off()

## Choose 18 again
## Get TOM and blockwise mods

rm(list=ls())
load("datExprForBlckMod_CBPure.rda")

net = blockwiseModules(datExpr=tdatExpr, maxBlockSize = 15000,
                       power = 18, TOMType = "signed", minModuleSize = 200,
                       reassignThreshold = 0, mergeCutHeight = 0.2,
                       numericLabels = TRUE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "ArgMs_CBPure_TOM",
                       verbose = 3, deepSplit = 4, pamStage = FALSE, networkType = "signed",
                       corType = "bicor", nThreads = 24)

save(net,file="blckMod_CBPure.rda")

## See what happened!!

load("blckMod_CBPure.rda")
load("ArgMs_CBPure_TOM-block.1.RData")

consMEs = net$multiMEs;
moduleLabels = net$colors;
moduleColors = (moduleLabels)
consTree = net$dendrograms[[1]];

pdf(file="ArgMs_PrelimDendro_CBPure.pdf")

plotDendroAndColors(consTree, moduleColors,
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "dendrogram and module colors - All ArgMs")
dev.off()


## Recut Tree to see if that will help with module identification

geneTree = hclust(1-TOM, method = "average")
colors = vector(mode="list")
labels = vector(mode="list")

load("datExprForBlckMod_CBPure.RData") ## nrow=samples ncol=genes
#load("ArgMs_TOM-block.1.RData")

for (pam in c(FALSE)) {
  for (minModSize in c(50,100, 200)) {
    for (dthresh in c(0.1, 0.2)) {
      for(ds in c(1,2,3,4)) { 
        print(paste("DS=", ds, ",MMS=", minModSize, ",DCOR=",dthresh,",PAM=",pam,sep=""))
        tree = cutreeHybrid(dendro = geneTree, minClusterSize= minModSize, pamStage=pam, cutHeight = 0.999, deepSplit=ds, distM=as.matrix(1-TOM))
        merged = mergeCloseModules(exprData= tdatExpr, colors = tree$labels, cutHeight=dthresh)
        colors = cbind(colors, labels2colors(merged$colors))
        
        labels = c(labels, paste("DS=", ds, ",MMS=", minModSize, ",DCOR=",dthresh,sep=""))
      }
    }
  }
}

pdf("ArgMsDendrogram_DiffParams_CBPure.pdf")
plotDendroAndColors(geneTree, colors, groupLabels=labels, addGuide= TRUE, 
                    dendroLabels=FALSE, main="Dendrogram", cex.colorLabels=0.5)
dev.off()

## ds=2 , mms=200, dcor=0.1      choose for test3 and test4


## Tree recut with chosen parameters

#load("datExprForBlckMod.rda") ## nrow=samples ncol=genes
#load("Liver_TOM-block.1.RData")

#geneTree = hclust(1-TOM, method = "average")

## Match new modules to old modules

#newColors = colors
#newMerged = merged
#rm(colors,merged)

#load("ArgMsCols.RData")
#orgColors = colors
#rm(colors,merged)
colors = vector(mode="list")

minModSize = 200
pamStage = F
ds = 2                
dthresh = 0.1        

tree = cutreeHybrid(dendro = geneTree, minClusterSize= minModSize, 
                    pamStage=pamStage, cutHeight = 0.999, deepSplit=ds, distM=as.matrix(1-TOM)); 
merged = mergeCloseModules(exprData= tdatExpr, colors = tree$labels, cutHeight=dthresh)
#consColors <- matchLabels(merged$colors,orgColors)
#colors = labels2colors(consColors)
colors = labels2colors(merged$colors)
labels = paste("DS=", ds, ",MMS=", minModSize, ",DCOR=",dthresh,sep="")
table(colors)
length(table(colors))

pdf(file="ArgMsDendrogram_ReCut_CBPure_test4.pdf")
plotDendroAndColors(geneTree, colors, groupLabels=(labels), 
                    addGuide= TRUE, dendroLabels=FALSE, main="Dendrogram", cex.colorLabels=0.4)
dev.off()

save(colors,merged,file="ArgMsCols_CBPure_test4.RData")







