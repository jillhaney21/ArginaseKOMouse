## ArgMouse.R
## Script to analyze KO, WT, Treated KO, and Het arginase mice

options(stringsAsFactors = FALSE)
setwd("C:/Users/Jill/Dropbox/DHGLab/ArgMouse/")

library(affy)

data.affy = ReadAffy(celfile.path = "./rawDat/")

GSM = rownames(pData(data.affy))
GSM = substr(GSM,1,7)
GSM[7:9]=substr(rownames(pData(data.affy))[7:9],1,14)

datExpr = log2(exprs(data.affy))
dim(datExpr)

datMeta = cbind(GSM,substr(GSM,1,2))
rownames(datMeta) = datMeta[,1]
datMeta = cbind(datMeta,rep("Ms",6))
colnames(datMeta) = c("name","dx","species")
datMeta = as.data.frame(datMeta)

pdf(file = "ArgMouse_rawStats.pdf")

#boxplot

boxplot(datExpr,range=0, col = as.numeric(as.factor(datMeta$dx)), xaxt='n', xlab = "Array", main = "Boxplot Pre-Normalization", ylab = "Intensity")
legend("topright",legend=levels(as.factor(datMeta$dx)),fill=as.numeric(as.factor((levels(as.factor(datMeta$dx))))))

#histogram

i=1
plot(density((datExpr[,i]),na.rm=T),col = as.numeric(as.factor(datMeta[,2]))[i],
     main = "Hist of Log2 Exp", xlab="log2 exp",ylim=c(0,0.38))
for(i in 2:12){
  lines(density((datExpr[,i]),na.rm=T), col = as.numeric(as.factor(datMeta[,2]))[i],)
}
legend("topright",legend=levels(as.factor(datMeta$dx)),fill=as.numeric(as.factor((levels(as.factor(datMeta$dx))))))

#black = control, red = parkinson's

#mdsplot

mds = cmdscale(dist(t(datExpr)),eig=TRUE)
plot(mds$points,col=as.numeric(as.factor(datMeta$dx)),pch=19)
legend("bottomleft",legend=levels(as.factor(datMeta$dx)),fill=as.numeric(as.factor((levels(as.factor(datMeta$dx))))))


##Normalization
datExpr = rma(data.affy, background=T, normalize=T, verbose=T)
datExpr = exprs(datExpr)


## Get Batch
batch = protocolData(data.affy)$ScanDate
batch = substr(batch,1,10)
batch = as.factor(batch)
table(batch)
datMeta$Batch = batch

mds = cmdscale(dist(t(log2(exprs(data.affy)))),eig=TRUE)
plot(mds$points,col=as.numeric(as.factor(datMeta$Batch)),pch=19)
legend("bottom",legend=levels(as.factor(datMeta$Batch)),fill=as.numeric(as.factor((levels(as.factor(datMeta$Batch))))))

dev.off()

## Batch Correction 

#source("http://bioconductor.org/biocLite.R")
##biocLite("sva")
library(sva)

## there is a confound of Dx and batch - just have to let it be and hope for the best
## can't model Dx effects bc of the confound, so just give ComBat the datExpr and batch vector

batch = as.factor(datMeta$Batch)
datExpr.combat = ComBat(dat = datExpr,batch = batch)

datExpr_preComBat = datExpr
datExpr = datExpr.combat


## OUtlier Removal
tree = hclust(dist(t(datExpr)), method="average")
plot(tree)

##biocLite("WGCNA")
library(WGCNA)
normadj = (0.5 + 0.5*bicor(datExpr))^2
netsummary = fundamentalNetworkConcepts(normadj)
C = netsummary$Connectivity
Z.C = (C-mean(C))/sqrt(var(C))
to_keep = abs(Z.C) < 2

## all good!

## QC after normalization -- repeat all the same steps as above

dim(datExpr)

setwd("./final/")
pdf(file = "ArgMouse_rma_cbPure.pdf")

#boxplot

boxplot(datExpr,range=0, col = as.numeric(as.factor(datMeta$dx)), xaxt='n', xlab = "Array", main = "Boxplot Normalized", ylab = "Intensity")
legend("topright",legend=levels(as.factor(datMeta$dx)),fill=as.numeric(as.factor((levels(as.factor(datMeta$dx))))))


#histogram

i=1
plot(density((datExpr[,i]),na.rm=T),col = as.numeric(as.factor(datMeta$dx[i])),main = "Hist of Log2 Exp", xlab="log2 exp")
for(i in 2:12){
  lines(density((datExpr[,i]),na.rm=T), col = as.numeric(as.factor(datMeta$dx))[i])
}
legend("topright",legend=levels(as.factor(datMeta$dx)),fill=as.numeric(as.factor((levels(as.factor(datMeta$dx))))))


#mdsplot

mds = cmdscale(dist(t(datExpr)),eig=TRUE)
plot(mds$points,col=as.numeric(as.factor(datMeta$dx)),pch=19)
legend("bottomleft",legend=levels(as.factor(datMeta$dx)),fill=as.numeric(as.factor((levels(as.factor(datMeta$dx))))))
plot(mds$points,col=as.numeric(as.factor(datMeta$Batch)),pch=19)
legend("bottomleft",legend=levels(as.factor(datMeta$Batch)),fill=as.numeric(as.factor((levels(as.factor(datMeta$Batch))))))

dev.off()

datExpr = as.data.frame(datExpr)
#normArg1_1 = datExpr1[which(rownames(datExpr1) == "1419549_at"),]

#annotating probes

#biocLite("biomaRt")
library(biomaRt)

ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL",dataset="mmusculus_gene_ensembl",host = "jul2015.archive.ensembl.org")

f = listFilters(ensembl)
a = listAttributes(ensembl)

identifier <- "affy_mouse430_2"
getinfo <- c("affy_mouse430_2", "ensembl_gene_id", "entrezgene", "external_gene_name")
geneDat <- getBM(attributes = getinfo, filters=identifier, values = rownames(datExpr),mart=ensembl)

## If you want to examine probes before filtering:
##ENSMUSG00000019987 = Arg1

idx = match(rownames(datExpr),geneDat$affy_mouse430_2)
geneDat_filt = geneDat[idx,]

table(is.na(geneDat_filt$ensembl_gene_id))

to_keep = (is.na(geneDat_filt$ensembl_gene_id) == FALSE)
geneDat_filt = geneDat_filt[to_keep,]
datExpr_org = datExpr
datExpr = datExpr[to_keep,]

#collapse rows

library(WGCNA)

CR = collapseRows(datExpr, rowGroup = geneDat_filt$ensembl_gene_id, rowID = geneDat_filt$affy_mouse430_2)
datExpr = CR$datETcollapsed
datExpr = as.data.frame(datExpr)
idx = match(CR$group2row[,"selectedRowID"], geneDat_filt$affy_mouse430_2)
geneDat_filt = geneDat_filt[idx,]
rownames(geneDat_filt) = geneDat_filt$ensembl_gene_id

dim(datExpr)
dim(geneDat_filt)
dim(datMeta)

write.csv(datExpr,file="DatExpr_CBpure.csv")
write.csv(datMeta, file = "DatMeta_CBpure.csv")
write.csv(geneDat_filt, file = "geneDat_CBpure.csv")

#differential expression

options(stringsAsFactors = FALSE)

datExpr = read.csv(file="DatExpr_CBpure.csv",row.names=1)
datMeta = read.csv(file="DatMeta_CBpure.csv",row.names=1)
geneDat_filt = read.csv(file="geneDat_CBpure.csv",row.names=1)

## Note: OPTIONAL, try different gene filtering steps to remove low variance genes from analysis
## test3 -> try filtering genes to draw out variance from more conservatively corrected data, remove the lowest 5% sd genes
## test4 -> remove the lowest 25% sd genes

pres = apply(datExpr,1,sd) 
idx = which(pres > quantile(apply(datExpr,1,sd),0.25)) ## standard deviation greater than 0.05 for test three, sd greater than 0.077 for test 4
datExpr_unfilt = datExpr
datExpr = datExpr[idx,]

setwd("./test4/")

save(datExpr,file="datExpr_filt.csv")
datExpr = read.csv(file="datExpr_filt.csv",row.names=1)

#biocLite("limma")
library(limma)

datMeta$Dx = factor(datMeta$dx, levels= c("WT", "KO","HE","TR"))
datMeta$Age = rep(c(13,14,15),4)
mod = model.matrix(~datMeta$Dx) ##+datMeta$Age)
fit = lmFit(datExpr,mod)
fit = eBayes(fit)

idx = match(rownames(datExpr),rownames(geneDat_filt))
geneDat_file = geneDat_filt[idx,]

tt_KO = topTable(fit,coef = 2,n = Inf,genelist = geneDat_file)
tt_HE = topTable(fit,coef = 3,n = Inf,genelist = geneDat_file)
tt_TR = topTable(fit,coef = 4,n = Inf,genelist = geneDat_file)

save(tt_KO,tt_HE,tt_TR,file="ttables_CBpure.RData")
load("ttables_CBpure.RData")

### do WGCNA on the first dataset

#rm(list=ls())
library(WGCNA)
options(stringsAsFactors =F)

## First, Calculate Soft-Threshold

bsize = 5000
powers = c(seq(1,9,by=1),seq(10,30,by=2))
enableWGCNAThreads()
allowWGCNAThreads()

tdatExpr = t(datExpr)
save(tdatExpr,file="datExprForBlckMod_CBPure.RData")

pdf("softThresh_CBPure.pdf")  ## *do on the server

sft = pickSoftThreshold(data= tdatExpr, networkType = "signed", corFnc="bicor",verbose=5,powerVector=powers,blockSize = bsize)

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Thresh Power", ylab="Scale free R^2",type="n")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels = powers, cex = 0.7, col="red",  xlab="Soft Thresh Power", ylab="Scale free R^2",main=names(datExpr))
abline(h=0.8, col="black")
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab = "Soft threshold power", ylab = "Mean connectivity", type = "n")
text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, cex = 0.7, col="black")

dev.off()

## Choose soft-threshold of 18 - very high, but we didn't have that many samples

## Time for blockwiseModules - do this on the server, since it's such a large dataset - see WGCNA.R on Orion

## Now, have colors for modules

rm(list=ls())
setwd("C:/Users/jillh/Dropbox/DHGLab/ArgMouse/test3/")
library(WGCNA)
load(file = "ArgMsCols_CBPure_test3.RData")  ## using test3 variability filter for WGCNA and DGE
datExpr = read.csv(file = "datExpr_filt.csv",row.names=1)
datMeta = read.csv(file = "C:/Users/Jill/Dropbox/DHGLab/ArgMouse/proc/datMeta.csv")
geneDat = read.csv("C:/Users/Jill/Dropbox/DHGLab/ArgMouse/test/geneDat_CBpure.csv",row.names=1)
idx = match(rownames(datExpr),geneDat$ensembl_gene_id)
geneDat = geneDat[idx,]

mods = cbind(rownames(datExpr),colors)
mods[grep("ENSMUSG00000019987",mods[,1]),]
## Arg1 in pink
idx=match(geneDat$ensembl_gene_id,mods[,1])
geneDat = cbind(geneDat,mods[idx,])

write.csv(geneDat,file="geneDat_CBPure_wMods.csv")
geneDat = read.csv("geneDat_CBPure_wMods.csv")

MEs = merged$newMEs
oldNames = colnames(MEs)
oldNames = as.numeric(gsub("^ME","",oldNames))
modNames = labels2colors(oldNames)
colnames(MEs) = modNames
colnames(mods) = c("ENSG","Module")

save(MEs,mods,file="MEs_CBPure.RData")
load("MEs_CBPure.RData")

# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
#Order MEs
MEs = orderMEs(MEs)


# Now background list

background = rownames(datExpr)

# Now module lists

modules = list()
for(i in 1:length(unique(mods[,2]))){
  modules[[i]] = subset(background,colors==modNames[i])
  names(modules)[i] = paste("MM",modNames[i],sep="")
}

save(modules,background,file="ArgMouse_Analysis_CBPure.RData")
load("ArgMouse_Analysis_CBPure.RData")

## Now GOElite

for (i in 1:18){
  write.table(modules[[i]],file=paste("./Modules/",names(modules)[i],".txt",sep=""),row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")
}

write.table(background,file="./Background/background.txt",row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")

### Run GOElite on PC (download the software online)

####Plotting the GO Output

uniquemodcolors = unique(colors)
uniquemodcolors=uniquemodcolors[-which(uniquemodcolors=="grey")]

pdf("GOElite_plot_Modules_altFilt.pdf",height=8,width=8)

for(i in 1:length(uniquemodcolors)){
  thismod= uniquemodcolors[i]
  nameThisMod = paste("MM",uniquemodcolors[i],sep="")
  
  tmp=read.csv(file=paste("./GOelite/GO-Elite_results/CompleteResults/ORA_pruned/",nameThisMod,"-GO_z-score_elite.txt",sep=""),sep="\t")
  tmp=subset(tmp,Ontology.Type=="biological_process")
  tmp=tmp[,c(2,4,7,9)] ## Select GO-terms and Z-score
  tmp=tmp[order(tmp$Z.Score,decreasing=T),] #
  
  if(length(which(tmp[,2] >= 5 | tmp[,3] >= 10)) != 0 ){         ## before just number changed >= 5 filtering
    tmp=tmp[which(tmp[,2] >= 5 | tmp[,3] >= 10),]
  } else
    next
  
  if (nrow(tmp)<10){
    tmp1=tmp ## Take top 10 Z-score
    tmp1 = tmp1[order(tmp1$Z.Score),] ##Re-arrange by increasing Z-score
    par(mar=c(4,22,4,4))
    barplot(tmp1$Z.Score,horiz=T,col="blue",names.arg= tmp1$Ontology.Name,cex.names=1.2,las=1,main=paste("Gene Ontology Plot of",thismod,"Module"),xlab="Z-Score")
    abline(v=2,col="red")
    
  } else
  
  tmp1=tmp[c(1:10),] ## Take top 10 Z-score
  tmp1 = tmp1[order(tmp1$Z.Score),] ##Re-arrange by increasing Z-score
  par(mar=c(4,22,4,4))
  barplot(tmp1$Z.Score,horiz=T,col="blue",names.arg= tmp1$Ontology.Name,cex.names=1.2,las=1,main=paste("Gene Ontology Plot of",thismod,"Module"),xlab="Z-Score")
  abline(v=2,col="red")
  
  cat('Done ...',thismod,'\n')
  
}

dev.off()


### Try out gProfileR as alternate annotation program

biocLite("gProfileR")
library(gProfileR)

for (mod in unique(colors)){
  print(mod)
  temp = gprofiler(modules[[paste("MM",mod,sep="")]],organism="mmusculus",ordered_query=T,correction_method="bonferroni",max_p_value = 0.05,max_set_size=500)
  print(paste(mod,"gProfileR complete!",sep=" "))
  temp = temp[order(temp$p.value),]
  BP = temp[which(temp$domain == "BP"),]
  write.csv(BP,file=paste("./gProfileR/BP/",mod,".csv"))
  write.csv(temp,file=paste("./gProfileR/all/",mod,".csv"))
}

## Look more closely at network dynamics, for KO and WT only

#datExpr_KOWT = datExpr[,-which(datMeta$dx=="HE" | datMeta$dx=="TR")]
#MEs_KOWT = MEs[-which(datMeta$dx=="HE" | datMeta$dx=="TR"),]

#kME = signedKME(t(datExpr_KOWT),MEs_KOWT,corFnc="bicor")  ## bicor of each gene with module eigengene, comparing WT and KO only
#kMEPvalue = as.data.frame(corPvalueStudent(as.matrix(kME), 6))

kME = signedKME(t(datExpr),MEs,corFnc="bicor")  ## bicor of each gene with module eigengene
kMEPvalue = as.data.frame(corPvalueStudent(as.matrix(kME), 12))
kIN = intramodularConnectivity.fromExpr(t(datExpr),colors,scaleByMax=TRUE,
                                                  corFnc="bicor",networkType="signed",power=18,getWholeNetworkConnectivity = FALSE)
colnames(kME) = colnames(kMEPvalue) = paste("kME",colnames(MEs),sep="")

save(kME,kMEPvalue,kIN,file="ModCharacteristics_CBPure.RData")
#save(kME,kMEPvalue,kIN,file="ModCharacteristics_KOWT.RData")
load("ModCharacteristics_CBPure.RData")

colors[grep("ENSMUSG00000019987",rownames(datExpr))]  ## Arg1 in brown for test3

## Compare module dysregulation of KO, HE, and TR

geneTraitSignificance = list()
GSPvalue = list()

pdf(file = "WTComparisons.pdf")
for (pheno in c("KO","HE","TR")){

dxWGCNA2 = datMeta$dx[which(datMeta$dx == "WT" | datMeta$dx == pheno)]
dxWGCNA2 = gsub("WT",0,dxWGCNA2)
dxWGCNA2 = gsub(pheno,1,dxWGCNA2)
dxWGCNA2 = as.factor(dxWGCNA2)

moduleTraitCor2 = cor(MEs[which(datMeta$dx == "WT" | datMeta$dx == pheno),],as.numeric(as.factor(dxWGCNA2)))
moduleTraitCor2_s = cor(MEs[which(datMeta$dx == "WT" | datMeta$dx == pheno),],as.numeric(as.factor(dxWGCNA2)),method="spearman")
moduleTraitPvalue2 = corPvalueStudent(moduleTraitCor2,6)
moduleTraitPvalue2_s = corPvalueStudent(moduleTraitCor2_s,6)

#sizeGrWindow(10,6)
# Will display correlations and their p-values - Pearson
textMatrix = paste(signif(moduleTraitCor2, 2), " (",
                   signif(moduleTraitPvalue2, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor2)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
#pdf(file=paste("ModuleHM_",pheno,"vWT.pdf",sep=""))
labeledHeatmap(Matrix = moduleTraitCor2,
               xLabels = "Dx",
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships Pearson \n where ",pheno," = 1 > WT = 0, no TR/HE",sep=""))

# Will display correlations and their p-values - Spearman
textMatrix_s = paste(signif(moduleTraitCor2_s, 2), " (",
                   signif(moduleTraitPvalue2_s, 1), ")", sep = "");
dim(textMatrix_s) = dim(moduleTraitCor2_s)
par(mar = c(2, 6, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor2_s,
               xLabels = "Dx",
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix_s,
               setStdMargins = FALSE,
               cex.text = 0.75,
               zlim = c(-1,1),
               main = paste("Module-trait relationships Spearman \n where ",pheno," = 1 > WT = 0, no TR/HE",sep=""))

geneTraitSignificance[[pheno]] = as.data.frame(cor(t(datExpr[,which(datMeta$dx == "WT" | datMeta$dx == pheno)]),(as.numeric(dxWGCNA2)), use = "p"));
GSPvalue[[pheno]] = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance[[pheno]]), 6));
names(geneTraitSignificance[[pheno]]) = "GS.dx";
names(GSPvalue[[pheno]]) = "p.GS.dx";

}

dev.off()
## For Pure phenotype (KO v. WT only) - test3
## Negative Significant Module Eigengene Correlation with KO - yellow, midnightblue, tan, salmon, red
## Positive Significant Module Eigengene Correlation with KO - brown, magenta, pink, lightyellow, cyan

## Modules of interest to investigate, (KO and WT only):
pdf(file="ModsPositiveCorrWithKO.pdf")

for(i in c("brown","magenta","pink","lightyellow","cyan")){

if (i != "lightyellow"){
  
module = i
column = match(module, colnames(MEs));
moduleGenes = colors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(kME[moduleGenes, column]),
                   abs(geneTraitSignificance[["KO"]][moduleGenes, 1]),
                   xlab = paste("Module Membership (kME) in", module, "module"),
                   ylab = "Gene significance for Dx",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
}
else{
  module = i
  column = match(module, colnames(MEs));
  moduleGenes = colors==module;
  sizeGrWindow(7, 7);
  par(mfrow = c(1,1));
  verboseScatterplot(abs(kME[moduleGenes, column]),
                     abs(geneTraitSignificance[["KO"]][moduleGenes, 1]),
                     xlab = paste("Module Membership (kME) in", module, "module"),
                     ylab = "Gene significance for Dx",
                     main = paste("Module membership vs. gene significance\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "gold")
}
}
dev.off()

pdf(file="ModsNegativeCorrWithKO.pdf")

for(i in c("yellow","midnightblue","tan","salmon","red")){

if (i != "yellow"){
    module = i
    column = match(module, colnames(MEs));
    moduleGenes = colors==module;
    sizeGrWindow(7, 7);
    par(mfrow = c(1,1));
    verboseScatterplot(abs(kME[moduleGenes, column]),
                       abs(geneTraitSignificance[["KO"]][moduleGenes, 1]),
                       xlab = paste("Module Membership (kME) in", module, "module"),
                       ylab = "Gene significance for Dx",
                       main = paste("Module membership vs. gene significance\n"),
                       cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = i)
}

else{
  module = i
  column = match(module, colnames(MEs));
  moduleGenes = colors==module;
  sizeGrWindow(7, 7);
  par(mfrow = c(1,1));
  verboseScatterplot(abs(kME[moduleGenes, column]),
                     abs(geneTraitSignificance[["KO"]][moduleGenes, 1]),
                     xlab = paste("Module Membership (kME) in", module, "module"),
                     ylab = "Gene significance for Dx",
                     main = paste("Module membership vs. gene significance\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "yellow2")
}

}
dev.off()

geneTraitSigAll = cor()

## Make meaningful lists of interesting modules

intMods = c("brown","magenta","pink","lightyellow","cyan","yellow","midnightblue","tan","salmon","red")

mods_kME_GS=list()
for(i in intMods){
  print(i)
  idx = which(colors == i)
  
  mods_kME_GS[[i]]=list(data = as.data.frame(sort(kME[idx,paste("kME",i,sep="")],decreasing=T,index.return=TRUE)))
  
  names = rownames(kME)[idx]
  names = names[mods_kME_GS[[i]]$data[,2]]
  mods_kME_GS[[i]]$data = as.data.frame(mods_kME_GS[[i]]$data[,-2])
  rownames(mods_kME_GS[[i]]$data) = names
  
  idxP = match(rownames(mods_kME_GS[[i]]$data),rownames(kMEPvalue))
  mods_kME_GS[[i]]$data = cbind(mods_kME_GS[[i]]$data,kMEPvalue[idxP,paste("kME",i,sep="")])  
  
  GSK = geneTraitSignificance[["KO"]][match(rownames(mods_kME_GS[[i]]$data),rownames(geneTraitSignificance[["KO"]])),]
  mods_kME_GS[[i]]$data = cbind(mods_kME_GS[[i]]$data,GSK)
  
  GSKP = GSPvalue[["KO"]][match(rownames(mods_kME_GS[[i]]$data),rownames(GSPvalue[["KO"]])),]
  mods_kME_GS[[i]]$data = cbind(mods_kME_GS[[i]]$data,GSKP)
  
  GST = geneTraitSignificance[["TR"]][match(rownames(mods_kME_GS[[i]]$data),rownames(geneTraitSignificance[["TR"]])),]
  mods_kME_GS[[i]]$data = cbind(mods_kME_GS[[i]]$data,GST)
  
  GSTP = GSPvalue[["TR"]][match(rownames(mods_kME_GS[[i]]$data),rownames(GSPvalue[["TR"]])),]
  mods_kME_GS[[i]]$data = cbind(mods_kME_GS[[i]]$data,GSTP)
  
  GSH = geneTraitSignificance[["HE"]][match(rownames(mods_kME_GS[[i]]$data),rownames(geneTraitSignificance[["HE"]])),]
  mods_kME_GS[[i]]$data = cbind(mods_kME_GS[[i]]$data,GSH)
  
  GSHP = GSPvalue[["HE"]][match(rownames(mods_kME_GS[[i]]$data),rownames(GSPvalue[["HE"]])),]
  mods_kME_GS[[i]]$data = cbind(mods_kME_GS[[i]]$data,GSHP)
  
  kINMod = kIN[match(rownames(mods_kME_GS[[i]]$data),rownames(datExpr))]
  mods_kME_GS[[i]]$data = cbind(mods_kME_GS[[i]]$data,kINMod)
  
  colnames(mods_kME_GS[[i]]$data) = c("kME","kME_p","GS_KO","GS_KO_p","GS_TR","GS_TR_p","GS_HE","GS_HE_p","kIN")

  write.csv(mods_kME_GS[[i]]$data,file=paste("./intMods/",i,"_kME_GS.csv"))
}


############ ANOVA

plotAnova = function(mod, nm,ctype) {
  s = summary(mod)
  a = anova(mod)
  anova.pval = a$"Pr(>F)"[grep("est",rownames(a))]
  title = paste(nm,ctype, ": p=",signif(anova.pval,2),sep=" ")
  
  data=t(as.data.frame(s$coefficients[-1,]))
  #data=data[grep("Dx",rownames(data)),]
  #data$Dx = substring(rownames(data),11)
  colnames(data) = c("Beta", "SEM", "T", "p") #, "Group")
  col = rep("black",times = nrow(data))
  col[data[,4]<0.05] = "red"
  rownames(data) = "KO"
  
  ylim = 1.2*c((data[,1]-2*data[,2]),(data[,1]+2*data[,2]))
  
  bp = barplot(data[,1],names.arg = rownames(data), main=title, xlab="", ylab="beta",col = col) #,ylim=ylim)
  segments(bp, data[,1]-data[,2], bp, data[,1]+data[,2], lwd=1)
  segments(bp-0.1, data[,1]+data[,2], bp+0.1, data[,1]+data[,2], lwd=1); segments(bp-0.1, data[,1]-data[,2], bp+0.1, data[,1]-data[,2], lwd=1)
  
  #   bp = ggplot(data, aes(x=Group, y=Beta, fill=col)) +
  #     geom_bar(stat="identity") +   
  #     scale_fill_manual(values=c("black", "red")) +  
  #     geom_errorbar(aes(ymin=(Beta - SEM), ymax=(Beta + SEM)), 
  #                   position=position_dodge(width=0.8), width=0.25,size=0.25) +   
  #     ggtitle(title) +   
  #     theme(plot.title = element_text(size=20, face="bold", vjust=2)) +   
  #     labs(x="", y="beta") +   
  #     theme(axis.text.x=element_text(angle=50, size=10, hjust=1)) +   
  #     theme(axis.text.y=element_text(size=10, vjust=0.5)) +   
  #     theme(legend.position="none",   
  #           plot.title = element_text(size=20, face="bold"),     
  #           axis.title.x = element_text(size=14, vjust=-0.35, face="bold"),     
  #           axis.title.y = element_text(size=14, vjust=0.5)
  #     )
  plotAnova = bp
}

pdf(file="ME_Dx_Correlation.pdf")
par(mfrow = c(2,2))
for (i in 1:23){
  plotAnova(lm(MEs[-which(datMeta$dx == "HE" | datMeta$dx == "TR"),i] ~ dxWGCNA2),colnames(MEs)[i])
}
dev.off()

#### ME and Cell Type correlation?

pdf(file="ME v. Cell Type.pdf")
par(mfrow = c(2,2))
for (i in 1:23){
  for (ctype in rownames(est.prop.3.expr)){
  plotAnova(lm(MEs[-which(datMeta$dx == "HE" | datMeta$dx == "TR"),i] ~ est.prop.3.expr[ctype,]),colnames(MEs)[i],ctype)
  }
}
dev.off()

############################################ END OF HOLD OFF FOR NOW (may not make it to final)

### Match modules to significantly differentially expressed genes

tt_KO_sig = subset(tt_KO,tt_KO$adj.P.Val<0.05)
tt_KO_sig$module = colors[match(rownames(tt_KO_sig),rownames(datExpr))]
write.csv(tt_KO_sig,file="tt_KO_sig.csv")
table(tt_KO_sig$module)
tt_KO_sig=read.csv(file = "tt_KO_sig.csv",row.names=1)

tt_TR_sig = subset(tt_TR,tt_TR$adj.P.Val<0.05)
## only one gene - Arginase 1. and Vnn1 almost significantly different (in midnight blue)

tt_HE_sig = subset(tt_HE,tt_HE$adj.P.Val<0.05)
## only one gene - Arginase 1.

### Try PSI instead of CellMix for cell type enrichment in modules

#source('http://www.bioconductor.org/biocLite.R')
#biocLite("pSI")
library(pSI)
zhang.datExpr = read.csv("C:/Users/jillh/Dropbox/DHGLab/ArgMouse/proc/zhang_datExpr.cell.avg.csv")
rownames(zhang.datExpr) = zhang.datExpr$X
zhang.datExpr = zhang.datExpr[,-1]

pSI.output = specificity.index(pSI.in=zhang.datExpr,bts=50,p_max=.1, e_min=quantile(as.matrix(zhang.datExpr),0.05)); pSI.count(pSI.output)

cell.p.zhang = matrix(NA, length(unique(colors)),7);  rownames(cell.p.zhang) = unique(colors)#c("cyan","darkred","greenyellow","green","lightgreen","lightyellow","yellow","pink", "purple"); 
colnames(cell.p.zhang) = colnames(pSI.output)
cell.p.zhang.05 = cell.p.zhang
cell.p.zhang.01 = cell.p.zhang
cell.p.zhang.001 = cell.p.zhang

for(mod in rownames(cell.p.zhang)) {
  f = fisher.iteration(pSI.output, candidate.genes = mods[which(mods[,2]==mod),1],p.adjust = T)   ## mods created at beginning of analysis
  cell.p.zhang.05[mod,] = f[,1]  
  cell.p.zhang.01[mod,] = f[,2]  
  cell.p.zhang.001[mod,] = f[,3]  
}

library(gplots)

pdf("./ModuleAnnotation-CellType-PSI_0.001_all.pdf",width=6,height=5)
par(mar=c(4,2,2,0), oma=c(2,2,2,0))

cell.p.zhang = cell.p.zhang.001; main="pSI 0.001"
heatmap.2(-log10((cell.p.zhang)),col=blueWhiteRed(1000,1)[500:1000], 
          scale="none",trace="none",cexRow = 0.8,cexCol = .8, density.info = "none", 
          colsep=0:15,rowsep=0:15,sepcolor="grey",sepwidth=c(0.02,0.02),
          srtCol=45,offsetRow=0,offsetCol=-0.5, Colv = as.dendrogram(hclust(dist(t(zhang.datExpr)))),
          Rowv=as.dendrogram(hclust(dist(cell.p.zhang))),
          key=T,key.xlab="-log10(P)", cellnote=signif(cell.p.zhang,1), notecex=.8, notecol="black",main=main)


dev.off()

## CBPure dataset ((test2)
## new red = old greenyellow (Arg1)
## new pink = old yellow      (downregulated myelinating oligodendrocyte)
## new cyan = old cyan        (downregulated endothelial)

## CBPure dataset (test3 - removed low variance genes)
## new brown = old greenyellow (Arg1) - endothelial correlation
## new yellow = old cyan (downregulated endothelial)
## new midnightblue also correlated with endothelial (downregulated)
## new red = old yellow (downregulated myelinating oligodendrocyte)
## new salmon also correlated with m.o. (downregulated)
## pink and lightyellow - upregulated neuron (pink moreso than lightyellow)

## Angiogenesis genes: Bmp4|Casp8|Col8a1|Elk3|Ets1|Itgb3|Plxnd1|Pten|Shh|Tal1|Tbx4|Tie1
## ENSMUSG00000021835, ENSMUSG00000030123, ENSMUSG00000028717,   ENSMUSG00000026029 ,   ENSMUSG00000068196, ENSMUSG00000000094, ENSMUSG00000002633, ENSMUSG00000032035,   ENSMUSG00000008398,  ENSMUSG00000013663,   ENSMUSG00000020689,ENSMUSG00000033191    

angio = c("ENSMUSG00000021835", "ENSMUSG00000030123", "ENSMUSG00000028717",   "ENSMUSG00000026029" ,   "ENSMUSG00000068196", "ENSMUSG00000000094", "ENSMUSG00000002633", "ENSMUSG00000032035",   "ENSMUSG00000008398",  "ENSMUSG00000013663",   "ENSMUSG00000020689","ENSMUSG00000033191")
idx = which(mods[,1]%in%angio)
angio_mods = mods[idx,]

## yellow module definitely enriched for angiogenesis genes (midnightblue has them as well)

## Myelination genes: Amigo1|Aspa|Fgfr3|Hexb|Mbp|Mpdz|Myo5a|Myoc|Nkx6-2|Plp1|Tsc1|Ugt8a
##   ENSMUSG00000021665, ENSMUSG00000031425, ENSMUSG00000026812,   ENSMUSG00000020774, ENSMUSG00000041607,   ENSMUSG00000034593, ENSMUSG00000028402,   ENSMUSG00000054252,  ENSMUSG00000026697,  ENSMUSG00000041309, ENSMUSG00000032854, ENSMUSG00000050947   

myelin = c("ENSMUSG00000021665","ENSMUSG00000031425","ENSMUSG00000026812","ENSMUSG00000020774","ENSMUSG00000041607","ENSMUSG00000034593","ENSMUSG00000028402","ENSMUSG00000054252","ENSMUSG00000026697","ENSMUSG00000041309","ENSMUSG00000032854","ENSMUSG00000050947")
idx = which(mods[,1]%in%myelin)
myelin_mods = mods[idx,]

## salmon contains myelin genes (tan has one, grey has several)


angio_2 = c("ENSMUSG00000052336","ENSMUSG00000042903","ENSMUSG00000024986","ENSMUSG00000042286","ENSMUSG00000033191","ENSMUSG00000031250")
idx = which(rownames(datExpr)%in%angio_2)
angio_2dat = datExpr[idx,]
barplot(as.matrix(angio_2dat[6,]))

##### Update 05/02/18 - Assess some myelination genes of interest

#CGT        UDP-galactose:ceramide galactosyltransferase; search UGT8
#MOG      myelin oligoglycoprotein
#PLP         myelin proteolipid protein (lipophilin); search PLP1
#MBP       myelin basic protein
#MAG      myelin-associated glycoprotein
#CNPase  2',3'-cyclic nucleotide 3'-phosphodiesterase; search CNP

rm(list=ls())

setwd("C:/Users/jillh/Dropbox/DHGLab/ArgMouse/final_test3/")

ttable=read.csv("tt_KO_sig.csv")

grep("UGT8",ttable$external_gene_name)
grep("ENSMUSG00000032854",ttable$ensembl_gene_id)

grep("MOG",ttable$external_gene_name)
grep("ENSMUSG00000076439",ttable$ensembl_gene_id)

grep("PLP1",ttable$external_gene_name)
grep("ENSMUSG00000031425",ttable$ensembl_gene_id)

grep("MBP",ttable$external_gene_name)
grep("ENSMUSG00000041607",ttable$ensembl_gene_id)

grep("MAG",ttable$external_gene_name)
grep("ENSMUSG00000036634",ttable$ensembl_gene_id)

grep("CNP",ttable$external_gene_name)
grep("ENSMUSG00000006782",ttable$ensembl_gene_id)

### not in ttable, in datExpr?

datExpr=read.csv("datExpr_filt.csv",row.names=1)

grep("ENSMUSG00000032854",rownames(datExpr)) #8245
grep("ENSMUSG00000076439",rownames(datExpr)) #15876
grep("ENSMUSG00000031425",rownames(datExpr)) # 7507
grep("ENSMUSG00000041607",rownames(datExpr)) # 11007
grep("ENSMUSG00000036634",rownames(datExpr)) # 9348
grep("ENSMUSG00000006782",rownames(datExpr)) #971

### hm, maybe wrong data?

setwd("../")
load("test3/ttables_CBpure.RData")

dat_logFC=data.frame("Gene"=c("UGT8","MOG","PLP1","MBP","MAG","CNP"),
                     "logFC_KO"=c(tt_KO$logFC[which(tt_KO$ensembl_gene_id=="ENSMUSG00000032854")],
                       tt_KO$logFC[which(tt_KO$ensembl_gene_id=="ENSMUSG00000076439")],
                       tt_KO$logFC[which(tt_KO$ensembl_gene_id=="ENSMUSG00000031425")],
                       tt_KO$logFC[which(tt_KO$ensembl_gene_id=="ENSMUSG00000041607")],
                       tt_KO$logFC[which(tt_KO$ensembl_gene_id=="ENSMUSG00000036634")],
                       tt_KO$logFC[which(tt_KO$ensembl_gene_id=="ENSMUSG00000006782")]),
                     "p"=c(tt_KO$P.Value[which(tt_KO$ensembl_gene_id=="ENSMUSG00000032854")],
                           tt_KO$P.Value[which(tt_KO$ensembl_gene_id=="ENSMUSG00000076439")],
                           tt_KO$P.Value[which(tt_KO$ensembl_gene_id=="ENSMUSG00000031425")],
                           tt_KO$P.Value[which(tt_KO$ensembl_gene_id=="ENSMUSG00000041607")],
                           tt_KO$P.Value[which(tt_KO$ensembl_gene_id=="ENSMUSG00000036634")],
                           tt_KO$P.Value[which(tt_KO$ensembl_gene_id=="ENSMUSG00000006782")]),
                     "adj_p"=c(tt_KO$adj.P.Val[which(tt_KO$ensembl_gene_id=="ENSMUSG00000032854")],
                               tt_KO$adj.P.Val[which(tt_KO$ensembl_gene_id=="ENSMUSG00000076439")],
                               tt_KO$adj.P.Val[which(tt_KO$ensembl_gene_id=="ENSMUSG00000031425")],
                               tt_KO$adj.P.Val[which(tt_KO$ensembl_gene_id=="ENSMUSG00000041607")],
                               tt_KO$adj.P.Val[which(tt_KO$ensembl_gene_id=="ENSMUSG00000036634")],
                               tt_KO$adj.P.Val[which(tt_KO$ensembl_gene_id=="ENSMUSG00000006782")]))

ens_int_genes=c("ENSMUSG00000032854","ENSMUSG00000076439",
  "ENSMUSG00000031425","ENSMUSG00000041607","ENSMUSG00000036634","ENSMUSG00000006782") 

datExpr_intGenes=datExpr[match(ens_int_genes,rownames(datExpr)),]
rownames(datExpr_intGenes)=c("UGT8","MOG","PLP1","MBP","MAG","CNP")

dx=factor(c(rep("HET",3),rep("KO",3),rep("TR_KO",3),rep("WT",3)),levels=c("WT","TR_KO","HET","KO"))

mods=read.csv("final_test3/Genes_Modules.csv")

mods_int=mods[match(ens_int_genes,mods[,2]),]

### interestingly, all in the grey module

pdf("final_test3/MyelinGenes_boxplots.pdf")

for(gene in dat_logFC$Gene){
  idx=which(dat_logFC$Gene==gene)
  boxplot(as.numeric(datExpr_intGenes[idx,])~dx,col=c(4,2,3,"grey35"),
        ylab="Normalized Gene Expression",xlab="Condition",
        main=paste(dat_logFC$Gene[idx],": logFC_KO=",signif(dat_logFC$logFC_KO[idx],3),
                   ", p=",signif(dat_logFC$p[idx],3),", fdr_p_adj=",signif(dat_logFC$adj_p[idx],3),sep=""))
}
dev.off()
