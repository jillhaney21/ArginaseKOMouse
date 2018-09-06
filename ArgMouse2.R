## ArgMouse2
## Script to further analyze the ArgMouse dataset

options(stringsAsFactors = FALSE)
library(WGCNA)

## Plot ANOVA of different dx MEs

#i=5
#idx_anova = c(10,11,12,4,5,6,7,8,9,1,2,3)
mod=lm(MEs[,i] ~ datMeta$Dx)
nm=colnames(MEs)[i]
mod2 = MEs[,i] ~ datMeta$Dx
#a = aov(mod2)
#posthoc = TukeyHSD(x=a,'datMeta$Dx[idx_anova]',conf.level=0.95)

idx_anova = c(10,11,12,4,5,6,1,2,3,7,8,9)
datMeta = datMeta[idx_anova,]
datExpr = datExpr[,idx_anova]
MEs = MEs[idx_anova,]

plotAnova = function(mod,mod2,nm) {
  s = summary(mod)
  a = anova(mod)
  a2 = aov(mod2)
  posthoc = TukeyHSD(x=a2,'datMeta$Dx',conf.level=0.95)
  anova.pval = a$"Pr(>F)"[grep("Dx",rownames(a))]
  tukey.pval = posthoc$`datMeta$Dx`[1:3,4]
  title = paste(nm, ": ANOVA p =",signif(anova.pval,2),sep=" ")
  
  data=as.data.frame(s$coefficients[-1,-4])
  data = cbind(data,tukey.pval)
  data=data[grep("Dx",rownames(data)),]
  data$Dx = substring(rownames(data),11,15)
  colnames(data) = c("Beta", "SEM", "T", "p", "Group")
  col = rep("black",times = nrow(data))
  col[data[,4]<=0.1] = "pink"
  col[data[,4]<0.05] = "red"
  rownames(data) = data$Group
  
  ylim = c(min((data[,1]-4*data[,2])),max(data[,1]+4*data[,2]))
  
  bp = barplot(data[,1],names.arg = paste(rownames(data),paste("tukey p=",signif(data[,4],3),sep=""),sep="\n"), main=title, xlab="", ylab="beta",col = col,ylim=ylim,cex.names=0.75)
  segments(bp, data[,1]-data[,2], bp, data[,1]+data[,2], lwd=1)
  segments(bp-0.1, data[,1]+data[,2], bp+0.1, data[,1]+data[,2], lwd=1); segments(bp-0.1, data[,1]-data[,2], bp+0.1, data[,1]-data[,2], lwd=1)
  
  plotAnova = bp
}

datMeta$Dx = factor(datMeta$dx, levels= c("WT", "KO","HE","TR"))
pdf(file="ME_Dx_ANOVA.pdf")
par(mfrow = c(2,2))
for (i in 1:20){
  plotAnova(lm(MEs[,i] ~ datMeta$Dx),MEs[,i] ~ datMeta$Dx,colnames(MEs)[i])
}
dev.off()

## ORA: Are differentially expressed genes significantly enriched in modules?
## Run ORA.R first, with alternative="greater" within OR function fisher.test()

tt_KO_sig = read.csv("tt_KO_sig.csv",row.names=1)
ORAmat = matrix(NA,20,9)
modNames = unique(colors)
rownames(ORAmat) = modNames

for (i in 1:20){

module = modNames[i]
testpath = mods[which(mods[,2]==module),1] ## all genes in the module
refpath = tt_KO_sig$ensembl_gene_id  ## all differentially expressed genes
testbackground = rownames(datExpr)  
refbackground = rownames(datExpr)

oraout=ORA(testpath,refpath,testbackground,refbackground)
ORAmat[i,]=oraout

if(i==20){
  colnames(ORAmat)=names(oraout)
}

}

ORAmat = as.data.frame(ORAmat)
ORAmat$Fisher_log10 = -log10(as.numeric(as.character(ORAmat$`Fisher p`)))
ORAmat$p.adjust = p.adjust(ORAmat$`Fisher p`,method="fdr")
save(ORAmat,file="ORAmat_DiffExpr.RData")

### Plot module eigengenes of all four phenotypes for each module (boxplot)

dx = c("KO","TR","HE","WT")

pdf(file = "ME_Boxplot.pdf")
par(mfrow=c(2,2))
for (i in modNames){
print(i)
tempDat = MEs[,i]
tempDat = rbind(tempDat[4:6],tempDat[7:9],tempDat[1:3],tempDat[10:12])
tempDat = t(tempDat)
colnames(tempDat) = dx

boxplot(tempDat, range=0,col = c(2,3,4,5),names=dx, xlab = "Array", main = paste("Boxplot ME",i,sep=" "), ylab = "Intensity")
}
dev.off()

## Investigate ubiquitin pathway highlighted genes

setwd("C:/Users/Jill/Dropbox/DHGLab/ArgMouse/test3/")
load("temp_332016.RData")
load(file = "ttables_CBpure.RData")

ubiq = c('Bcl2','Hace1','Itch','Klhl42','March8','Rc3h2','Rffl','Rnf111','Rnf6','Shprh','Ube2q2','Clock','Crbn','Dtl','Myo5a','Rnf11','Sdf2l1','Spop','Tbl1xr1','Uba6','Uhrf2','Usp13','Usp21','Usp32','Usp45',' Arrb2','Kctd13','Nr1d1','Rbck1','Siah1a','Ube2s','Ube4b','Ufd1l','Vimp','Znrf1')

geneDat = read.csv("geneDat_CBPure_wMods.csv")

idx = match(geneDat$ensembl_gene_id,rownames(tt_KO))
tt_KO = tt_KO[idx,]

idx = which(geneDat$external_gene_name %in% ubiq)
tt_KO_ubiq = tt_KO[idx,]
tt_KO_ubiq = tt_KO_ubiq[order(tt_KO_ubiq$adj.P.Val),]

## Look at Arg1 in developing human brain

rm(list=ls())
options(stringsAsFactors = FALSE)
setwd("C:/Users/jillh/Dropbox/DHGLab/SCS_X")

expr = read.csv("genes_matrix_csv/expression_matrix.csv") ## RPKM of developing brain counts
datMeta = read.csv("genes_matrix_csv/columns_metadata.csv")
rownames = read.csv("genes_matrix_csv/rows_metadata.csv")

idx = which(rownames[,1] %in% expr[,1])
rownames = rownames[idx,]
idx = match(expr[,1],rownames[,1])
rownames(expr) = rownames$ensembl_gene_id[idx]
expr = expr[,-1]

regions = unique(datMeta$structure_acronym)
table(datMeta$structure_name)
ages_unq = unique(datMeta$age)

datMeta$age = factor(datMeta$age, levels= unique(datMeta$age))

setwd("../")
setwd("ArgMouse/test3/")

ages_bin = as.character(datMeta$age)
unique(ages_bin)
################################################# original
ages_bin = gsub("^8 pcw",1,ages_bin)
ages_bin = gsub("^9 pcw",1,ages_bin)
ages_bin = gsub("^12 pcw",2,ages_bin)
ages_bin = gsub("^13 pcw",3,ages_bin)
ages_bin = gsub("^16 pcw",4,ages_bin)
ages_bin = gsub("^17 pcw",5,ages_bin)
ages_bin = gsub("^19 pcw",5,ages_bin)
ages_bin = gsub("^21 pcw",5,ages_bin)
ages_bin = gsub("^24 pcw",5,ages_bin)
ages_bin = gsub("^25 pcw",6,ages_bin)
ages_bin = gsub("^26 pcw",6,ages_bin)
ages_bin = gsub("^35 pcw",6,ages_bin)
ages_bin = gsub("^37 pcw",6,ages_bin)
ages_bin = gsub("^4 mos",7,ages_bin)
ages_bin = gsub("^10 mos",8,ages_bin)
ages_bin = gsub("^1 yrs",9,ages_bin)
ages_bin = gsub("^2 yrs",10,ages_bin)
ages_bin = gsub("^3 yrs",11,ages_bin)
ages_bin = gsub("^4 yrs",11,ages_bin)
ages_bin = gsub("^8 yrs",12,ages_bin)
ages_bin = gsub("^11 yrs",12,ages_bin)
ages_bin = gsub("^13 yrs",12,ages_bin)
ages_bin = gsub("^15 yrs",13,ages_bin)
ages_bin = gsub("^18 yrs",13,ages_bin)
ages_bin = gsub("^19 yrs",13,ages_bin)
ages_bin = gsub("^21 yrs",14,ages_bin)
ages_bin = gsub("^23 yrs",14,ages_bin)
ages_bin = gsub("^30 yrs",14,ages_bin)
ages_bin = gsub("^36 yrs",14,ages_bin)
ages_bin = gsub("^37 yrs",14,ages_bin)
ages_bin = gsub("^40 yrs",14,ages_bin)
############################################################

############################################################ revised

ages_bin = gsub("^8 pcw",1,ages_bin)
ages_bin = gsub("^9 pcw",1,ages_bin)
ages_bin = gsub("^12 pcw",2,ages_bin)
ages_bin = gsub("^13 pcw",3,ages_bin)
ages_bin = gsub("^16 pcw",4,ages_bin)
ages_bin = gsub("^17 pcw",5,ages_bin)
ages_bin = gsub("^19 pcw",6,ages_bin)
ages_bin = gsub("^21 pcw",7,ages_bin)
ages_bin = gsub("^24 pcw",8,ages_bin)
ages_bin = gsub("^25 pcw",8,ages_bin)
ages_bin = gsub("^26 pcw",9,ages_bin)
ages_bin = gsub("^35 pcw",9,ages_bin)
ages_bin = gsub("^37 pcw",9,ages_bin)
ages_bin = gsub("^4 mos",10,ages_bin)
ages_bin = gsub("^10 mos",11,ages_bin)
ages_bin = gsub("^1 yrs",12,ages_bin)
ages_bin = gsub("^2 yrs",13,ages_bin)
ages_bin = gsub("^3 yrs",14,ages_bin)
ages_bin = gsub("^4 yrs",15,ages_bin)
ages_bin = gsub("^8 yrs",16,ages_bin)
ages_bin = gsub("^11 yrs",17,ages_bin)
ages_bin = gsub("^13 yrs",17,ages_bin)
ages_bin = gsub("^15 yrs",18,ages_bin)
ages_bin = gsub("^18 yrs",18,ages_bin)
ages_bin = gsub("^19 yrs",18,ages_bin)
ages_bin = gsub("^21 yrs",19,ages_bin)
ages_bin = gsub("^23 yrs",19,ages_bin)
ages_bin = gsub("^30 yrs",19,ages_bin)
ages_bin = gsub("^36 yrs",19,ages_bin)
ages_bin = gsub("^37 yrs",19,ages_bin)
ages_bin = gsub("^40 yrs",19,ages_bin)

############################################################ revised 2


ages_bin = gsub("^17 pcw",1,ages_bin)
ages_bin = gsub("^19 pcw",1,ages_bin)
ages_bin = gsub("^21 pcw",1,ages_bin)
ages_bin = gsub("^24 pcw",1,ages_bin)
ages_bin = gsub("^25 pcw",1,ages_bin)
ages_bin = gsub("^26 pcw",1,ages_bin)
ages_bin = gsub("^35 pcw",1,ages_bin)
ages_bin = gsub("^37 pcw",1,ages_bin)
ages_bin = gsub("^4 mos",1,ages_bin)
ages_bin = gsub("^10 mos",1,ages_bin)
ages_bin = gsub("^1 yrs",1,ages_bin)
ages_bin = gsub("^2 yrs",1,ages_bin)

ages_bin = gsub("^21 yrs",2,ages_bin)
ages_bin = gsub("^23 yrs",2,ages_bin)
ages_bin = gsub("^30 yrs",2,ages_bin)
ages_bin = gsub("^36 yrs",2,ages_bin)
ages_bin = gsub("^37 yrs",2,ages_bin)
ages_bin = gsub("^40 yrs",2,ages_bin)

idx1 = grep("pcw",ages_bin)
idx2 = grep("yrs",ages_bin)
ages_bin = ages_bin[-c(idx1,idx2)]

to_remove = c(idx1,idx2)

unique(ages_bin)

dataplot = cbind(log2(as.numeric(expr[4867,-to_remove])+1),ages_bin)
colnames(dataplot) = c("log2_Arg1_Normalized_Expression","Ages_Years")
dataplot = as.data.frame(dataplot)

#install.packages("scales")
library(scales)

pdf(file="Allen_Arg1_SmoothScatter_altplot.pdf")
par(mar=c(8,5,4,2))
scatter.smooth(dataplot$Ages_Years,as.numeric(dataplot$log2_Arg1_Normalized_Expression),span=0.66,
                pch=16,cex=2,col=alpha("gray50",0.25),lpars = list(lwd = 3.5),
                ylab="Log2 of Normalized Expression",xaxt='n',xlab="")
## org span = 0.2
## rev span = 0.33
#axis(1,at=seq(1:14),labels = c("8-9 pcw","12 pcw","13 pcw","16 pcw","17-24 pcw","25-37 pcw",
#             "4 mos","10 mos","1 yrs","2 yrs","3-4 yrs","8-13 yrs","15-19 yrs",
#              "21-40 yrs"),las=2)
axis(1,at=seq(1:19),labels = c("8-9 pcw","12 pcw","13 pcw","16 pcw","17 pcw","19 pcw","21 pcw",
                               "24-25 pcw","26-37 pcw","4 mos","10 mos","1 yrs","2 yrs","3 yr",
                               "4 yr","8 yr","11-13 yr","15-19 yr","21-40 yr"),las=2)
mtext("Age",side=1,line=5,cex.lab=30)
title("Arg1 Expression Across Human Development")
abline(v=9.5,col="red",lwd=2)
abline(h=0.1789978,col="blue",lwd=1.5)
dev.off()

## for scatter lowess fitting: "gaussian" fitting is by least-squares, and if family = "symmetric" a re-descending M estimator is used. 
## log scale isn't meaningful for this data

## Now make Allen region line plot

brain = data.frame(matrix(NA,nrow=26,ncol=3))
colnames(brain) = c("acronym","name","area")
brain[,2]=unique(datMeta$structure_name)
brain[,1]=unique(datMeta$structure_acronym)
brain[24:25,3]="cortical"
brain[26,3]="subcortical"## go through all 26 to create brain
save(brain,file="Brain_Regions_Allen.RData")
load("Brain_Regions_Allen.RData")
## add 'cortical' or 'subcortical' label to datMeta

datMeta.rev = datMeta[-to_remove,]
brain_region = datMeta.rev$structure_acronym
idx = which(brain[,3]=="cortical")
cortical = brain[idx,1]
idx = which(brain[,3]=="subcortical")
subcortical = brain[idx,1]

for(reg in cortical){
  brain_region = gsub(paste("\\b",reg,"\\b",sep=""),"cortical",brain_region)
}
for(reg in subcortical){
  brain_region = gsub(paste("\\b",reg,"\\b",sep=""),"subcortical",brain_region)
}

#idx.pre = which(ages_bin == 6)
#idx.post = which(ages_bin == 7 | ages_bin == 8 | ages_bin == 9 | ages_bin == 10)

idx.pre = which(ages_bin == 1)
idx.post = which(ages_bin == 2)

expr.rev = expr[,-to_remove]

dataplot.pre=cbind(log2(as.numeric(expr.rev[4867,idx.pre])+1),brain_region[idx.pre],
                   datMeta.rev$structure_acronym[idx.pre])
colnames(dataplot.pre) = c("log2_Arg1_Normalized_Expression","Region_Broad","Region_Specific")
dataplot.pre = as.data.frame(dataplot.pre)
dataplot.pre = dataplot.pre[order(dataplot.pre$Region_Broad),]
regs = unique(dataplot.pre[,3])
dataplot.pre.unique=as.data.frame(matrix(NA,nrow=17,ncol=4))
dataplot.pre.unique[,4]=regs
colnames(dataplot.pre.unique) = c("mean_log2_Arg1","se_log2_Arg1","Region_Broad","Region_Specific")
j=1
for(i in regs){
  idx = which(dataplot.pre[,3]==i)
  mn = mean(as.numeric(dataplot.pre[idx,1]))
  sem = sd(as.numeric(dataplot.pre[idx,1]))/sqrt(length(as.numeric(dataplot.pre[idx,1])))
  dataplot.pre.unique[j,1]=mn
  dataplot.pre.unique[j,2]=sem
  dataplot.pre.unique[j,3]=dataplot.pre[idx,2][1]
  j=j+1
}

dataplot.post=cbind(log2(as.numeric(expr.rev[4867,idx.post])+1),brain_region[idx.post],
                   datMeta.rev$structure_acronym[idx.post])
colnames(dataplot.post) = c("log2_Arg1_Normalized_Expression","Region_Broad","Region_Specific")
dataplot.post = as.data.frame(dataplot.post)
dataplot.post = dataplot.post[order(dataplot.post$Region_Broad),]
regs = unique(dataplot.post[,3])
dataplot.post.unique=as.data.frame(matrix(NA,nrow=16,ncol=4))
dataplot.post.unique[,4]=regs
colnames(dataplot.post.unique) = c("mean_log2_Arg1","se_log2_Arg1","Region_Broad","Region_Specific")
j=1
for(i in regs){
  idx = which(dataplot.post[,3]==i)
  mn = mean(as.numeric(dataplot.post[idx,1]))
  sem = sd(as.numeric(dataplot.post[idx,1]))/sqrt(length(as.numeric(dataplot.post[idx,1])))
  dataplot.post.unique[j,1]=mn
  dataplot.post.unique[j,2]=sem
  dataplot.post.unique[j,3]=dataplot.post[idx,2][1]
  j=j+1
}

idx = match(dataplot.post.unique[,4],dataplot.pre.unique[,4])
dataplot.pre.unique = dataplot.pre.unique[idx,]

brain2 = read.csv("brain_regs_allen.csv")[,-1]
xNames = brain2$figName[match(dataplot.post.unique[,4],brain2$acronym)]

pdf(file="Allen_Arg1_LinePlot_Age_Region.pdf")
par(mar=c(12,5,4,2))
plot(c(1:16),dataplot.pre.unique[,1],type="l",ylim=c(0,0.9),
     ylab="Mean Log2 of Normalized Expression",xaxt='n',xlab="",lwd=3,col="gray50")
arrows(c(1:16),dataplot.pre.unique[,1]-dataplot.pre.unique[,2],c(1:16),
       dataplot.pre.unique[,1]+dataplot.pre.unique[,2],length=0.04, angle=90, code=3,lwd=2,col="gray50")
points(c(1:16),dataplot.post.unique[,1],type="l",lwd=3,col="blue")
arrows(c(1:16),dataplot.post.unique[,1]-dataplot.post.unique[,2],c(1:16),
       dataplot.post.unique[,1]+dataplot.post.unique[,2],length=0.04, angle=90, code=3,col="blue",lwd=2)
axis(1,at=seq(1:16),labels = xNames,las=2,cex.axis=0.85)
mtext("Region",side=1,line=6,cex.lab=30,adj=0.85)
title("Arg1 Expression Across Human Development and Brain Region",cex.main=1)
legend("topright",inset=0.02,title="Age",c("17 pcw - 2 yrs","21 - 40 yrs"),fill=c("gray50","blue"))
dev.off()

## Now make a barplot with ggplot, grouped in twos, each region, both timepoints

library(ggplot2)

colnames(dataplot.pre.unique)[c(1,2)]=c("Mean.log2.Arg1","log2.Arg1.sem")
colnames(dataplot.post.unique)[c(1,2)]=c("Mean.log2.Arg1","log2.Arg1.sem")
dataplot = rbind(dataplot.pre.unique[,1:2],dataplot.post.unique[,1:2])
dataplot = cbind(dataplot,rbind(dataplot.pre.unique[,3:4],dataplot.post.unique[,3:4]))
dataplot = cbind(dataplot,c(rep("25 - 37 pcw",16),rep("4 mos - 2 yrs",16)))
colnames(dataplot)[c(4,5)] = c("Region","Age")

dataplot$Region <- factor(dataplot$Region,levels=dataplot$Region[1:16])

pdf(file="Allen_Arg1_BarPlot_Age_Region.pdf")
ggplot(dataplot,aes(x=Region, y=Mean.log2.Arg1,fill=Age)) +
                       geom_bar(stat="identity",position="dodge") + 
  geom_errorbar(aes(ymax=Mean.log2.Arg1+log2.Arg1.sem,ymin=Mean.log2.Arg1-log2.Arg1.sem),
                position=position_dodge(0.9),data=dataplot) +
  ggtitle("Mean Log2 of Normalized Arg1 Expression\nAcross Brain Regions and Ages")
dev.off()

pdf(file="Boxplot_BrainspanArg1.pdf")
par(mar=c(8,3,2,2))
boxplot(as.numeric(expr[4867,])~datMeta$age,range=0,main="Arg1 in HS by age",las=2,
        names=paste(names(table(datMeta$age)),table(datMeta$age),sep=",n="))
boxplot(as.numeric(expr[4867,])~datMeta$structure_acronym,range=0,main="Arg1 in HS by region",las=2,
        names=paste(names(table(datMeta$structure_acronym)),table(datMeta$structure_acronym),sep=",n="))
dev.off()

## now do a scatter instead of a boxplot

ages = unique(datMeta$age)

age_sum = rep(NA,31)
for(i in 1:31){
    age_sum[i] = mean(as.numeric(expr[4867,which(datMeta$age==ages[i])]))
}

pdf(file="Scatterplot_BrainspanArg1.pdf")
par(mar=c(5,5,4,2))
plot(as.numeric(datMeta$age),as.numeric(expr[4867,]),main="Arg1 in HS by age",ylab="Arg1 Expression",
     xaxt="n",xlab="Age")
axis(side=1,at=c(1:31),labels=ages,las=2,cex.axis=0.75)
lines(age_sum)
dev.off()

## make a region by age matrix

ageReg = matrix(NA,nrow=31,ncol=26)
rownames(ageReg)=ages
colnames(ageReg)=regions

for(i in 1:length(ages)){
  for (j in 1:length(regions)){
    idx = which(datMeta$age==ages[i] & datMeta$structure_acronym==regions[j])
    mn = median(as.matrix(expr[4867,idx]))
    ageReg[i,j]=mn
  }
}

#colors=colorRampPalette(c("red","yellow","green"))(n=600)
#heatmap(ageReg,na.rm=TRUE,Rowv = NA,Colv = NA,col=colors)

library(gplots)
library(WGCNA)

# following code limits the lowest and highest color to 5%, and 95% of your range, respectively
quantile.range <- quantile(ageReg, probs = seq(0, 1, 0.01),na.rm=TRUE)
palette.breaks <- seq(quantile.range["5%"], quantile.range["95%"], 0.1)

# use http://colorbrewer2.org/ to find optimal divergent color palette (or set own)
color.palette  <- colorRampPalette(c("green", "yellow", "red"))(length(palette.breaks) - 1)

pdf(file="HumanBrainMedianArg1Expr.pdf")

heatmap.2(ageReg,
  dendrogram = "none",
  scale      = "none",
  trace      = "none",
  key        = TRUE,
#  labRow     = NA,
#  labCol     = NA,
  col    = color.palette,
  breaks = palette.breaks,
  na.rm=TRUE,
  Rowv=FALSE,Colv=FALSE,
  cexRow = 0.8,cexCol = .8, density.info = "none", lhei = c(1,12),lwid=c(1,8),
  colsep=0:26,rowsep=0:31,sepcolor="grey",sepwidth=c(0.02,0.02),na.color="darkgrey",
  srtCol=45,offsetRow=0,offsetCol=-0.5,cellnote=signif(ageReg,1),notecol="black",notecex=0.5
#  main="Human Brain Median Arg1 Expression"
)

dev.off()

## no developing mouse brain data for Arg1

## from human brain data, seem to have increase in Arg1 expression in later pcw and early mos
## like between 24 pcw and 2 yrs old
## cerebellum, cerebellar cortex, visual cortex, dosolateral prefrontal cortex, and ventrolateral prefrontal cortex
## these areas visually seem to show more expression, especially in time periods identified above
## HOWEVER expression is definitely small - most normalized below 0.5, highest was ~ 1.5

## We see gene expression in adult mouse in isocortex+olfactory cortex most heavily, and throughout the brain (but still low throughout)
## ^ From Allen Brain Adult Mouse Brain Atlas 3D viewer

##### Look for module preservation in Allen developing mouse data
### Can't do because we would need the same strain of mouse PFC transcriptome - not publicly available

################## Look for enrichment of ID genes 

rm(list=ls())
options(stringsAsFactors = FALSE)
datMeta = read.csv("C:/Users/jillh/Dropbox/DHGLab/ArgMouse/proc/DatMeta.csv",row.names=1)
datMeta$Dx = factor(datMeta$dx, levels= c("WT", "KO", "HE", "TR"))
setwd("C:/Users/jillh/Dropbox/DHGLab/ArgMouse/test3/")
load("MEs_CBPure.RData")
datExpr = read.csv("datExpr_filt.csv",row.names=1)
load("ArgMsCols_CBPure_test3.RData")

id = read.csv("IDgenes_Nature2016.csv")

library(biomaRt)

ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL",dataset="mmusculus_gene_ensembl",host = "jul2015.archive.ensembl.org")

identifier <- "external_gene_name"
getinfo <- c("ensembl_gene_id", "entrezgene", "external_gene_name")
geneDat <- getBM(attributes = getinfo, filters=identifier, values = id[,1],mart=ensembl)
## 682 genes identified

id_ens = geneDat$ensembl_gene_id
table(is.na(id_ens)) # no NA

idx = match(id_ens,mods[,1])
mods_id = mods[idx,]
mods_id = mods_id[-which(is.na(mods_id[,1])),]
#605 genes total
table(mods_id[,2])

## now get list of all hugo gene id's that match ensembl gene ids

hugo = getBM(attributes=getinfo,mart=ensembl)
table(is.na(hugo[,3]))
hugo_genes = hugo[,3]
ens_idRef_genes = hugo[,1]


## ORA again but for ID genes: Are ID genes significantly enriched in modules?
## Run ORA.R first, with alternative="greater" within OR function fisher.test()

ORAmat_ID = matrix(NA,20,9)
modNames = unique(colors)
rownames(ORAmat_ID) = modNames

for (i in 1:20){
  
  module = modNames[i]
  testpath = mods[which(mods[,2]==module),1] ## all genes in the module
  refpath = id_ens  ## all id genes
  testbackground = rownames(datExpr)  
  refbackground = ens_idRef_genes
  
  oraout=ORA(testpath,refpath,testbackground,refbackground)
  ORAmat_ID[i,]=oraout
  
  if(i==20){
    colnames(ORAmat_ID)=names(oraout)
  }
  
}

ORAmat_ID = as.data.frame(ORAmat_ID)
ORAmat_ID$Fisher_log10 = -log10(as.numeric(as.character(ORAmat_ID$`Fisher p`)))
ORAmat_ID$p.adjust = p.adjust(ORAmat_ID[,2],method="fdr")
save(ORAmat_ID,file="ORAmat_ID.RData")
## tan is the only significantly enriched module for ID genes, it is significantly downregulated in KO

## for ANOVA, find adjusted p value

pANOVA = c(0.55,0.81,0.003,0.033,0.66,0.54,0.67,0.7,0.52,0.27,0.0069,0.0013,0.023,0.9,0.1,0.0045,0.018,0.084,0.015,0.13)
pANOVA.adjust = p.adjust(pANOVA,method="fdr")
pANOVA.adjust


## Look for blood brain barrier enriched genes in modules - genes expressed more in brain endothelial than lung and liver endothelial

bbb = read.csv(file="BBBgenes_Barres2010.csv",row.names=1)   ## blood brain barrier genes enriched in CNS endothelial

bbbprobes = rownames(bbb)

bbbgeneDat = getBM(attributes = c("affy_mouse430_2",getinfo), filters="affy_mouse430_2", values = bbbprobes,mart=ensembl)
which(is.na(bbbgeneDat$ensembl_gene_id))  ## no NA's

bbbgenes = unique(bbbgeneDat$ensembl_gene_id)
idx = which(rownames(datExpr) %in% bbbgenes)
mods_bbb = mods[idx,]
table(mods_bbb[,2])

bbbAllgenes = getBM(attributes = c("affy_mouse430_2",getinfo),mart=ensembl)
which(is.na(bbbAllgenes$ensembl_gene_id))
bbb_idRef_genes = unique(bbbAllgenes$ensembl_gene_id)

## Run ORA.R first, with alternative="greater" within OR function fisher.test()

ORAmat_BBB = matrix(NA,20,9)
modNames = unique(colors)
rownames(ORAmat_BBB) = modNames

for (i in 1:20){
  
  module = modNames[i]
  testpath = mods[which(mods[,2]==module),1] ## all genes in the module
  refpath = bbbgenes  ## all bbb genes
  testbackground = rownames(datExpr)  
  refbackground = bbb_idRef_genes
  
  oraout=ORA(testpath,refpath,testbackground,refbackground)
  ORAmat_BBB[i,]=oraout
  
  if(i==20){
    colnames(ORAmat_BBB)=names(oraout)
  }  
}

ORAmat_BBB = as.data.frame(ORAmat_BBB)
ORAmat_BBB$Fisher_log10 = -log10(as.numeric(as.character(ORAmat_BBB$`Fisher p`)))
ORAmat_BBB$p.adjust = p.adjust(ORAmat_BBB[,2],method="fdr")
save(ORAmat_BBB,file="ORAmat_BBB_CNS-specific.RData")
## brown, purple, pink, and black enriched - almost midnight blue

##### now look at bbbUp - CNS vascular genes upregulated in development

bbbup = read.csv(file="BBBupgenes_Barres2010.csv",row.names=1)   

bbbupprobes = rownames(bbbup)

bbbupgeneDat = getBM(attributes = c("affy_mouse430_2",getinfo), filters="affy_mouse430_2", values = bbbupprobes,mart=ensembl)
which(is.na(bbbupgeneDat$ensembl_gene_id))  ## no NA's

bbbUPgenes = unique(bbbupgeneDat$ensembl_gene_id)
idx = which(rownames(datExpr) %in% bbbUPgenes)
mods_bbbup = mods[idx,]
table(mods_bbbup[,2])

ORAmat_BBB_Up = matrix(NA,20,9)
modNames = unique(colors)
rownames(ORAmat_BBB_Up) = modNames

for (i in 1:20){
  
  module = modNames[i]
  testpath = mods[which(mods[,2]==module),1] ## all genes in the module
  refpath = bbbUPgenes  ## all bbb genes
  testbackground = rownames(datExpr)  
  refbackground = bbb_idRef_genes
  
  oraout=ORA(testpath,refpath,testbackground,refbackground)
  ORAmat_BBB_Up[i,]=oraout
  
  if(i==20){
    colnames(ORAmat_BBB_Up)=names(oraout)
  }  
}

ORAmat_BBB_Up = as.data.frame(ORAmat_BBB_Up)
ORAmat_BBB_Up$Fisher_log10 = -log10(as.numeric(as.character(ORAmat_BBB_Up$`Fisher p`)))
ORAmat_BBB_Up$p.adjust = p.adjust(ORAmat_BBB_Up[,2],method="fdr")
save(ORAmat_BBB_Up,file="ORAmat_BBB_CNS-Up.RData")
## purple, salmon, tan, brown, midnightblue, lightcyan, pink


##### now lood at bbbDown - CNS vascular genes downregulated in development

bbbdn = read.csv(file="BBBdowngenes_Barres2010.csv",row.names=1)   

bbbdnprobes = rownames(bbbdn)

bbbdngeneDat = getBM(attributes = c("affy_mouse430_2",getinfo), filters="affy_mouse430_2", values = bbbdnprobes,mart=ensembl)
which(is.na(bbbdngeneDat$ensembl_gene_id))  ## no NA's

bbbDNgenes = unique(bbbdngeneDat$ensembl_gene_id)
idx = which(rownames(datExpr) %in% bbbDNgenes)
mods_bbbdn = mods[idx,]
table(mods_bbbdn[,2])

ORAmat_BBB_Dn = matrix(NA,20,9)
modNames = unique(colors)
rownames(ORAmat_BBB_Dn) = modNames

for (i in 1:20){
  
  module = modNames[i]
  testpath = mods[which(mods[,2]==module),1] ## all genes in the module
  refpath = bbbDNgenes  ## all bbb genes
  testbackground = rownames(datExpr)  
  refbackground = bbb_idRef_genes
  
  oraout=ORA(testpath,refpath,testbackground,refbackground)
  ORAmat_BBB_Dn[i,]=oraout
  
  if(i==20){
    colnames(ORAmat_BBB_Dn)=names(oraout)
  }  
}

ORAmat_BBB_Dn = as.data.frame(ORAmat_BBB_Dn)
ORAmat_BBB_Dn$Fisher_log10 = -log10(as.numeric(as.character(ORAmat_BBB_Dn$`Fisher p`)))
ORAmat_BBB_Dn$p.adjust = p.adjust(ORAmat_BBB_Dn[,2],method="fdr")
save(ORAmat_BBB_Dn,file="ORAmat_BBB_CNS-down.RData")
## lightgreen, yellow, cyan, magenta, pink


################## now genes that are upregulated during bbb development AND CNS specific

bbbUp_CNS = intersect(bbbUPgenes,bbbgenes)
idx = which(rownames(datExpr)%in%bbbUp_CNS)
mods_bbbUp_CNS = mods[idx,]
table(mods_bbbUp_CNS[,2])

ORAmat_BBB_UpCN = matrix(NA,20,9)
modNames = unique(colors)
rownames(ORAmat_BBB_UpCN) = modNames

for (i in 1:20){
  
  module = modNames[i]
  testpath = mods[which(mods[,2]==module),1] ## all genes in the module
  refpath = bbbUp_CNS  ## all bbb genes
  testbackground = rownames(datExpr)  
  refbackground = bbb_idRef_genes
  
  oraout=ORA(testpath,refpath,testbackground,refbackground)
  ORAmat_BBB_UpCN[i,]=oraout
  
  if(i==20){
    colnames(ORAmat_BBB_UpCN)=names(oraout)
  }  
}

ORAmat_BBB_UpCN = as.data.frame(ORAmat_BBB_UpCN)
ORAmat_BBB_UpCN$Fisher_log10 = -log10(as.numeric(as.character(ORAmat_BBB_UpCN$`Fisher p`)))
ORAmat_BBB_UpCN$p.adjust = p.adjust(ORAmat_BBB_UpCN[,2],method="fdr")
save(ORAmat_BBB_UpCN,file="ORAmat_BBB_Up_CNS-specific.RData")
###brown the only module enriched for genes that are upregulated during bbb development AND CNS specific

################### now gene that are downregulated during bbb development and CNS specific

bbbDn_CNS = intersect(bbbDNgenes,bbbgenes)
idx = which(rownames(datExpr)%in%bbbDn_CNS)
mods_bbbDn_CNS = mods[idx,]
table(mods_bbbDn_CNS[,2])

ORAmat_BBB_DnCN = matrix(NA,20,9)
modNames = unique(colors)
rownames(ORAmat_BBB_DnCN) = modNames

for (i in 1:20){
  
  module = modNames[i]
  testpath = mods[which(mods[,2]==module),1] ## all genes in the module
  refpath = bbbDn_CNS  ## all bbb genes
  testbackground = rownames(datExpr)  
  refbackground = bbb_idRef_genes
  
  oraout=ORA(testpath,refpath,testbackground,refbackground)
  ORAmat_BBB_DnCN[i,]=oraout
  
  if(i==20){
    colnames(ORAmat_BBB_DnCN)=names(oraout)
  }  
}

ORAmat_BBB_DnCN = as.data.frame(ORAmat_BBB_DnCN)
ORAmat_BBB_DnCN$Fisher_log10 = -log10(as.numeric(as.character(ORAmat_BBB_DnCN$`Fisher p`)))
ORAmat_BBB_DnCN$p.adjust = p.adjust(ORAmat_BBB_DnCN[,2],method="fdr")
## none significant, not suprising

####### Check out the bbbDiff_Endo genes which are differentially expressed in our group

tt_KO_sig = read.csv("tt_KO_sig.csv",row.names=1)
idx = which( rownames(tt_KO_sig) %in% bbbUp_CNS)
tt_KO_bbbUp_CNS = tt_KO_sig[idx,]

idx = which( rownames(tt_KO_sig) %in% bbbDn_CNS)
tt_KO_bbbDn_CNS = tt_KO_sig[idx,]

idx = which( rownames(tt_KO_sig) %in% bbbUPgenes)
tt_KO_bbbUp = tt_KO_sig[idx,]

idx = which( rownames(tt_KO_sig) %in% bbbDNgenes)
tt_KO_bbbDN = tt_KO_sig[idx,]

idx = which( rownames(tt_KO_sig) %in% bbbgenes)
tt_KO_bbbCNS = tt_KO_sig[idx,]

## identify the genes in the brown, midnightblue, salmon, and magenta modules that are bbb enriched

# First, the CNS bbb specific genes in brown and midnightblue

idx = which(mods[,2] == "brown")
brown = mods[idx,]
idx = which(mods[,2] == "midnightblue")
midnightblue = mods[idx,]
idx = which(mods[,2] == "salmon")
salmon = mods[idx,]
idx = which(mods[,2] == "magenta")
magenta = mods[idx,]

idx = which(brown[,1] %in% bbbgenes)
brown_bbbCNS = brown[idx,1]
idx = which(midnightblue[,1] %in% bbbgenes)
midnightblue_bbbCNS = midnightblue[idx,1]

# Now the 'supposed to be' upregulated bbb genes in midnightblue, salmon, and brown

idx = which(brown[,1] %in% bbbUPgenes)
brown_bbbUp = brown[idx,1]
idx = which(midnightblue[,1] %in% bbbUPgenes)
midnightblue_bbbUp = midnightblue[idx,1]
idx = which(salmon[,1] %in% bbbUPgenes)
salmon_bbbUp = salmon[idx,1]

# Now the 'supposed to be' downregulated bbb genes in magenta

idx = which(magenta[,1] %in% bbbDNgenes)
magenta_bbbDn = magenta[idx,1]

# Finally, the bbb genes that are CNS specific and upregulated

idx = which(brown[,1] %in% bbbUp_CNS)
brown_bbbUp_CNS = brown[idx,1]

outersect <- function(x, y) {
  sort(c(setdiff(x, y),
         setdiff(y, x)))
}

brown_CNSonly = outersect(brown_bbbCNS,brown_bbbUp_CNS)
brown_bbbUponly = outersect(brown_bbbUp,brown_bbbUp_CNS)

length(brown_bbbUp_CNS) <- length(brown_CNSonly) <- length(brown_bbbUponly) <- 54

all_brown = cbind(brown_bbbUp_CNS,brown_CNSonly,brown_bbbUponly)
colnames(all_brown) = c("Up_CNS","CNS_only","Up_only")


midnightblue_bbbUp_CNS = intersect(midnightblue_bbbCNS,midnightblue_bbbUp)
length(midnightblue_bbbUp_CNS) <- length(midnightblue_bbbCNS) <- length(midnightblue_bbbUp) <- 26

all_midnightblue = cbind(midnightblue_bbbUp_CNS,midnightblue_bbbCNS,midnightblue_bbbUp)
colnames(all_midnightblue) = c("Up_CNS","CNS_only","Up_only")

save(ORAmat_BBB,ORAmat_BBB_Dn,ORAmat_BBB_Up,ORAmat_BBB_UpCN,all_brown,all_midnightblue,salmon_bbbUp,magenta_bbbDn,tt_KO_bbbCNS,tt_KO_bbbUp,
     tt_KO_bbbDN,tt_KO_bbbDn_CNS,tt_KO_bbbUp_CNS,file="BBB_enrichment.RData")
## also created BBB_enrichmend.xlxs with the four differentially expressed modules and their enriched bbb genes 

write.csv(tt_KO_bbbCNS,file="tt_KO_BBB_CNS.csv")
write.csv(tt_KO_bbbDN,file="tt_KO_BBB_Down.csv")
write.csv(tt_KO_bbbUp,file="tt_KO_BBB_Up.csv")
write.csv(tt_KO_bbbDn_CNS,file = "tt_KO_BBB_DownCNS.csv")
write.csv(tt_KO_bbbUp_CNS,file="tt_KO_BBB_UpCNS.csv")

## make a boxplot of Arg1 expression across phenotypes

arg1 = which(rownames(datExpr)=="ENSMUSG00000019987")

arg1_KO = mean(as.matrix(datExpr[arg1,which(datMeta$dx == "KO")]))
arg1_KO_sem = sem<-sd(as.matrix(datExpr[arg1,which(datMeta$dx == "KO")]))/sqrt(length(as.matrix(datExpr[arg1,which(datMeta$dx == "KO")])))
arg1_WT = mean(as.matrix(datExpr[arg1,which(datMeta$dx == "WT")]))
arg1_WT_sem = sem<-sd(as.matrix(datExpr[arg1,which(datMeta$dx == "WT")]))/sqrt(length(as.matrix(datExpr[arg1,which(datMeta$dx == "WT")])))
arg1_TR = mean(as.matrix(datExpr[arg1,which(datMeta$dx == "TR")]))
arg1_TR_sem = sem<-sd(as.matrix(datExpr[arg1,which(datMeta$dx == "TR")]))/sqrt(length(as.matrix(datExpr[arg1,which(datMeta$dx == "TR")])))
arg1_HE = mean(as.matrix(datExpr[arg1,which(datMeta$dx == "HE")]))
arg1_HE_sem = sem<-sd(as.matrix(datExpr[arg1,which(datMeta$dx == "HE")]))/sqrt(length(as.matrix(datExpr[arg1,which(datMeta$dx == "HE")])))

library(ggplot2)
dataplot = as.data.frame(matrix(NA,nrow=4,ncol=3))
dataplot[,1] = unique(datMeta$dx)
dataplot[,2] = c(arg1_WT,arg1_KO,arg1_HE,arg1_TR)
dataplot[,3] = c(arg1_WT_sem,arg1_KO_sem,arg1_HE_sem,arg1_TR_sem)
colnames(dataplot)=c("Phenotype","Normalized_Gene_Expression","SEM")
dataplot = within(dataplot,Phenotype <- factor(Phenotype,levels=c("WT","KO","TR","HE")))

library(grid)
pdf(file="Arg1_Mean_Expression.pdf")
ggplot(dataplot,aes(x=Phenotype,y=Normalized_Gene_Expression))+
  geom_bar(position="dodge",stat="identity",fill = c("blue3","gold2","magenta3","green3"))+
  geom_errorbar(aes(ymin=Normalized_Gene_Expression-SEM,ymax=Normalized_Gene_Expression+SEM),width=0.25)+
  labs(title="Mean Arginase 1 Gene Expression",y = "Mean of Normalized Expression",
       x="Genotype") +
  theme(plot.title = element_text(face="bold",margin = margin(t=10,b=25)),
        axis.title.y = element_text(margin = margin(0,15,0,0)),
        axis.title.x = element_text(margin = margin(t=15,b=15)))
dev.off()

## human orthologs of all modules for Dapple PPI analysis

rm(list=ls())
setwd("C:/Users/Jill/Dropbox/DHGLab/ArgMouse/final/DAPPLE/")

brown = read.csv("brown.csv",header=FALSE)
magenta = read.csv("magenta.csv",header=FALSE)
midnightblue = read.csv("midnightblue.csv",header=FALSE)
salmon = read.csv("salmon.csv",header=FALSE)

mods = data.frame(matrix(NA,1186,2))
mods[,1] = c(brown[,1],magenta[,1],midnightblue[,1],salmon[,1])
mods[,2] = c(rep("brown",300),rep("magenta",300),rep("midnightblue",286),rep("salmon",300))

library(biomaRt)

mouse = useMart(biomart="ENSEMBL_MART_ENSEMBL",dataset="mmusculus_gene_ensembl",host = "jul2015.archive.ensembl.org")
a = listAttributes(mouse)

attributes =c("ensembl_gene_id","hsapiens_homolog_ensembl_gene")
orth.human=getBM(attributes,filters="with_homolog_hsap",values=TRUE,mart=mouse,bmHeader=FALSE)
dim(orth.human)

idx = match(mods[,1],orth.human[,1])
orth.human=orth.human[idx,]

table(is.na(orth.human[,2]))
mods_human = mods[!is.na(orth.human[,2]),]
orth.human = orth.human[!is.na(orth.human[,2]),]
mods_human[,1] = orth.human[,2]

###
human = useMart(biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl",host = "jul2015.archive.ensembl.org")
attributes = c("ensembl_gene_id","hgnc_symbol")
human.hgnc=getBM(attributes,filters="ensembl_gene_id",values=mods_human[,1],mart=human)

idx = match(mods_human[,1],human.hgnc$ensembl_gene_id)
human.hgnc = human.hgnc[idx,]

mods_human = mods_human[!is.na(human.hgnc$hgnc_symbol),]
human.hgnc = human.hgnc[!is.na(human.hgnc$hgnc_symbol),]

mods_human[,1] = human.hgnc[,2]
idx = which(mods_human == "")
mods_human = mods_human[-idx,]
idx = duplicated(mods_human[,1])
mods_human = mods_human[!idx,]

# Now human module lists

modules = list()
modNames=unique(mods_human[,2])
background=mods_human[,1]
for(i in 1:4){
  mod = subset(background,mods_human[,2]==modNames[i])
  mod_an = gsub("[^[:alnum:] ]", "REMOVE", mod)
  idx = grep("REMOVE",mod_an)
  mod_an = mod_an[-idx]
  modules[[i]] = mod_an
  names(modules)[i] = paste("MM",modNames[i],sep="")
}

save(mods_human,background,file="ArgMouse_Analysis_CBPure_humanOrthologs_hgnc300kIN.RData")
#load("ArgMouse_Analysis_CBPure.RData")

## Now csv's for Dapple

for (i in 1:4){
  write.table(modules[[i]],file=paste("./Modules_Human_HGNC_300kIN/",names(modules)[i],".txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
}

#### Make a cell type smooth scatter plot of Arg1, human, different ages

options(stringsAsFactors = FALSE)
setwd("C:/Users/jillh/Dropbox/DHGLab/ArgMouse/GSE67835_HS_Cells/")

## Get datMeta

datMeta = t(read.csv("datMeta_raw.csv",header=FALSE))
colnames(datMeta) = datMeta[1,]
rownames(datMeta) = datMeta[,2]
datMeta = datMeta[-1,]
datMeta = as.data.frame(datMeta)

idx = grep("years",datMeta$Age)
datMeta$Pre_Post_Age[idx] = "Post"
idx = grep("W",datMeta$Age)
datMeta$Pre_Post_Age[idx] = "Pre"
datMeta$Age=gsub(" years","",datMeta$Age)
datMeta$Age=gsub(" W","",datMeta$Age)

write.csv(datMeta,file="datMeta_proc.csv")

## Get counts

raw = list()

samples = list.files("./GSE67835_RAW/")  ## A survey of human brain transcriptome diversity at the single cell level., ye zhang and barres 2015
idx = grep(".gz",samples)
samples = samples[-idx]

for(samp in samples){
  name = substr(samp,1,10)
  raw[[name]] = read.csv(paste("./GSE67835_RAW/",samp,"/",samp,sep=""),sep="\t")
}

arg1 = data.frame(matrix(NA,nrow=466,ncol=2))
arg1[,1] = substr(samples,1,10)

for(samp in arg1[,1]){
  temp = raw[[samp]][grep("ARG1",raw[[samp]][,1]),]
  arg1[grep(samp,arg1[,1]),2] = temp[,2]
}

save(arg1,file="Arg1_Counts.RData")
load("Arg1_Counts.RData")

## Start compiling scatter plot data

unique(datMeta$Cell_Type) ## don't include hybrid type

to_remove = grep("hybrid",datMeta$Cell_Type)
datMeta = datMeta[-to_remove,]
arg1 = arg1[-to_remove,]

arg1[,2] = log2(arg1[,2]+1)
quantile(arg1[,2],c(0,.25,.5,.75,.8,.85,.9,.95,1))

unique(datMeta$Age) ## Going to look at pre and post age group averages

## Smooth scatter - black for preconception weeks, blue for adult

# remove this sample, probably incorrectly labeled cell type

to_remove = which(rownames(datMeta)=="GSM1658003")
datMeta = datMeta[-to_remove,]
arg1 = arg1[-to_remove,]

datMeta$Cell_Type = factor(datMeta$Cell_Type,levels = c("fetal_replicating","fetal_quiescent",
                              "astrocytes","endothelial","microglia","neurons","OPC","oligodendrocytes"))
datMeta$CellBin = as.numeric(datMeta$Cell_Type)
idxPreC = which(datMeta$Pre_Post_Age == "Pre" & datMeta$Tissue == "cortex")
idxPreC = idxPreC[-which(as.numeric(arg1$X2)[idxPreC]==0)]
idxPostC = which(datMeta$Pre_Post_Age == "Post" & datMeta$Tissue == "cortex")
idxPostC = idxPostC[-which(as.numeric(arg1$X2)[idxPostC]==0)]
#idxPreH = which(datMeta$Pre_Post_Age == "Pre" & datMeta$Tissue == "hippocampus") #none, expected
idxPostH = which(datMeta$Pre_Post_Age == "Post" & datMeta$Tissue == "hippocampus")
idxPostH = idxPostH[-which(as.numeric(arg1$X2)[idxPostH]==0)]

mns = vector()
sem = vector()
for(i in 1:8){
x=as.numeric(arg1$X2)[which(datMeta$CellBin==i)]
mns = c(mns,mean(x))
sem = c(sem,sd(x)/sqrt(length(x)))
}

library(scales)

pdf(file="Barres_CellType_Arg1_SmoothScatter.pdf")
par(mar=c(9,5,4,9),xpd=TRUE)
bar <- barplot(mns,col=c("gray","gray",rep("blue",6)),ylim=c(0,7),ylab="Mean of Log2 of Expression")

xVal <- as.character(datMeta$CellBin)
j=1
for (i in c("1","2","3","4","5","6","7","8")){
  xVal <- gsub(paste("^",i,"$",sep=""),bar[j],xVal)
  j=j+1
}

points(xVal[idxPreC],as.numeric(arg1$X2)[idxPreC],col=alpha("gray50",0.25),
     pch=16,cex=2,xaxt='n',xlab="",ylab="Log2.Arg1")
points(xVal[idxPostC],as.numeric(arg1$X2)[idxPostC],col=alpha("blue",0.25),
     pch=16,cex=2)
points(xVal[idxPostH],as.numeric(arg1$X2)[idxPostH],col=alpha("blue",0.25),
       pch=17,cex=2)
arrows(bar,mns-sem,bar,mns+sem,length=0.04, angle=90,code=3,col="black",lwd=1.5)
axis(1,at=bar,labels = levels(as.factor(datMeta$Cell_Type)),las=2)
mtext("Cell.Type",side=1,line=7,cex.lab=30)
title("Human Arg1 Expression Across\nAge, Cell Type, and Brain Region")
legend(9,6.5,c("Cortex\n16-18 pcw","Cortex\n21-63 yrs","Hippocampus\n21-63 yrs"),
       pch=c(16,16,17),cex=1.1,col=c("gray50","blue","blue"),y.intersp = 1.5)
dev.off()

###################### Plot salmon module

setwd("C:/Users/jillh/Dropbox/DHGLab/ArgMouse/test3")

datMeta = read.csv("C:/Users/jillh/Dropbox/DHGLab/ArgMouse/proc/DatMeta.csv",row.names=1)
datMeta$Dx = factor(datMeta$dx, levels= c("WT", "KO", "HE", "TR"))
datExpr = read.csv("datExpr_filt.csv",row.names=1)

load("MEs_CBPure.RData")
load("ModCharacteristics_CBPure.RData")
load("ArgMouse_Analysis_CBPure.RData")
load("ArgMsCols_CBPure_test3.RData")

geneDat = read.csv("geneDat_CBPure_wMods.csv",row.names=1)
softThresh = 18

for(mod in moduleColors) {
  mod="salmon"
  module_genes = rownames(datExpr)[labels2colors(merged$colors) == mod]
  gene_idx = order(kME[module_genes,paste("kME", mod, sep="")],decreasing = T)
  
  module_kme = kME[module_genes, paste("kME", mod, sep="")]
  
  range(module_kme[gene_idx[1:100]])
  
  #numgenesingraph = 29;
  numconnections2keep = 500;
  
  genesinMod = module_genes[gene_idx[1:numgenesingraph]]
  genes = geneDat$external_gene_name[match(module_genes[gene_idx],geneDat$ensembl_gene_id)]
  genes = genes[1:numgenesingraph]
  modTOM = TOMsimilarityFromExpr(t(datExpr[genesinMod,]), power = softThresh);
  dimnames(modTOM) = list(genes, genes)
  
  reducedTOM = modTOM
  orderind = order(reducedTOM,decreasing=TRUE);
  connections2keep = orderind[(1+nrow(modTOM)):(numconnections2keep+nrow(modTOM))];
  reducedTOM = matrix(0,nrow(reducedTOM),ncol(reducedTOM));
  reducedTOM[connections2keep] = 1;
  
  g0 <- graph.adjacency(as.matrix(reducedTOM[1:3,1:3]),mode="undirected",weighted=TRUE,diag=FALSE)
  layoutMata <- layout.circle(g0)
  
  g0 <- graph.adjacency(as.matrix(reducedTOM[4:29,4:29]),mode="undirected",weighted=TRUE,diag=FALSE)
  layoutMatb <- layout.circle(g0)
  
  #g0 <- graph.adjacency(as.matrix(reducedTOM[30:ncol(reducedTOM),30:ncol(reducedTOM)]),mode="undirected",weighted=TRUE,diag=FALSE)
  #layoutMatc <- layout.circle(g0)
  
  g1 <- graph.adjacency(as.matrix(reducedTOM),mode="undirected",weighted=TRUE,diag=FALSE)
  layoutMat <- rbind(layoutMata*0.15,layoutMatb*0.5)
  
  CairoFonts(regular="Arial:style=Regular",bold="Arial:style=Bold",italic="Arial:style=Italic",bolditalic="Arial:style=Bold Italic,BoldItalic",symbol="Symbol")
  plot.igraph(g1,edge.color="grey",vertex.color=mod,vertex.label=as.character(genes),
              vertex.label.cex=0.75,vertex.label.color="black",
              layout= layoutMat,vertex.size=20,edge.width = .5*E(g1)$weight^10,
              main = "Salmon Module")
  
  #   plot(g1,edge.color="grey",vertex.color=mod,vertex.label=as.character(genes),
  #        vertex.label.cex=0.7,vertex.label.dist=0.45,vertex.label.degree=-pi/4,
  #        vertex.label.color="black",layout= layoutMat,
  #        vertex.size=module[1:numgenesingraph]^3*4,
  #        main=paste(mod,"module"))
  
  # plot.igraph(g1, vertex.label = geneSymbols,
  #             vertex.label.dist=0.3, 
  #             vertex.size=4,
  #             vertex.label.color="black",
  #             vertex.color = labels2colors(merged$colors)[gene_idx],
  #             vertex.label.cex=0.4,
  #             vertex.frame.color="black",
  #             layout=layoutFR,
  #             edge.color="green")# 
  
}













