#Load required R packages
library(DESeq2)
library(ggrepel)
library(limma)
#Set working directory
setwd("/share/crsp/lab/demetriou/rnaseq_Nov2018/")
#Read counts table
countsTable <- read.delim("all_merged.txt", header=T, stringsAsFactors=FALSE,check.names=F)
rownames(countsTable) <- countsTable$gene
countsTable <- countsTable[,-1]
#create metadata table
x=unlist(strsplit(colnames(countsTable)," "))
condition=factor(x[c(T,F)])
age=factor(substr(x[c(F,T)],1,1))
sex=factor(substr(x[c(F,T)],2,2))
colData=data.frame(age=age,condition=condition,sex=sex, condsSex=condition:sex)
#########################################
#Generate PCA using DESeq2
ddsFull<-DESeqDataSetFromMatrix(countsTable,colData,design=~1)
dds<-DESeq(ddsFull)
rld<-rlog(dds)
colnames(rld)=colnames(countsTable)
d=plotPCA(rld,intgroup=c("age","condition"),returnData=T)
percentVar=round(100*attr(d,"percentVar"))
ggplot(d, aes(PC1,PC2,color=age,shape=condition))+geom_point(size=4)+theme(text = element_text(size = 20))+
  xlab(paste0("PC1: ",percentVar[1],"% variance"))+ 
  theme(axis.text = element_text(size = 16))+ 
  ylab(paste0("PC2: ",percentVar[2],"% variance"))+geom_text_repel(aes(label=colnames(rld)),size=5)
write.csv(d[,1:2],"PCAcoordinates.csv",quote=T)
#########################################
#Generate MDS plot
rownames(colData)=colnames(countsTable)
plotMDS(assay(rld),top=500,pch=5,labels=colnames(rld),col=c(rep("purple",6),rep("blue",6),rep("green",6), rep("red",6)))
plotMDS(assay(rld),top=500,pch=5,col=c(rep("purple",6),rep("blue",6),rep("green",6), rep("red",6)))
##############
#Use DESeq2 for standard differential analysis
ct=countsTable[,condition=="Memory" & sex =="F"]
cD=colData[condition=="Memory" & sex =="F",]
ddsFull<-DESeqDataSetFromMatrix(data.matrix(ct),cD,design=~ age)
dds<-DESeq(ddsFull)
res<-results(dds)
res<-res[order(res$padj),]
x1=subset(res,padj<0.05)
write.csv(as.data.frame(x1),file="YoungvsOld_MemoryFemale.csv")
ct=countsTable[,condition=="Memory" & sex =="M"]
cD=colData[condition=="Memory" & sex =="M",]
ddsFull<-DESeqDataSetFromMatrix(data.matrix(ct),cD,design=~ age)
dds<-DESeq(ddsFull)
res<-results(dds)
res<-res[order(res$padj),]
x1=subset(res,padj<0.05)
write.csv(as.data.frame(x1),file="YoungvsOld_MemoryMale.csv")
ct=countsTable[,condition !="Memory" & sex =="F"]
cD=colData[condition!="Memory" & sex =="F",]
ddsFull<-DESeqDataSetFromMatrix(data.matrix(ct),cD,design=~ age)
dds<-DESeq(ddsFull)
res<-results(dds)
res<-res[order(res$padj),]
x1=subset(res,padj<0.05)
write.csv(as.data.frame(x1),file="YoungvsOld_NaïveFemale.csv")
ct=countsTable[,condition!="Memory" & sex =="M"]
cD=colData[condition!="Memory"& sex =="M",]
ddsFull<-DESeqDataSetFromMatrix(data.matrix(ct),cD,design=~ age)
dds<-DESeq(ddsFull)
res<-results(dds)
res<-res[order(res$padj),]
x1=subset(res,padj<0.05)
write.csv(as.data.frame(x1),file="YoungvsOld_NaïveMale.csv")
