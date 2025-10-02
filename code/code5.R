library(stringr)

setwd("D:\\data\\4_Transcriptionfactor&logics")
mrna <- read.table("mRNA.txt",header=T,row.names = 1,check.names = F) 
lncrna <- read.table("lncRNA.txt",header=T,row.names = 1,check.names = F)  
mrnat <- t(mrna)  
lncrnat <- t(lncrna)

cancer_samples_mrna <- subset(mrnat, str_sub(rownames(mrnat),start=14,end=15) < 10) 
cancer_samples_lncrna <- subset(lncrnat, str_sub(rownames(lncrnat),start=14,end=15) < 10)  
rowdata_mrna <- str_sub(rownames(cancer_samples_mrna),end=12)  
rowdata_lncrna <- str_sub(rownames(cancer_samples_lncrna),end=12)
row.names(cancer_samples_mrna) <- rowdata_mrna
row.names(cancer_samples_lncrna) <- rowdata_lncrna
write.csv(cancer_samples_mrna,"mrna.csv",row.names = T)
write.csv(cancer_samples_lncrna,"lncrna.csv",row.names = T)
mrna1 <- read.csv("mrna.csv",header = T)
lncrna1 <- read.csv("lncrna.csv",header = T)

dup_rows_mrna <- duplicated(mrna1[, 1])
dup_rows_lncrna <- duplicated(lncrna1[, 1])
mrna2 <- mrna1[!dup_rows_mrna, ]
lncrna2 <- mrna1[!dup_rows_lncrna, ]
rownames(mrna2) <- mrna2[,1]
rownames(lncrna2) <- lncrna2[,1]
mrna2 <- mrna2[-1]
lncrna2 <- lncrna2[-1]

cluster <- read.table("ICIcluster2.txt", header=F, sep="\t", check.names=F, row.names=1)
colnames(cluster)[1] <- "subtype"

sameSample=intersect(row.names(cluster), row.names(mrna2))
mrna3=mrna2[sameSample,]
lncrna3=lncrna2[sameSample,]
write.csv(mrna3,"mrna.csv",row.names = T)
write.csv(lncrna3,"lncrna.csv",row.names = T)

mrna4=cbind(mrna2[sameSample,], subtype=cluster[sameSample,])
mrna5 <- mrna4[order(mrna4$subtype),]
mrna5 <- mrna5[, -ncol(mrna5)]

lncrna4=cbind(lncrna2[sameSample,], subtype=cluster[sameSample,])
lncrna5 <- lncrna4[order(lncrna4$subtype),]
lnrna5 <- lncrna5[, -ncol(lncrna5)]

normal=114 
tumor=50   

library("edgeR")
foldChange=1
padj=0.001
rt=as.matrix(t(lncrna5))
exp=rt[,1:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>1,]

group=c(rep("subtype1",normal),rep("subtype2",tumor))             
design <- model.matrix(~group)
y <- DGEList(counts=data,group=group)
y <- calcNormFactors(y)
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)
et <- exactTest(y,pair = c("subtype1","subtype2"))
topTags(et)
ordered_tags <- topTags(et, n=100000)

allDiff=ordered_tags$table
allDiff=allDiff[is.na(allDiff$FDR)==FALSE,]
diff=allDiff
newData=y$pseudo.counts

write.table(diff,file="edgerOut.xls",sep="\t",quote=F)
diffSig = diff[(diff$FDR < padj & (diff$logFC>foldChange | diff$logFC<(-foldChange))),]
write.table(diffSig, file="diffSig.xls",sep="\t",quote=F)
diffUp = diff[(diff$FDR < padj & (diff$logFC>foldChange)),]
write.table(diffUp, file="up.xls",sep="\t",quote=F)
diffDown = diff[(diff$FDR < padj & (diff$logFC<(-foldChange))),]
write.table(diffDown, file="down.xls",sep="\t",quote=F)

normalizeExp=rbind(id=colnames(newData),newData)
write.table(normalizeExp,file="normalizeExp.txt",sep="\t",quote=F,col.names=F)  
diffExp=rbind(id=colnames(newData),newData[rownames(diffSig),])
write.table(diffExp,file="diffmRNAExp.txt",sep="\t",quote=F,col.names=F)        

library("edgeR")
foldChange=1
padj=0.001
rt=as.matrix(t(mrna5))
exp=rt[,1:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>1,]

group=c(rep("subtype1",normal),rep("subtype2",tumor))            
design <- model.matrix(~group)
y <- DGEList(counts=data,group=group)
y <- calcNormFactors(y)
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)
et <- exactTest(y,pair = c("subtype1","subtype2"))
topTags(et)
ordered_tags <- topTags(et, n=100000)

allDiff=ordered_tags$table
allDiff=allDiff[is.na(allDiff$FDR)==FALSE,]
diff=allDiff
newData=y$pseudo.counts

write.table(diff,file="edgerOut.xls",sep="\t",quote=F)
diffSig = diff[(diff$FDR < padj & (diff$logFC>foldChange | diff$logFC<(-foldChange))),]
write.table(diffSig, file="diffSig_mrna.xls",sep="\t",quote=F)
diffUp = diff[(diff$FDR < padj & (diff$logFC>foldChange)),]
write.table(diffUp, file="up.xls",sep="\t",quote=F)
diffDown = diff[(diff$FDR < padj & (diff$logFC<(-foldChange))),]
write.table(diffDown, file="down.xls",sep="\t",quote=F)

normalizeExp=rbind(id=colnames(newData),newData)
write.table(normalizeExp,file="normalizeExp.txt",sep="\t",quote=F,col.names=F) 
diffExp=rbind(id=colnames(newData),newData[rownames(diffSig),])
write.table(diffExp,file="diffmRNAExp.txt",sep="\t",quote=F,col.names=F)        




