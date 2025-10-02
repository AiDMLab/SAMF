library(stringr)
library("limma")
normalCount=114               
tumorCount=50              

setwd("D:\\data\\3_immunecells")         
mrna=read.table("mRNA.txt",sep="\t",header=T,check.names=F,row.names = 1)  
mrna <- t(mrna)
cancer_samples <- subset(mrna, str_sub(rownames(mrna),start=14,end=15) < 10)  
row.names(cancer_samples) <- rowdata
mrna <- as.data.frame(cancer_samples)

cluster <- read.table("ICIcluster2.txt", header=F, sep="\t", check.names=F, row.names=1)
colnames(cluster)[1] <- "subtype"
new_row_names <- substr(row.names(cluster), 1, 4) 
new_row_names <- paste0(new_row_names, ".", substr(row.names(cluster), 6, nchar(row.names(cluster)))) 
new_row_names1 <- paste0(new_row_names1, ".", substr(new_row_names, 9, nchar(new_row_names))) 
rownames(cluster) <- new_row_names1

sameSample=intersect(row.names(cluster), row.names(mrna))
mrna1 <- cbind(mrna[sameSample,], subtype=cluster[sameSample,])
mrna2 <- mrna1[order(mrna1$subtype),]
mrna3 <- mrna2[, -ncol(mrna2)]

rt <- as.matrix(t(mrna3))
exp=rt[,1:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

group=c(rep("subtype1",normalCount),rep("subtype2",tumorCount))
design <- model.matrix(~factor(group))
colnames(design)=levels(factor(group))
rownames(design)=colnames(data)
v <-voom(data, design = design, plot = F, save.plot = F)
out=v$E
out=rbind(ID=colnames(out),out)
write.table(out,file="uniq.symbol.txt",sep="\t",quote=F,col.names=F)        #?????ļ?

inputFile="uniq.symbol.txt"      
source("ICI12.CIBERSORT.R")     

seed <- 123  
set.seed(seed)
outTab=CIBERSORT("ref.txt", inputFile, perm=100, QN=TRUE)

outTab=outTab[outTab[,"P-value"]<0.05,]
outTab=as.matrix(outTab[,1:(ncol(outTab)-3)])
outTab=rbind(id=colnames(outTab), outTab)
outputFileName = paste0("CIBERSORT-Results_", seed,".txt")
write.table(outTab, file=outputFileName, sep="\t", quote=F, col.names=F)
