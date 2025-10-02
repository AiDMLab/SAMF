library(caret)     
library(edgeR)
library(glmnet)
setwd("D:\\data\\4_Transcriptionfactor&logics")

mirna <- read.csv("mirna_removal.csv",row.names = 1)
tmb <- read.csv("TMB_removal.csv",row.names = 1)
tmb <- as.matrix(tmb)
tmb <- avereps(tmb)
tmb2 <- ifelse(tmb>=10,1,0)

sameSample <- intersect(row.names(mirna), row.names(tmb2))
data <- cbind(mirna[sameSample,], tmb=tmb2[sameSample,])
data <- data[,c(ncol(data),1:(ncol(data)-1))]

cli=read.csv("clinical.csv",row.names = 1)     
sameSample=intersect(row.names(cli),row.names(data))
data=data[sameSample,]
cli=cli[sameSample,]

inTrain=createDataPartition(y=data[,1],p=0.6,list=F)     
train=data[inTrain,]
test=data[-inTrain,]

trainOut=cbind(id=row.names(train),train)
testOut=cbind(id=row.names(test),test)
write.table(trainOut,file="train.txt",sep="\t",quote=F,row.names=F)
write.table(testOut,file="test.txt",sep="\t",quote=F,row.names=F)

trainCli=cli[row.names(train),]
testCli=cli[row.names(test),]
trainFlag=cbind(trainCli,flag="Train")
testFlag=cbind(testCli,flag="Test")
newTable=rbind(trainFlag,testFlag)
cliStatOut=data.frame()
for(i in 1:(ncol(newTable)-1)){
  nameStat=colnames(newTable)[i]
  tableStat=table(newTable[,c(nameStat,"flag")])
  tableStatSum=cbind(Total=rowSums(tableStat),tableStat)
  tableStatRatio=prop.table(tableStatSum,2)
  tableStatRatio=round(tableStatRatio*100,2)
  tableStatPaste=paste(tableStatSum,"(",tableStatRatio,"%)",sep="")
  tableStatOut=matrix(tableStatPaste,ncol=3,dimnames=dimnames(tableStatSum))
  pStat=chisq.test(tableStat[row.names(tableStat)!="unknow",])
  pValueStat=round(pStat$p.value,4)
  pValueCol=c(pValueStat,rep(" ",(nrow(tableStatOut)-1)) )
  tableStatOut=cbind(Covariates=nameStat,Type=row.names(tableStatOut),tableStatOut,Pvalue=pValueCol)
  cliStatOut=rbind(cliStatOut,tableStatOut)
}
write.table(cliStatOut,file="cliStat1.xls",sep="\t",quote=F,row.names=F)

library(limma)
library(pheatmap)

inputFile="train_NEW.txt"                                        
fdrFilter=0.01                                                
foldChange=1.5                                                 
logFCfilter=log2(foldChange)                                 

rt=read.table(inputFile,sep="\t",header=T,check.names=F,row.names=1)
rtLow=rt[rt$tmb==0,]
rtHigh=rt[rt$tmb==1,]
conNum=nrow(rtLow)
treatNum=nrow(rtHigh)
data=rbind(rtLow,rtHigh)
data=t(data[,(2:ncol(data))])


outTab=data.frame()
grade=c(rep(1,conNum),rep(2,treatNum))
samplePercent=ncol(data)*0.1
for(i in row.names(data)){
  if(length(data[i,(data[i,]==0)])>samplePercent){
    next
  }
  geneName=unlist(strsplit(i,"\\|",))[1]
  geneName=gsub("\\/", "_", geneName)
  rt=rbind(expression=data[i,],grade=grade)
  rt=as.matrix(t(rt))
  wilcoxTest<-wilcox.test(expression ~ grade, data=rt)
  conGeneMeans=mean(data[i,1:conNum])
  treatGeneMeans=mean(data[i,(conNum+1):ncol(data)])
  logFC=log2(treatGeneMeans)-log2(conGeneMeans)
  pvalue=wilcoxTest$p.value
  conMed=median(data[i,1:conNum])
  treatMed=median(data[i,(conNum+1):ncol(data)])
  diffMed=treatMed-conMed
  if( ((logFC>0) & (diffMed>0)) | ((logFC<0) & (diffMed<0)) ){  
    outTab=rbind(outTab,cbind(gene=i,conMean=conGeneMeans,treatMean=treatGeneMeans,logFC=logFC,pValue=pvalue))
  }
}
pValue=outTab[,"pValue"]
fdr=p.adjust(as.numeric(as.vector(pValue)),method="fdr")
outTab=cbind(outTab,fdr=fdr)

write.table(outTab,file="all.xls",sep="\t",row.names=F,quote=F)
outDiff=outTab[( abs(as.numeric(as.vector(outTab$logFC)))>logFCfilter & as.numeric(as.vector(outTab$fdr))<fdrFilter),]
write.table(outDiff,file="diff.xls",sep="\t",row.names=F,quote=F)
heatmap=data[as.vector(outDiff[,1]),]
heatmap=t(heatmap)
heatmap=log2(heatmap+1)
TMB=c(rep(0,conNum),rep(1,treatNum))
heatmap=cbind(TMB,heatmap)
heatmap=rbind(ID=colnames(heatmap),heatmap)
write.table(heatmap,file="diffMirnaExp.txt",sep="\t",col.names=F,quote=F)


library(glmnet)
library(pROC)

train=read.table("diffMirnaExp.txt",header=T,sep="\t",check.names=F,row.names=1)  
test=read.table("test.txt",header=T,sep="\t",check.names=F,row.names=1)            
test[,c(2:ncol(test))]=log2(test[,c(2:ncol(test))]+1)
test=test[,colnames(train)]
rt=rbind(train,test)

x=as.matrix(train[,c(2:ncol(train))])
y=train[,1]
fit=glmnet(x, y, family = "binomial", alpha=1)
cvfit=cv.glmnet(x, y, family="binomial", alpha=1,type.measure='auc',nfolds = 10)

coef=coef(fit, s = cvfit$lambda.min)
index=which(coef != 0)
actCoef=coef[index]
lassoGene=row.names(coef)[index]
coef=cbind(Gene=lassoGene,Coef=actCoef)
write.table(coef,file="coef.txt",sep="\t",quote=F,row.names=F)
write.table(lassoGene[-1],file="miRNAlist.txt",sep="\t",quote=F,row.names=F,col.names=F)

threshold=0.5
bioStat=function(x=null, y=null, dataType=null){
  pred=predict(cvfit,newx=x,s=cvfit$lambda.min,type = 'response')
  outTab=cbind(id=row.names(pred),TMB=as.character(y),pred)
  colnames(outTab)=c("id","TMB","value")
  sameGene=intersect(colnames(x),lassoGene)
  outTab=cbind(outTab,x[,sameGene])
  write.table(outTab,file=paste0("value.",dataType,".txt"),sep="\t",quote=F,row.names=F)
  
  pred_new=as.integer(pred>threshold) 
  tab=table(pred_new,y)
  TP=tab[2,2];TN=tab[1,1];FP=tab[2,1];FN=tab[1,2]
  Accuracy=round((TP+TN)/(TP+FN+FP+TN),4)
  SE=round(TP/(TP+FN),4)
  SP=round(TN/(TN+FP),4)
  PPV=round(TP/(TP+FP),4)
  NPV=round(TN/(TN+FN),4)
  
  rocobj1=roc(y, as.vector(pred))
  AUC=auc(rocobj1)
  AUC=round(AUC,4)
  return(c(SE,SP,PPV,NPV,Accuracy,AUC))
}

trainValue=bioStat(x=x,y=y,dataType="train")
testX=as.matrix(test[,c(2:ncol(test))])
testY=test[,1]
testValue=bioStat(x=testX,y=testY,dataType="test")
totalX=as.matrix(rt[,c(2:ncol(rt))])
totalY=rt[,1]
totalValue=bioStat(x=totalX,y=totalY,dataType="all")

statTab=rbind(Train=trainValue,Test=testValue,Total=totalValue)
colnames(statTab)=c("SE","SP","PPV","NPV","Accuracy","AUC")
statTab=rbind(ID=colnames(statTab),statTab)
write.table(statTab,file="accuracyStat.xls",sep="\t",col.names=F,quote=F)