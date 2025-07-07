library(survival)
library(survminer)
clusterFile="ICIcluster2.txt"     #code5生成的免疫聚类文件
#clusterFile="ICIcluster2_SNF_Lung_NEW.txt"     #code7生成的免疫聚类文件
cliFile="clinical.csv"               #生存数据文件
setwd("D:\\data\\6_cluster_survival")      #设置工作目录

#读取输入文件
cluster=read.table(clusterFile, header=F, sep="\t", check.names=F, row.names=1)
rownames(cluster)=gsub("(.*?)\\_(.*?)", "\\2", rownames(cluster))
cli=read.table(cliFile, header=T, sep=",", check.names=F, row.names=1)
#colnames(cli)=c("age","futime","sex","site","fustat")
#cli[,c(2)]<-as.numeric(unlist(cli[,c(2)]),na.rm=T)
cli$futime=cli$futime/365

#数据合并
sameSample=intersect(row.names(cluster), row.names(cli))
rt=cbind(cli[sameSample,], ICIcluster=cluster[sameSample,])
letter=c("A","B","C","D","E","F","G")
uniqClu=levels(factor(rt$ICIcluster))
rt$ICIcluster=letter[match(rt$ICIcluster, uniqClu)]
print(rt$ICIcluster)

diff=survdiff(Surv(futime, fustat) ~ICIcluster,data = rt)    #survdiff函数计算高低风险组之间的差异
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)     #p保留4位有效数字
pValue=format(pValue, scientific = TRUE)    #p值用科学记数法表示

fit <- survfit(Surv(futime, fustat) ~ ICIcluster, data = rt)    #每个时间段的病人有多少个
summary(fit)    #查看五年生存率

#绘制生存曲线

pdf(file="survival.pdf",onefile = FALSE, #file="survival_Lung_NEW.pdf"或file="survival_SNF_NEW.pdf"
    width = 5.5,             #图片的宽度
    height =5)             #图片的高度
ggsurvplot(fit, 
           data=rt,
           conf.int=FALSE,
           pval=paste0("p=",pValue),
           pval.size=4,
           risk.table=TRUE,
           legend.labs=c("Cluster 1", "Cluster 2"),
           legend.title="Cluster",
           xlab="Time(years)",
           break.time.by = 1,
           risk.table.title="",
           palette=c("red", "blue"),
           risk.table.height=.25)
dev.off()


#绘制code7SNF算法的生存曲线
# pdf(file="survival.pdf",onefile = FALSE, #file="survival_Lung_NEW.pdf"或file="survival_SNF_NEW.pdf"
#     width = 5.5,             #图片的宽度
#     height =5)             #图片的高度
# ggsurvplot(fit, 
#            data=rt,
#            conf.int=FALSE,
#            pval=paste0("p=",pValue),
#            pval.size=4,
#            risk.table=TRUE,
#            legend.labs=c("Cluster 1", "Cluster 2"),
#            legend.title="Cluster",
#            xlab="Time(years)",
#            break.time.by = 1,
#            risk.table.title="",
#            palette=c("red", "blue"),
#            risk.table.height=.25)
# dev.off()

