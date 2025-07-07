
library("iClusterPlus")

setwd("D:\\data\\6_compare")
# 读取数据
dna <- read.csv("gene_expression_removal.csv",row.names = 1)
rna <- read.csv("mirna_removal.csv",row.names = 1)
tmb <- read.csv("TMB_removal.csv",row.names = 1)
wei <- read.csv("wei.csv",row.names = 1)
# 处理缺失值
wei[is.na(wei)] <- 0

rownames_dna <- rownames(dna)
rownames_rna <- rownames(rna)
rownames_tmb <- rownames(tmb)
rownames_wei <- rownames(wei)

common_rows <- intersect(intersect(intersect(rownames_dna, rownames_rna), rownames_tmb),rownames_wei)
dna_common <- dna[common_rows, ]
rna_common <- rna[common_rows, ]
tmb_common <- as.data.frame(tmb[common_rows, ])
rownames(tmb_common) <- common_rows
wei_common <- wei[common_rows, ]
#数据标准化
dna1 <- as.data.frame(apply(dna_common, 2, function(x) (x - min(x)) / (max(x) - min(x))))
rna1 <- as.data.frame(apply(rna_common, 2, function(x) (x - min(x)) / (max(x) - min(x))))
tmb1 <- as.data.frame(apply(tmb_common, 2, function(x) (x - min(x)) / (max(x) - min(x))))
wei1 <- as.data.frame(apply(wei_common, 2, function(x) (x - min(x)) / (max(x) - min(x))))
#计算欧氏距离矩阵
dna2 <- as.matrix(dist(dna1))
rna2 <- as.matrix(dist(rna1))
tmb2 <- as.matrix(dist(tmb1))
wei2 <- as.matrix(dist(wei1))

# 删除全是NA的列
dna_1 <- dna1[, apply(dna1, 2, function(x) !all(is.na(x)))]
rna_1 <- rna1[, apply(rna1, 2, function(x) !all(is.na(x)))]
tmb_1 <- tmb1
wei_1 <- wei1[, apply(wei1, 2, function(x) !all(is.na(x)))]

# 将数据转换为矩阵或数据框
dna_1 <- as.matrix(dna_1)
rna_1 <- as.matrix(rna_1)
tmb_1 <- as.matrix(tmb_1)
wei_1 <- as.matrix(wei_1)


iClusterBayes_result <- iClusterBayes(dna_1,rna_1,tmb_1,wei_1,
              type = c("gaussian","gaussian","gaussian","gaussian"),K=1,n.burnin=1000,n.draw=1200,
              prior.gamma=rep(0.1,6),sdev=0.5,beta.var.scale=1,thin=1,pp.cutoff=0.5)
iClusterBayes_result
cluster <- iClusterBayes_result$cluster
result_matrix <- cbind(common_rows, cluster)
result_matrix
write.table(result_matrix,file="ICIcluster2_iClusterBayes.txt",sep="\t",quote=F,col.names=F, row.names = FALSE)

library(survival)
library(survminer)
clusterFile="ICIcluster2_iClusterBayes.txt"     #iClusterBayes算法生成的免疫聚类文件
cliFile="clinical.csv"               #生存数据文件
setwd("D:\\data\\6_compare")      #设置工作目录

#读取输入文件
cluster=read.table(clusterFile, header=F, sep="\t", check.names=F, row.names=1)
rownames(cluster)=gsub("(.*?)\\_(.*?)", "\\2", rownames(cluster))
cli=read.table(cliFile, header=T, sep=",", check.names=F, row.names=1)
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

pdf(file="survival_iClusterBayes.pdf",onefile = FALSE, 
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
