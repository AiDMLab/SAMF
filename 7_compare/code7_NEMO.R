library(NEMO)
library(SNFtool)

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
# 对数据框进行转置
tdna <- t(dna_1)
trna <- t(rna_1)
ttmb <- t(tmb_1)
twei <- t(wei_1)
#进行NEMO算法聚类
omics.list = list(tdna, trna, ttmb, twei)
clustering = nemo.clustering(omics.list)
cluster = nemo.clustering(omics.list, num.clusters=2, num.neighbors=50)
write.table(cluster,file="ICIcluster2_NEMO.txt",sep="\t",quote=F,col.names=F)

affinity.graph = nemo.affinity.graph(omics.list, k = 20)

# 将每个元素替换为它的倒数
W_1 <- 1 / affinity.graph

# 找到所有有限值的最大值
max_value <- max(W_1[is.finite(W_1)], na.rm = TRUE)

# 将 Inf 替换为最大值加 10
W_1[is.infinite(W_1)] <- max_value + 10

library(survival)
library(survminer)
clusterFile="ICIcluster2_NEMO.txt"     #NEMO算法生成的免疫聚类文件
cliFile="clinical.csv"               #生存数据文件
setwd("D:\\data\\6_compare")     #设置工作目录

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

pdf(file="survival_NEMO.pdf",onefile = FALSE, 
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