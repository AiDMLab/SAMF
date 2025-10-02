library(SNFtool)
library(ConsensusClusterPlus)
setwd("D:\\data\\2_SAMF")

dna <- read.csv("gene_expression_removal.csv",row.names = 1)
rna <- read.csv("mirna_removal.csv",row.names = 1)
tmb <- read.csv("TMB_removal.csv",row.names = 1)
wei <- read.csv("wei.csv",row.names = 1)
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

dna11 <- read.csv("dna_aligned.csv",row.names = 1)
rna11 <- read.csv("rna_aligned.csv",row.names = 1)
tmb11 <- read.csv("tmb_aligned.csv",row.names = 1)
wei11 <- read.csv("wei_aligned.csv",row.names = 1)
cross11 <- read.csv("fus_aligned.csv",row.names = 1)

dna2 <- as.matrix(dist(dna11))
rna2 <- as.matrix(dist(rna11))
tmb2 <- as.matrix(dist(tmb11))
wei2 <- as.matrix(dist(wei11))

similarity_matrix_dna <- affinityMatrix(dna2,K=20,sigma=0.5)
similarity_matrix_rna <- affinityMatrix(rna2,K=20,sigma=0.5)
similarity_matrix_tmb <- affinityMatrix(tmb2,K=20,sigma=0.5)
similarity_matrix_wei <- affinityMatrix(wei2,K=20,sigma=0.5)

W = SNF(list(similarity_matrix_dna,similarity_matrix_rna,similarity_matrix_tmb,similarity_matrix_wei), 20, 20)

maxK=9
results=ConsensusClusterPlus(W,
                             maxK=maxK,
                             reps=50,
                             pItem=0.8,
                             pFeature=1,
                             title="D:\\data\\2_SAMF\\picture",
                             seed=123456,
                             plot="png")

clusterNum=2        
cluster=results[[clusterNum]][["consensusClass"]]
write.table(cluster,file="ICIcluster2.txt",sep="\t",quote=F,col.names=F)

library(survival)
library(survminer)
clusterFile="ICIcluster2.txt"    
cliFile="clinical.csv"            
setwd("D:\\data\\2_SAMF")     

cluster=read.table(clusterFile, header=F, sep="\t", check.names=F, row.names=1)
rownames(cluster)=gsub("(.*?)\\_(.*?)", "\\2", rownames(cluster))
cli=read.table(cliFile, header=T, sep=",", check.names=F, row.names=1)
cli$futime=cli$futime/365

sameSample=intersect(row.names(cluster), row.names(cli))
rt=cbind(cli[sameSample,], ICIcluster=cluster[sameSample,])
letter=c("A","B","C","D","E","F","G")
uniqClu=levels(factor(rt$ICIcluster))
rt$ICIcluster=letter[match(rt$ICIcluster, uniqClu)]
print(rt$ICIcluster)

diff=survdiff(Surv(futime, fustat) ~ICIcluster,data = rt)   
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)   
pValue=format(pValue, scientific = TRUE)   

fit <- survfit(Surv(futime, fustat) ~ ICIcluster, data = rt)    #
summary(fit) 

pdf(file="survival.pdf",onefile = FALSE, 
    width = 5.5,          
    height =5)            
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
