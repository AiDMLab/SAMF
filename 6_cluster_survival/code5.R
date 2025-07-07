library(SNFtool)
library(ConsensusClusterPlus)
setwd("D:\\data\\6_cluster_survival")

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

dna1 <- as.data.frame(apply(dna_common, 2, function(x) (x - min(x)) / (max(x) - min(x))))
rna1 <- as.data.frame(apply(rna_common, 2, function(x) (x - min(x)) / (max(x) - min(x))))
tmb1 <- as.data.frame(apply(tmb_common, 2, function(x) (x - min(x)) / (max(x) - min(x))))
wei1 <- as.data.frame(apply(wei_common, 2, function(x) (x - min(x)) / (max(x) - min(x))))

dna2 <- as.matrix(dist(dna1))
rna2 <- as.matrix(dist(rna1))
tmb2 <- as.matrix(dist(tmb1))
wei2 <- as.matrix(dist(wei1))

similarity_matrix_dna <- affinityMatrix(dna2,K=20,sigma=0.5)
similarity_matrix_rna <- affinityMatrix(rna2,K=20,sigma=0.5)
similarity_matrix_tmb <- affinityMatrix(tmb2,K=20,sigma=0.5)
similarity_matrix_wei <- affinityMatrix(wei2,K=20,sigma=0.5)

n <- nrow(dna2)
sample_similarity <- apply(similarity_matrix_dna, 1, function(x) sum(x) / (n-1))
sample_weights <- 1 / sample_similarity
total_weight <- sum(sample_weights)
weights <- sample_weights/total_weight
dna3 <- similarity_matrix_dna
for (i in 1:n) {
  for (j in 1:n) {
    dna3[i,j] <- weights[i] * weights[j] * similarity_matrix_dna[i,j]
  }
}

n <- nrow(rna2)
sample_similarity <- apply(similarity_matrix_rna, 1, function(x) sum(x) / (n-1))
sample_weights <- 1 / sample_similarity
total_weight <- sum(sample_weights)
weights <- sample_weights/total_weight
rna3 <- similarity_matrix_rna
for (i in 1:n) {
  for (j in 1:n) {
    rna3[i,j] <- weights[i] * weights[j] * similarity_matrix_rna[i,j]
  }
}

n <- nrow(tmb2)
sample_similarity <- apply(similarity_matrix_tmb, 1, function(x) sum(x) / (n-1))
sample_weights <- 1 / sample_similarity
total_weight <- sum(sample_weights)
weights <- sample_weights/total_weight
tmb3 <- similarity_matrix_tmb
for (i in 1:n) {
  for (j in 1:n) {
    tmb3[i,j] <- weights[i] * weights[j] * similarity_matrix_tmb[i,j]
  }
}

n <- nrow(wei2)
sample_similarity <- apply(similarity_matrix_wei, 1, function(x) sum(x) / (n-1))
sample_weights <- 1 / sample_similarity
total_weight <- sum(sample_weights)
weights <- sample_weights/total_weight
wei3 <- similarity_matrix_wei
for (i in 1:n) {
  for (j in 1:n) {
    wei3[i,j] <- weights[i] * weights[j] * similarity_matrix_wei[i,j]
  }
}

W = SNF(list(dna3,rna3,tmb3,wei3), 20, 20)


maxK=9
results=ConsensusClusterPlus(W,
                             maxK=maxK,
                             reps=50,
                             pItem=0.8,
                             pFeature=1,
                             title="D:\\悦\\晓萌师姐论文\\2020级张晓萌\\data\\test_NEW",
                             seed=123456,
                             plot="png")

#输出结果
clusterNum=2        #分几类，根据判断标准判断
cluster=results[[clusterNum]][["consensusClass"]]
write.table(cluster,file="ICIcluster2.txt",sep="\t",quote=F,col.names=F)