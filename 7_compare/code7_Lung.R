library(SNFtool)
setwd("D:\\data\\7_compare")

dna=read.table("LUNG_Gene_Expression.txt")  #需要转置dna <- t(dna)
dna <- t(dna) 
rna=read.table("LUNG_Mirna_Expression.txt")
rna <- t(rna)
tmb=read.table("LUNG_Methy_Expression.txt")
tmb <- t(tmb)

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

dna3 <- affinityMatrix(dna2,K=20,sigma=0.5)
rna3 <- affinityMatrix(rna2,K=20,sigma=0.5)
tmb3 <- affinityMatrix(tmb2,K=20,sigma=0.5)
wei3 <- affinityMatrix(wei2,K=20,sigma=0.5)

W = SNF(list(dna3,rna3,tmb3,wei3), 20, 20)
labels = spectralClustering(W,2)
results <- list(name=common_rows,type=labels)
write.table(results,file="ICIcluster2_Lung.txt",sep="\t",quote=F,col.names=F)

cluster=read.table("ICIcluster2_Lung.txt", header=F, sep="\t", check.names=F, row.names=2)
cluster <- cluster[,-1,drop=F]
write.table(cluster,file="ICIcluster2_Lung.txt",sep="\t",quote=F,col.names=F)