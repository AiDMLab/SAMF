library(SNFtool)
library(ConsensusClusterPlus)
setwd("D:\\data\\1_data_processing")

dna <- read.csv("gene_expression_removal_NEW.csv",row.names = 1)
rna <- read.csv("mirna_removal_NEW.csv",row.names = 1)
tmb <- read.csv("TMB_removal_NEW.csv",row.names = 1)
wei <- read.csv("wei_NEW.csv",row.names = 1)
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

write.csv(dna1, file = "dna1.csv")
write.csv(rna1, file = "rna1.csv")
write.csv(tmb1, file = "tmb1.csv")
write.csv(wei1, file = "wei1.csv")
