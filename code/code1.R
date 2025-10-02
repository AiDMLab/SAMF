library(stringr)

#rna
setwd("D:\\data\\1_data_processing")
data <- read.table("symbol_1.txt",header=T,row.names = 1,check.names = F)  
datat <- t(data)  

cancer_samples <- subset(datat, str_sub(rownames(datat),start=14,end=15) < 10) 
rowdata <- str_sub(rownames(cancer_samples),end=12)  
row.names(cancer_samples) <- rowdata
write.csv(cancer_samples,"gene_expression.csv",row.names = T)
data_csv <- read.csv("gene_expression.csv",header = T)

dup_rows <- duplicated(data_csv[, 1])
data <- data_csv[!dup_rows, ]

rname <- rownames(cancer_samples)
dup_rows <- duplicated(rname)
data_removal <- cancer_samples[!dup_rows, ] 

missing_ratio <- apply(data_removal, 1, function(x) sum(is.na(x))/length(x))
data_cleaned <- data_removal[missing_ratio <= 0.2, ] 

missing_ratio <- apply(data_cleaned, 2, function(x) sum(is.na(x))/length(x))
data_cleaned2 <- data_cleaned[,missing_ratio <= 0.2 ] 
write.csv(data_cleaned2,"gene_expression_removal.csv",row.names = T)

#mirna
data <- read.table("miRNAmatrix.txt",header=T,row.names = 1,check.names = F) 
datat <- t(data)
cancer_samples <- subset(datat, str_sub(rownames(datat),start=14,end=15) < 10)  
rowdata <- str_sub(rownames(cancer_samples),end=12)  
row.names(cancer_samples) <- rowdata

rname <- rownames(cancer_samples)
dup_rows <- duplicated(rname)
data_removal <- cancer_samples[!dup_rows, ]  

missing_ratio <- apply(data_removal, 1, function(x) sum(is.na(x))/length(x))
data_cleaned <- data_removal[missing_ratio <= 0.2, ]    

missing_ratio <- apply(data_cleaned, 2, function(x) sum(is.na(x))/length(x))
data_cleaned2 <- data_cleaned[,missing_ratio <= 0.2 ]    

write.csv(data_cleaned2,"mirna_removal.csv",row.names = T)

#tmb
data <- read.table("TMB.txt",header=T,row.names = 1,check.names = F)  

ddata <- as.matrix(data)
rname <- rownames(ddata)
dup_rows <- duplicated(rname)
data_removal <- ddata[!dup_rows, ] 

write.csv(data_removal,"TMB_removal.csv",row.names = T)

#microbial
data1 <- read.csv("data1.csv",header=T,row.names = 1,check.names = F) 
data1t <- t(data1)  
data2 <- read.csv("data3.csv",header=T,check.names = F)  

dup_rows <- duplicated(data2[, 1])
data2n <- data2[!dup_rows, ]  
rownames(data2n) <- data2n[,1]
data2n <- data2n[-1]

data <- cbind(data1t, data2n)    

missing_ratio <- apply(data, 2, function(x) sum(is.na(x))/length(x))
data_cleaned <- data[,missing_ratio <= 0.2]   

missing_ratio <- apply(data_cleaned, 1, function(x) sum(is.na(x))/length(x))
data_cleaned2 <- data_cleaned[missing_ratio <= 0.2, ]    
write.csv(data_cleaned2,"wei.csv",row.names = T)