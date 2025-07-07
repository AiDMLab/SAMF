library(stringr)

setwd("D:\\data\\3_mirna_NEW")
data <- read.table("miRNAmatrix.txt",header=T,row.names = 1,check.names = F)  #2234�л���1039������
datat <- t(data)  #1039��������2234�л���

cancer_samples <- subset(datat, str_sub(rownames(datat),start=14,end=15) < 10)  #��������������985��������2234�л���
rowdata <- str_sub(rownames(cancer_samples),end=12)  
row.names(cancer_samples) <- rowdata
#write.csv(cancer_samples,"gene_expression.csv",row.names = T)
#data_csv <- read.csv("gene_expression.csv",header = T)

# �жϵ�һ���Ƿ����ظ�ֵ������һ���߼�����
#dup_rows <- duplicated(data_csv[, 1])
# ɸѡ���������ظ�ֵ����
#data <- data_csv[!dup_rows, ]

rname <- rownames(cancer_samples)
dup_rows <- duplicated(rname)
data_removal <- cancer_samples[!dup_rows, ]  #962������,2234

missing_ratio <- apply(data_removal, 1, function(x) sum(is.na(x))/length(x))
data_cleaned <- data_removal[missing_ratio <= 0.2, ]    #962������,2234

missing_ratio <- apply(data_cleaned, 2, function(x) sum(is.na(x))/length(x))
data_cleaned2 <- data_cleaned[,missing_ratio <= 0.2 ]    #962������,2234

write.csv(data_cleaned2,"mirna_removal.csv",row.names = T)

