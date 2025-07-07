library(stringr)

setwd("D:\\data\\1_microorganism")
data1 <- read.csv("data1.csv",header=T,row.names = 1,check.names = F)  #2852�л���203������
data1t <- t(data1)  #203��������2852��ϸ��
data2 <- read.csv("data3.csv",header=T,check.names = F)  #208��������68�л���


# �жϵ�һ���Ƿ����ظ�ֵ������һ���߼�����
dup_rows <- duplicated(data2[, 1])
# ɸѡ���������ظ�ֵ����
data2n <- data2[!dup_rows, ]   #203��68����
rownames(data2n) <- data2n[,1]
data2n <- data2n[-1]


data <- cbind(data1t, data2n)    #203��������2920��΢����


missing_ratio <- apply(data, 2, function(x) sum(is.na(x))/length(x))
data_cleaned <- data[,missing_ratio <= 0.2]    #203������,16��΢����

missing_ratio <- apply(data_cleaned, 1, function(x) sum(is.na(x))/length(x))
data_cleaned2 <- data_cleaned[missing_ratio <= 0.2, ]    #193������,16

write.csv(data_cleaned2,"wei.csv",row.names = T)

