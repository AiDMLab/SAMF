library(stringr)

setwd("D:\\data\\4_TMB_NEW")
data <- read.table("TMB.txt",header=T,row.names = 1,check.names = F)  #941������

ddata <- as.matrix(data)
rname <- rownames(ddata)
dup_rows <- duplicated(rname)
data_removal <- ddata[!dup_rows, ]  #941

write.csv(data_removal,"TMB_removal.csv",row.names = T)

