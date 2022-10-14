library(dplyr)
library(stringr)

##############
#MouseBrainV3
##############
rm(list = ls())
#Step1: create a barcode data
setwd("D:/9.CellSeg_MouseBrainV3_0519/")
base_calling_data<-read.table('1.Raw/all_basecalling_data.txt',sep ='\t' ,
                              comment.char='', quote='')
base_calling_data <- base_calling_data[,c(2,4,5)]
names(base_calling_data)=c("barcode","x","y")

barcode<-read.table('1.Raw/mouse brain barcode2.txt',sep ='\t' ,
                    comment.char='', quote='')
barcode <- barcode[,1:2]
names(barcode)=c("barcode","gene")
barcode$gene <- str_to_title(barcode$gene) #小鼠基因，首字母大写
head(barcode)

Bar_data<-merge(base_calling_data,barcode,by="barcode")
head(Bar_data)
table(Bar_data$gene)
write.csv(Bar_data,'Barcode_data.csv',row.names = F, quote = F)

#Step2: create a Freq data(作为ClusterMap的输入)
df <- data.frame(
  Index = 1:length(Bar_data[,1]),
  spot_location_1 = Bar_data$y,
  spot_location_2 = Bar_data$x,
  spot_location_3 = rep(0,length(Bar_data[,1])),
  gene = factor(Bar_data$gene,labels = 1:36),
  gene_name = Bar_data$gene
) %>% arrange(gene) #根据基因编号从小到大排序
sort(table(df$gene_name),decreasing = T) #看看哪个基因最多
write.csv(df,'Freq.csv',row.names = F, quote = F)


