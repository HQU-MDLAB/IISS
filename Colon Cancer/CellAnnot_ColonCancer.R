library(dplyr)

rm(list = ls())
options(stringsAsFactors = F)
setwd("D:/11.CellSeg_ColonV2_0613/3.cellpose/2.概率模型_0722/")

##########################函数主体##############################################
Find_Multiple_Maximum_Index <- function(result){
  index <- c()
  for (c in 1:ncol(result)){
    cell <- result[,c]
    GreaterZeroNum <- length(cell[which(cell > 0)])
    MaxValueNum <- length(which(cell == max(cell)))
    if (MaxValueNum == GreaterZeroNum) {
      if (max(cell) != 1) {index=c(index,colnames(result[c]))}
    }
  }
  return(index)
}

Scoring <- function(mx_path){
  
  mx <- read.csv(mx_path,header = T,check.names = F)
  rownames(mx) <- mx$gene
  mx <- mx[,3:ncol(mx)]
  
  ### 第一步：定义细胞类型marker
  Stromal <- c("ACTA2","CAV1","ENG","TIMP2","DES","TNC","THY1","ITGB1","TIMP1")
  Epithelial <- c("PPIB","REG4","LDHA","S100A4","SOX4","CD63","HSPB1","PIGR")
  Normal <- c("TAGLN")
  Malignancy <- c("CD24","CEACAM5")
  B_Cell <- c("CD79A","CD81")
  T_Cell <- c("CD3D","CD8A","CD7","CD3E","CD3G","CD8B","CD56")
  Mast_Cell <- c("CPA3","MS4A2","CTSG")
  Myeloids <- c("SPP1","CD68","CD80","S100A8")
  NfkB_Pathway_related <- c("NFKBIA")
  RTK_Pathway_related <- c("PDGFRB")
  Wnt_Pathway_related <- c("OLFM4","WNT5A")
  
  CellTypeList <- list(Stromal,Epithelial,Normal,Malignancy,B_Cell,T_Cell,Mast_Cell,Myeloids,NfkB_Pathway_related,RTK_Pathway_related,Wnt_Pathway_related)
  
  ### 第二步：计算细胞类型得分
  mx <- as.data.frame(apply(mx,2,prop.table)) #将数值转换成每列百分比
  result <- data.frame() #创建一个空data frame
  
  for (i in 1:length(CellTypeList)){
    gene <- CellTypeList[[i]] #i这个细胞类型对应的基因；字符型
    mx_subset <- mx[gene,] #只保留上述基因
    mx_subset <- apply(mx_subset,2,sum) #细胞类型占比总和
    mx_subset[is.na(mx_subset)] <-  0 #缺失值设置为0
    mx_subset <- as.data.frame(t(mx_subset)) #data frame转置以方便下述result表格生成
    result <-  rbind(result,mx_subset) #result为结果打分结果表格
  }
  rownames(result) <- c("Stromal","Epithelial","Normal","Malignancy","B_Cell","T_Cell","Mast_Cell","Myeloids","NfkB_Pathway_related","RTK_Pathway_related","Wnt_Pathway_related")
  
  ### 第三步：找出打分不确定（即在这个细胞内，同时有多个类型打分值相同）的细胞
  index <- Find_Multiple_Maximum_Index(result = result)
  
  ### 第四步：判断最大值
  CellType <- apply(result,2,FUN = function(t) rownames(result)[which.max(t)])
  CellType <- as.data.frame(CellType)
  
  ### 第五步：给第三步挑选出的细胞定义为"Unknown"类型
  CellType[index,] <- "Unknown"
  
  return(CellType)
}

###############################示例运行#########################################
mx_path <- file.path("Matrix.csv")
CellType <- Scoring(mx_path = mx_path)


##############空间可视化########################################################
library(ggplot2)

### 导入python生成的mask质心位置
centroids <- read.csv("centroid.csv",header = T) 
centroids <- centroids[,-1] #删除第一列索引值
head(centroids)

### visualizatio：用于细胞类型空间位置可视化的data frame
visualization = data.frame(centroids,celltype=0,cell=rownames(centroids)) #新建celltype列
visualization = visualization[,c('cell','row','column','celltype')]
visualization[rownames(CellType),]$celltype <- CellType$CellType #特定行重新赋值
visualization = filter(visualization, !(celltype == "0" | celltype == "Unknown"))
head(visualization,100)

write.csv(visualization, "CellMeta.csv", row.names = F)

### 0：没有信号点的细胞
### Unknown：打分不明确的细胞

### 总图
ggplot(data = visualization,mapping = aes(x=column, y=row, color = celltype))+
  geom_point(size=0.7)+
  scale_y_reverse()+
  theme_void()+
  theme(panel.background = element_rect(fill = 'black', color = 'black'),
        legend.position = "none",
        legend.text = element_text(face = 'italic', color = "white", size = 20),
        legend.title = element_text(face = 'italic', color = "white", size = 25, hjust = 0.5))+
  facet_wrap( ~ celltype)+
  theme(strip.text = element_text(face = 'bold', color = "black", size = rel(1.5)),
        strip.background = element_rect(fill = 'lightblue', color = 'black'))+
  scale_colour_manual(values = c("#ff0029",'#ff6000','#ffea00','#8aff00','#00ff00','#00ff89','#00ecff','#0061ff','#2a00ff','#b400ff','#ff00bf'),
                      limits = c("Malignancy","Epithelial","Normal","Stromal","B_Cell","T_Cell","Mast_Cell","Myeloids","NfkB_Pathway_related","RTK_Pathway_related","Wnt_Pathway_related"))
ggsave(last_plot(),filename = "CellTypeMap_facet.png", device = "png", height = 20,width = 25,dpi = 300)

ggplot(data = visualization,mapping = aes(x=column, y=row, color = celltype))+
  geom_point(size=0.7)+
  scale_y_reverse()+
  theme_void()+
  theme(panel.background = element_rect(fill = 'black', color = 'black'),
        legend.position = "none",
        legend.text = element_text(face = 'italic', color = "white", size = 10),
        legend.title = element_text(face = 'italic', color = "white", size = 20, hjust = 0.5))+
  scale_colour_manual(values = c("#ff0029",'#ff6000','#ffea00','#8aff00','#00ff00','#00ff89','#00ecff','#0061ff','#2a00ff','#b400ff','#ff00bf'),
                      limits = c("Malignancy","Epithelial","Normal","Stromal","B_Cell","T_Cell","Mast_Cell","Myeloids","NfkB_Pathway_related","RTK_Pathway_related","Wnt_Pathway_related"))
ggsave(last_plot(),filename = "CellTypeMap.png", device = "png", height = 10,width = 10,dpi = 300)

### 单类型分布
### Malignancy
visualization = data.frame(centroids,celltype=0) #新建celltype列
visualization[rownames(CellType),]$celltype <- CellType$CellType #特定行重新赋值
visualization = filter(visualization, !(celltype == "0" | celltype == "Unknown"))
visualization = filter(visualization, celltype == "Malignancy")
head(visualization,100)

ggplot(data = visualization,mapping = aes(x=column, y=row, color = celltype))+
  geom_point(size=1,colour="#ff0029")+
  scale_y_reverse()+
  theme_void()+
  theme(panel.background = element_rect(fill = 'black', color = 'black'))
ggsave(last_plot(),filename = "Malignancy.png", device = "png", height = 10,width = 10,dpi = 300)

### Epithelial
visualization = data.frame(centroids,celltype=0) #新建celltype列
visualization[rownames(CellType),]$celltype <- CellType$CellType #特定行重新赋值
visualization = filter(visualization, !(celltype == "0" | celltype == "Unknown"))
visualization = filter(visualization, celltype == "Epithelial")
head(visualization,100)

ggplot(data = visualization,mapping = aes(x=column, y=row, color = celltype))+
  geom_point(size=1,colour="#ff6000")+
  scale_y_reverse()+
  theme_void()+
  theme(panel.background = element_rect(fill = 'black', color = 'black'))
ggsave(last_plot(),filename = "Epithelial.png", device = "png", height = 10,width = 10,dpi = 300)

### Normal
visualization = data.frame(centroids,celltype=0) #新建celltype列
visualization[rownames(CellType),]$celltype <- CellType$CellType #特定行重新赋值
visualization = filter(visualization, !(celltype == "0" | celltype == "Unknown"))
visualization = filter(visualization, celltype == "Normal")
head(visualization,100)

ggplot(data = visualization,mapping = aes(x=column, y=row, color = celltype))+
  geom_point(size=1,colour="#ffea00")+
  scale_y_reverse()+
  theme_void()+
  theme(panel.background = element_rect(fill = 'black', color = 'black'))
ggsave(last_plot(),filename = "Normal.png", device = "png", height = 10,width = 10,dpi = 300)

### Stromal
visualization = data.frame(centroids,celltype=0) #新建celltype列
visualization[rownames(CellType),]$celltype <- CellType$CellType #特定行重新赋值
visualization = filter(visualization, !(celltype == "0" | celltype == "Unknown"))
visualization = filter(visualization, celltype == "Stromal")
head(visualization,100)

ggplot(data = visualization,mapping = aes(x=column, y=row, color = celltype))+
  geom_point(size=1,colour="#8aff00")+
  scale_y_reverse()+
  theme_void()+
  theme(panel.background = element_rect(fill = 'black', color = 'black'))
ggsave(last_plot(),filename = "Stromal.png", device = "png", height = 10,width = 10,dpi = 300)

### B_Cell
visualization = data.frame(centroids,celltype=0) #新建celltype列
visualization[rownames(CellType),]$celltype <- CellType$CellType #特定行重新赋值
visualization = filter(visualization, !(celltype == "0" | celltype == "Unknown"))
visualization = filter(visualization, celltype == "B_Cell")
head(visualization,100)

ggplot(data = visualization,mapping = aes(x=column, y=row, color = celltype))+
  geom_point(size=2,colour="#00ff00")+
  scale_y_reverse()+
  theme_void()+
  theme(panel.background = element_rect(fill = 'black', color = 'black'))
ggsave(last_plot(),filename = "B_Cell.png", device = "png", height = 10,width = 10,dpi = 300)

### T_Cell
visualization = data.frame(centroids,celltype=0) #新建celltype列
visualization[rownames(CellType),]$celltype <- CellType$CellType #特定行重新赋值
visualization = filter(visualization, !(celltype == "0" | celltype == "Unknown"))
visualization = filter(visualization, celltype == "T_Cell")
head(visualization,100)

ggplot(data = visualization,mapping = aes(x=column, y=row, color = celltype))+
  geom_point(size=2,colour="#00ff89")+
  scale_y_reverse()+
  theme_void()+
  theme(panel.background = element_rect(fill = 'black', color = 'black'))
ggsave(last_plot(),filename = "T_Cell.png", device = "png", height = 10,width = 10,dpi = 300)

### Mast_Cell
visualization = data.frame(centroids,celltype=0) #新建celltype列
visualization[rownames(CellType),]$celltype <- CellType$CellType #特定行重新赋值
visualization = filter(visualization, !(celltype == "0" | celltype == "Unknown"))
visualization = filter(visualization, celltype == "Mast_Cell")
head(visualization,100)

ggplot(data = visualization,mapping = aes(x=column, y=row, color = celltype))+
  geom_point(size=2,colour="#00ecff")+
  scale_y_reverse()+
  theme_void()+
  theme(panel.background = element_rect(fill = 'black', color = 'black'))
ggsave(last_plot(),filename = "Mast_Cell.png", device = "png", height = 10,width = 10,dpi = 300)

### Myeloids
visualization = data.frame(centroids,celltype=0) #新建celltype列
visualization[rownames(CellType),]$celltype <- CellType$CellType #特定行重新赋值
visualization = filter(visualization, !(celltype == "0" | celltype == "Unknown"))
visualization = filter(visualization, celltype == "Myeloids")
head(visualization,100)

ggplot(data = visualization,mapping = aes(x=column, y=row, color = celltype))+
  geom_point(size=2,colour="#0061ff")+
  scale_y_reverse()+
  theme_void()+
  theme(panel.background = element_rect(fill = 'black', color = 'black'))
ggsave(last_plot(),filename = "Myeloids.png", device = "png", height = 10,width = 10,dpi = 300)

### NfkB_Pathway_related
visualization = data.frame(centroids,celltype=0) #新建celltype列
visualization[rownames(CellType),]$celltype <- CellType$CellType #特定行重新赋值
visualization = filter(visualization, !(celltype == "0" | celltype == "Unknown"))
visualization = filter(visualization, celltype == "NfkB_Pathway_related")
head(visualization,100)

ggplot(data = visualization,mapping = aes(x=column, y=row, color = celltype))+
  geom_point(size=2,colour="#2a00ff")+
  scale_y_reverse()+
  theme_void()+
  theme(panel.background = element_rect(fill = 'black', color = 'black'))
ggsave(last_plot(),filename = "NfkB_Pathway_related.png", device = "png", height = 10,width = 10,dpi = 300)

### RTK_Pathway_related
visualization = data.frame(centroids,celltype=0) #新建celltype列
visualization[rownames(CellType),]$celltype <- CellType$CellType #特定行重新赋值
visualization = filter(visualization, !(celltype == "0" | celltype == "Unknown"))
visualization = filter(visualization, celltype == "RTK_Pathway_related")
head(visualization,100)

ggplot(data = visualization,mapping = aes(x=column, y=row, color = celltype))+
  geom_point(size=2,colour="#b400ff")+
  scale_y_reverse()+
  theme_void()+
  theme(panel.background = element_rect(fill = 'black', color = 'black'))
ggsave(last_plot(),filename = "RTK_Pathway_related.png", device = "png", height = 10,width = 10,dpi = 300)

### Wnt_Pathway_related
visualization = data.frame(centroids,celltype=0) #新建celltype列
visualization[rownames(CellType),]$celltype <- CellType$CellType #特定行重新赋值
visualization = filter(visualization, !(celltype == "0" | celltype == "Unknown"))
visualization = filter(visualization, celltype == "Wnt_Pathway_related")
head(visualization,100)

ggplot(data = visualization,mapping = aes(x=column, y=row, color = celltype))+
  geom_point(size=2,colour="#ff00bf")+
  scale_y_reverse()+
  theme_void()+
  theme(panel.background = element_rect(fill = 'black', color = 'black'))
ggsave(last_plot(),filename = "Wnt_Pathway_related.png", device = "png", height = 10,width = 10,dpi = 300)

##############QC小提琴图########################################################
rm(list = ls())
options(stringsAsFactors = F)
setwd("D:/11.CellSeg_ColonV2_0613/3.cellpose/2.概率模型_0722/")

mx <- read.csv("Matrix.csv",header = T,check.names = F)
rownames(mx) <- mx$gene
mx <- mx[,-c(1:2)] #去掉基因列、“0”列

meta <- read.csv("CellMeta.csv",header = T,check.names = F)
head(meta)
index <- as.character(meta$cell)
mx <- mx[,index] #矩阵只保留去除了“Unknown”的细胞

unique_gene <- data.frame()
for (i in 1:ncol(mx)){
  unique_gene[i,1] <- sum(mx[,i] > 0)
}
colnames(unique_gene)[1] <- 'gene'

result <- data.frame(
  cell = meta$cell,
  celltype = meta$celltype,
  count = apply(mx,2,sum),
  gene = unique_gene
)
head(result)

### QC1:count数目
p <- ggplot(result,aes(x=celltype,y=count))+
  geom_violin(aes(fill = factor(celltype)),adjust=2)+
  geom_boxplot(width = 0.1, fill = 'black')+
  stat_summary(fun = median, geom = 'point', fill = 'white', shape = 21, size = 2.5)+
  theme_classic()+
  scale_fill_manual(values = c(Malignancy = "#ff0029",
                               Epithelial = '#ff6000',
                               Normal = '#ffea00',
                               Stromal = '#8aff00',
                               B_Cell = '#00ff00',
                               T_Cell = '#00ff89',
                               Mast_Cell = '#00ecff',
                               Myeloids = '#0061ff',
                               NfkB_Pathway_related = '#2a00ff',
                               RTK_Pathway_related = '#b400ff',
                               Wnt_Pathway_related = '#ff00bf'))+ #颜色对应关系；图例顺序关系
  labs(fill='Cell Type')+
  scale_x_discrete(limits = c("Malignancy",
                              "Epithelial",
                              "Normal","Stromal",
                              "B_Cell",
                              "T_Cell",
                              "Mast_Cell",
                              "Myeloids",
                              "NfkB_Pathway_related",
                              "RTK_Pathway_related",
                              "Wnt_Pathway_related"))+ #X轴顺序
  theme(axis.text.x = element_text(face = 'bold.italic',size = rel(2),hjust = 0.5,vjust = 0.5,angle = 45),
        axis.text.y = element_text(face = 'bold.italic',size = rel(2),hjust = 0.5,vjust = 0.5),
        axis.title.x = element_text(face = 'bold.italic',size = rel(2),hjust = 0.5,vjust = 0.5,margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(face = 'bold.italic',size = rel(2),hjust = 0.5,vjust = 0.5,margin = margin(t = 0, r = 20, b = 0, l = 0)),
        legend.text = element_text(face = 'bold.italic',size = rel(1),hjust = 0.5,vjust = 0.5),
        legend.title = element_text(face = 'bold.italic',size = rel(1.2),hjust = 0.5,vjust = 0.5))+
  xlab("Cell Type")+
  ylab("Number of total reads per cell")
p
ggsave(p,filename = "QC1.png",width = 10,height = 8,device = 'png',dpi = 300)

### QC2:gene数目
p <- ggplot(result,aes(x=celltype,y=gene))+
  geom_violin(aes(fill = factor(celltype)),adjust=2)+
  geom_boxplot(width = 0.1, fill = 'black')+
  stat_summary(fun = median, geom = 'point', fill = 'white', shape = 21, size = 2.5)+
  theme_classic()+
  scale_fill_manual(values = c(Malignancy = "#ff0029",
                               Epithelial = '#ff6000',
                               Normal = '#ffea00',
                               Stromal = '#8aff00',
                               B_Cell = '#00ff00',
                               T_Cell = '#00ff89',
                               Mast_Cell = '#00ecff',
                               Myeloids = '#0061ff',
                               NfkB_Pathway_related = '#2a00ff',
                               RTK_Pathway_related = '#b400ff',
                               Wnt_Pathway_related = '#ff00bf'))+ #颜色对应关系；图例顺序关系
  labs(fill='Cell Type')+
  scale_x_discrete(limits = c("Malignancy",
                              "Epithelial",
                              "Normal","Stromal",
                              "B_Cell",
                              "T_Cell",
                              "Mast_Cell",
                              "Myeloids",
                              "NfkB_Pathway_related",
                              "RTK_Pathway_related",
                              "Wnt_Pathway_related"))+ #X轴顺序
  theme(axis.text.x = element_text(face = 'bold.italic',size = rel(2),hjust = 0.5,vjust = 0.5,angle = 45),
        axis.text.y = element_text(face = 'bold.italic',size = rel(2),hjust = 0.5,vjust = 0.5),
        axis.title.x = element_text(face = 'bold.italic',size = rel(2),hjust = 0.5,vjust = 0.5,margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(face = 'bold.italic',size = rel(2),hjust = 0.5,vjust = 0.5,margin = margin(t = 0, r = 20, b = 0, l = 0)),
        legend.text = element_text(face = 'bold.italic',size = rel(1),hjust = 0.5,vjust = 0.5),
        legend.title = element_text(face = 'bold.italic',size = rel(1.2),hjust = 0.5,vjust = 0.5))+
  xlab("Cell Type")+
  ylab("Number of unique genes per cell")
p
ggsave(p,filename = "QC2.png",width = 10,height = 8,device = 'png',dpi = 300)

##############基因热图###################
rm(list = ls())

mx_path <- file.path("Matrix.csv")
mx <- read.csv(mx_path,header = T,check.names = F)
rownames(mx) <- mx$gene
mx <- mx[,3:ncol(mx)] #第一列为gene，第二列为'0'，所以删掉前两列 

mx <- as.data.frame(apply(mx,2,prop.table))

### CellMeta.csv
meta <- read.csv("CellMeta.csv",header = T,check.names = F)
head(meta)

### 热图矩阵
heatmap <- function(type){
  index <- filter(meta, celltype == type)
  index <- as.character(index$cell)
  mx_filtered <- mx[,index]
  
  out <- apply(mx_filtered, 1, mean)
  return(out)
}

Malignancy <- heatmap(type = "Malignancy")
Epithelial <- heatmap(type = "Epithelial")
Normal <- heatmap(type = "Normal")
Stromal <- heatmap(type = "Stromal")
B_Cell <- heatmap(type = "B_Cell")
T_Cell <- heatmap(type = "T_Cell")
Mast_Cell <- heatmap(type = "Mast_Cell")
Myeloids <- heatmap(type = "Myeloids")
NfkB_Pathway_related <- heatmap(type = "NfkB_Pathway_related")
RTK_Pathway_related <- heatmap(type = "RTK_Pathway_related")
Wnt_Pathway_related <- heatmap(type = "Wnt_Pathway_related")

result <- data.frame(Malignancy,Epithelial,Normal,Stromal,B_Cell,T_Cell,Mast_Cell,Myeloids,NfkB_Pathway_related,RTK_Pathway_related,Wnt_Pathway_related)

### 调整细胞类型顺序(为了使热图呈阶梯状)
TypeOrder <- c("T_Cell",
               "Mast_Cell",
               "Wnt_Pathway_related",
               "B_Cell",
               "NfkB_Pathway_related",
               "Normal",
               "RTK_Pathway_related",
               "Stromal",
               "Myeloids",
               "Malignancy",
               "Epithelial")
### 调整基因顺序
GeneOrder <- c("CEACAM5","CD24", #2
               "PPIB","REG4","LDHA","S100A4","SOX4","CD63","HSPB1","PIGR", #8
               "TAGLN",
               "ACTA2","CAV1","ENG","TIMP2","DES","TNC","THY1","ITGB1","TIMP1", #9
               "CD79A","CD81", #2
               "CD3D","CD8A","CD7","CD3E","CD3G","CD8B","CD56", #7
               "CPA3","MS4A2","CTSG", #3
               "SPP1","CD68","CD80","S100A8", #4
               "NFKBIA",
               "PDGFRB",
               "OLFM4","WNT5A") #2
result <- result[GeneOrder,TypeOrder]
head(result)

###
library(pheatmap)
### 基因集注释
annotation_row = data.frame(
  Type = factor(rep(c("Malignancy",
                      "Epithelial",
                      "Normal","Stromal",
                      "B_Cell",
                      "T_Cell",
                      "Mast_Cell",
                      "Myeloids",
                      "NfkB_Pathway_related",
                      "RTK_Pathway_related",
                      "Wnt_Pathway_related"), c(2,8,1,9,2,7,3,4,1,1,2)))
)
rownames(annotation_row) <- rownames(result)
### 标签颜色指定
ann_colors = list(
  Type = c(Malignancy = "#ff0029",
           Epithelial = '#ff6000',
           Normal = '#ffea00',
           Stromal = '#8aff00',
           B_Cell = '#00ff00',
           T_Cell = '#00ff89',
           Mast_Cell = '#00ecff',
           Myeloids = '#0061ff',
           NfkB_Pathway_related = '#2a00ff',
           RTK_Pathway_related = '#b400ff',
           Wnt_Pathway_related = '#ff00bf')
)

pheatmap(result, 
         cluster_rows = T, 
         cluster_cols = F, 
         annotation_row = annotation_row,
         color = colorRampPalette(c("black","white","firebrick3"))(50), 
         annotation_colors = ann_colors,
         scale = "row",fontsize = 10,angle_col = 45,width = 10,height = 10,filename = "heatmap.png")



