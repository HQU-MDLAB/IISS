library(dplyr)

rm(list = ls())
options(stringsAsFactors = F)
setwd("D:/9.CellSeg_MouseBrainV3_0519/6.cellpose/")

##########################函数主体##############################################
### 提取打分分值相同的细胞编号
### 过滤出同时符合以下条件的细胞编号：(1)分值大于0的个数与分值最大值个数一样；(2)分值最大值不为1

Find_Multiple_Maximum_Index <- function(result){
  index <- c()
  for (c in 1:ncol(result)){
    cell <- result[,c] #对每一列进行操作，每一列为一个细胞
    GreaterZeroNum <- length(cell[which(cell > 0)]) #打分大于0的类型个数
    MaxValueNum <- length(which(cell == max(cell))) #打分最大值出现的个数
    if (MaxValueNum == GreaterZeroNum) { 
      if (max(cell) != 1) {index=c(index,colnames(result[c]))} #index：Unknown细胞类型对应的细胞编号
    }
  }
  return(index)
}

### 细胞类型打分
Scoring <- function(mx_path){
  
  ### 导入矩阵
  mx <- read.csv(mx_path,header = T,check.names = F)
  rownames(mx) <- mx$gene
  mx <- mx[,3:ncol(mx)] #第一列为gene，第二列为'0'，所以删掉前两列 
  
  ### 第一步：定义细胞类型marker
  Inhibitory_Neurons <- c("Gad2","Slc32a1","Crhbp","Pthlh","Cnr1","Vip","Cpne5","Crh","Kcnip2")
  Excitatory_Neurons <- c("Lamp5","Rorb","Syt6","Tbr1","Slc17a7")
  Astrocytes <- c("Gfap","Serpinf1","Mfge8","Aldoc")
  Oligodendrocytes <- c("Plp1","Pdgfra","Cemip2","Itpr2","Ctps","Bmp4","Anln","Sox10","Mbp")
  Immune_Cells <- c("Hexb","Mrc1")
  Vasculature <- c("Flt1","Apln","Acta2","Vtn")
  Ependymal <- c("Ttr","Foxj1","Cd24a")
  
  CellTypeList <- list(Inhibitory_Neurons,Excitatory_Neurons,Astrocytes,
                       Oligodendrocytes,Immune_Cells,Vasculature,Ependymal)
  
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
  rownames(result) <- c('Inhibitory_Neurons','Excitatory_Neurons','Astrocytes',
                        'Oligodendrocytes','Immune_Cells','Vasculature','Ependymal')
  
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
centroids <- read.csv("Centroid.csv",header = T) 
centroids <- centroids[,-1] #删除第一列索引值
head(centroids)

### visualization：用于细胞类型空间位置可视化的data frame
visualization <-  data.frame(centroids,cell = rownames(centroids),celltype=0) #新建celltype列
visualization[rownames(CellType),]$celltype <- CellType$CellType #特定行重新赋值
visualization <-  filter(visualization, !(celltype == "0" | celltype == "Unknown")) #去掉两种细胞类型
visualization <- visualization[,c('cell','row','column','celltype')] #列重排
head(visualization,100)
write.csv(visualization,"CellMeta.csv",row.names = F) #row.names = F，删除行索引所在的那一列
### 0：没有信号点的细胞
### Unknown：打分不明确的细胞

### 总图：七种细胞类型
ggplot(data = visualization,mapping = aes(x=column, y=row, color = celltype))+
  geom_point(size=0.7)+
  scale_y_reverse()+
  theme_void()+
  theme(panel.background = element_rect(fill = 'black', color = 'black'),
        legend.position = "none",
        legend.text = element_text(face = 'italic', color = "white", size = 10),
        legend.title = element_text(face = 'italic', color = "white", size = 13))+
  scale_colour_manual(values = c("#ff0029","#ffbc00",'#5cff00','#00ff89','#008fff','cyan','#ff00bf'),
                      limits = c("Inhibitory_Neurons","Excitatory_Neurons","Astrocytes","Oligodendrocytes","Immune_Cells","Vasculature","Ependymal"))
ggsave(last_plot(),filename = "CellTypeMap.png", device = "png", height = 8,width = 10,dpi = 300)

ggplot(data = visualization,mapping = aes(x=column, y=row, color = celltype))+
  geom_point(size=0.5)+
  scale_y_reverse()+
  theme_void()+
  facet_wrap( ~ celltype)+
  theme(strip.text = element_text(face = 'bold', color = "black", size = rel(1.5)),
        strip.background = element_rect(fill = 'lightblue', color = 'black'))+
  theme(panel.background = element_rect(fill = 'black', color = 'black'),
        legend.position = "none",
        legend.text = element_text(face = 'italic', color = "white", size = 10),
        legend.title = element_text(face = 'italic', color = "white", size = 13))+
  scale_colour_manual(values = c("#ff0029","#ffbc00",'#5cff00','#00ff89','#008fff','cyan','#ff00bf'),
                      limits = c("Inhibitory_Neurons","Excitatory_Neurons","Astrocytes","Oligodendrocytes","Immune_Cells","Vasculature","Ependymal"))
ggsave(last_plot(),filename = "CellTypeMap_facet.png", device = "png", height = 18,width = 20,dpi = 300)


### 单类型图
### Inhib

visualization = data.frame(centroids,celltype=0) #新建celltype列
visualization[rownames(CellType),]$celltype <- CellType$CellType #特定行重新赋值
visualization = filter(visualization, !(celltype == "0" | celltype == "Unknown"))
visualization = filter(visualization, celltype == "Inhibitory_Neurons")
head(visualization,100)

ggplot(data = visualization,mapping = aes(x=column, y=row, color = celltype))+
  geom_point(size=1,colour="#ff0029")+
  scale_y_reverse()+
  theme_void()+
  theme(panel.background = element_rect(fill = 'black', color = 'black'))
ggsave(last_plot(),filename = "Inhib.png", device = "png", height = 8,width = 10,dpi = 300)

### Exc
visualization = data.frame(centroids,celltype=0) #新建celltype列
visualization[rownames(CellType),]$celltype <- CellType$CellType #特定行重新赋值
visualization = filter(visualization, !(celltype == "0" | celltype == "Unknown"))
visualization = filter(visualization, celltype == "Excitatory_Neurons")
head(visualization,100)

ggplot(data = visualization,mapping = aes(x=column, y=row, color = celltype))+
  geom_point(size=1,colour="#ffbc00")+
  scale_y_reverse()+
  theme_void()+
  theme(panel.background = element_rect(fill = 'black', color = 'black'))
ggsave(last_plot(),filename = "Exc.png", device = "png", height = 8,width = 10,dpi = 300)

### Astro
visualization = data.frame(centroids,celltype=0) #新建celltype列
visualization[rownames(CellType),]$celltype <- CellType$CellType #特定行重新赋值
visualization = filter(visualization, !(celltype == "0" | celltype == "Unknown"))
visualization = filter(visualization, celltype == "Astrocytes")
head(visualization,100)

ggplot(data = visualization,mapping = aes(x=column, y=row, color = celltype))+
  geom_point(size=1,colour="#5cff00")+
  scale_y_reverse()+
  theme_void()+
  theme(panel.background = element_rect(fill = 'black', color = 'black'))
ggsave(last_plot(),filename = "Astro.png", device = "png", height = 8,width = 10,dpi = 300)

### Oligo
visualization = data.frame(centroids,celltype=0) #新建celltype列
visualization[rownames(CellType),]$celltype <- CellType$CellType #特定行重新赋值
visualization = filter(visualization, !(celltype == "0" | celltype == "Unknown"))
visualization = filter(visualization, celltype == "Oligodendrocytes")
head(visualization,100)

ggplot(data = visualization,mapping = aes(x=column, y=row, color = celltype))+
  geom_point(size=1,colour="#00ff89")+
  scale_y_reverse()+
  theme_void()+
  theme(panel.background = element_rect(fill = 'black', color = 'black'))
ggsave(last_plot(),filename = "Oligo.png", device = "png", height = 8,width = 10,dpi = 300)

### Immune
visualization = data.frame(centroids,celltype=0) #新建celltype列
visualization[rownames(CellType),]$celltype <- CellType$CellType #特定行重新赋值
visualization = filter(visualization, !(celltype == "0" | celltype == "Unknown"))
visualization = filter(visualization, celltype == "Immune_Cells")
head(visualization,100)

ggplot(data = visualization,mapping = aes(x=column, y=row, color = celltype))+
  geom_point(size=1,colour="#008fff")+
  scale_y_reverse()+
  theme_void()+
  theme(panel.background = element_rect(fill = 'black', color = 'black'))
ggsave(last_plot(),filename = "Immune.png", device = "png", height = 8,width = 10,dpi = 300)

### Vasculature
visualization = data.frame(centroids,celltype=0) #新建celltype列
visualization[rownames(CellType),]$celltype <- CellType$CellType #特定行重新赋值
visualization = filter(visualization, !(celltype == "0" | celltype == "Unknown"))
visualization = filter(visualization, celltype == "Vasculature")
head(visualization,100)

ggplot(data = visualization,mapping = aes(x=column, y=row, color = celltype))+
  geom_point(size=1,colour="cyan")+
  scale_y_reverse()+
  theme_void()+
  theme(panel.background = element_rect(fill = 'black', color = 'black'))
ggsave(last_plot(),filename = "Vas.png", device = "png", height = 8,width = 10,dpi = 300)

### Ependymal
visualization = data.frame(centroids,celltype=0) #新建celltype列
visualization[rownames(CellType),]$celltype <- CellType$CellType #特定行重新赋值
visualization = filter(visualization, !(celltype == "0" | celltype == "Unknown"))
visualization = filter(visualization, celltype == "Ependymal")
head(visualization,100)

ggplot(data = visualization,mapping = aes(x=column, y=row, color = celltype))+
  geom_point(size=1,colour="#ff00bf")+
  scale_y_reverse()+
  theme_void()+
  theme(panel.background = element_rect(fill = 'black', color = 'black'))
ggsave(last_plot(),filename = "Epen.png", device = "png", height = 8,width = 10,dpi = 300)



##############QC小提琴图########################################################
library(dplyr)
library(ggplot2)

rm(list = ls())
options(stringsAsFactors = F)
setwd("D:/9.CellSeg_MouseBrainV3_0519/6.cellpose/")

mx <- read.csv("Matrix.csv",header = T,check.names = F)
rownames(mx) <- mx$gene
mx <- mx[,-c(1:2)]

meta <- read.csv("CellMeta.csv",header = T,check.names = F)
head(meta)
index <- as.character(meta$cell)
mx <- mx[,index]

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
  scale_fill_manual(values=c("#ff0029","#ffbc00",'#5cff00','#00ff89','#008fff','cyan','#ff00bf'),
                    limits = c("Inhibitory_Neurons","Excitatory_Neurons","Astrocytes","Oligodendrocytes","Immune_Cells","Vasculature","Ependymal"),
                    breaks = c("Inhibitory_Neurons","Excitatory_Neurons","Astrocytes","Oligodendrocytes","Immune_Cells","Vasculature","Ependymal"),
                    labels = c("Inhib","Exc","Astro","Oligo","Immune","Vas","Epen"))+
  labs(fill='Cell Type')+
  scale_x_discrete(limits = c("Inhibitory_Neurons","Excitatory_Neurons","Astrocytes","Oligodendrocytes","Immune_Cells","Vasculature","Ependymal"),
                   breaks = c("Inhibitory_Neurons","Excitatory_Neurons","Astrocytes","Oligodendrocytes","Immune_Cells","Vasculature","Ependymal"),
                   labels = c("Inhib","Exc","Astro","Oligo","Immune","Vas","Epen"))+
  theme(axis.text.x = element_text(face = 'bold.italic',size = rel(2),hjust = 0.5,vjust = 0.5),
        axis.text.y = element_text(face = 'bold.italic',size = rel(2),hjust = 0.5,vjust = 0.5),
        axis.title.x = element_text(face = 'bold.italic',size = rel(2),hjust = 0.5,vjust = 0.5,margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(face = 'bold.italic',size = rel(2),hjust = 0.5,vjust = 0.5,margin = margin(t = 0, r = 20, b = 0, l = 0)),
        legend.text = element_text(face = 'bold.italic',size = rel(2),hjust = 0.5,vjust = 0.5),
        legend.title = element_text(face = 'bold.italic',size = rel(1.8),hjust = 0.5,vjust = 0.5))+
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
  scale_fill_manual(values=c("#ff0029","#ffbc00",'#5cff00','#00ff89','#008fff','cyan','#ff00bf'),
                    limits = c("Inhibitory_Neurons","Excitatory_Neurons","Astrocytes","Oligodendrocytes","Immune_Cells","Vasculature","Ependymal"),
                    breaks = c("Inhibitory_Neurons","Excitatory_Neurons","Astrocytes","Oligodendrocytes","Immune_Cells","Vasculature","Ependymal"),
                    labels = c("Inhib","Exc","Astro","Oligo","Immune","Vas","Epen"))+
  labs(fill='Cell Type')+
  scale_x_discrete(limits = c("Inhibitory_Neurons","Excitatory_Neurons","Astrocytes","Oligodendrocytes","Immune_Cells","Vasculature","Ependymal"),
                   breaks = c("Inhibitory_Neurons","Excitatory_Neurons","Astrocytes","Oligodendrocytes","Immune_Cells","Vasculature","Ependymal"),
                   labels = c("Inhib","Exc","Astro","Oligo","Immune","Vas","Epen"))+
  theme(axis.text.x = element_text(face = 'bold.italic',size = rel(2),hjust = 0.5,vjust = 0.5),
        axis.text.y = element_text(face = 'bold.italic',size = rel(2),hjust = 0.5,vjust = 0.5),
        axis.title.x = element_text(face = 'bold.italic',size = rel(2),hjust = 0.5,vjust = 0.5,margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(face = 'bold.italic',size = rel(2),hjust = 0.5,vjust = 0.5,margin = margin(t = 0, r = 20, b = 0, l = 0)),
        legend.text = element_text(face = 'bold.italic',size = rel(2),hjust = 0.5,vjust = 0.5),
        legend.title = element_text(face = 'bold.italic',size = rel(1.8),hjust = 0.5,vjust = 0.5))+
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

inhib <- heatmap(type = "Inhibitory_Neurons")
exc <- heatmap(type = "Excitatory_Neurons")
astro <- heatmap(type = "Astrocytes")
oligo <- heatmap(type = "Oligodendrocytes")
immune <- heatmap(type = "Immune_Cells")
vas <- heatmap(type = "Vasculature")
epen <- heatmap(type = "Ependymal")

result <- data.frame(
  Inhibitory_Neurons <- inhib,
  Excitatory_Neurons <- exc,
  Astrocytes <- astro,
  Oligodendrocytes <- oligo,
  Immune_Cells <- immune,
  Vasculature <- vas,
  Ependymal <- epen
)
### 修改列名
colnames(result) <- c('Inhib','Exc','Astro','Oligo','Immune','Vas','Epen')
TypeOrder <- c("Immune","Oligo","Vas","Astro","Inhib","Exc","Epen")
### 调整基因顺序
GeneOrder <- c("Gad2","Slc32a1","Crhbp","Pthlh","Cnr1","Vip","Cpne5","Crh","Kcnip2", #9
               "Lamp5","Rorb","Syt6","Tbr1","Slc17a7", #5
               "Gfap","Serpinf1","Mfge8","Aldoc", #4
               "Plp1","Pdgfra","Cemip2","Itpr2","Ctps","Bmp4","Anln","Sox10","Mbp", #9
               "Hexb","Mrc1", #2
               "Flt1","Apln","Acta2","Vtn", #4
               "Ttr","Foxj1","Cd24a") #3
result <- result[GeneOrder,TypeOrder]
head(result)

###
library(pheatmap)
### 基因集注释
annotation_row = data.frame(
  Type = factor(rep(c("Inhib", "Exc", "Astro","Oligo","Immune","Vas","Epen"), c(9,5,4,9,2,4,3)))
)
rownames(annotation_row) <- rownames(result)
### 标签颜色指定
ann_colors = list(
  Type = c(Inhib = "#ff0029",
           Exc = "#ffbc00",
           Astro = "#5cff00",
           Oligo = "#00ff89",
           Immune = "#008fff",
           Vas = "cyan",
           Epen = "#ff00bf")
)

pheatmap(result, 
         cluster_rows = T, 
         cluster_cols = F, 
         annotation_row = annotation_row,
         color = colorRampPalette(c("black","white","firebrick3"))(50), 
         annotation_colors = ann_colors,
         scale = "row",fontsize = 11,angle_col = 0,width = 10,height = 10, filename = "heatmap.png")
