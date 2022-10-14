# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 16:26:43 2022

@author: Administrator
"""

import os
os.getcwd() #获得当前工作路径
os.chdir("D:/11.CellSeg_ColonV2_0613/3.cellpose") #改变工作路径
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from skimage.segmentation import expand_labels
import tifffile
import skimage

### 示例教程
### https://stackoverflow.com/questions/68341726/how-to-test-if-coordinates-are-inside-a-label-or-mask

### 步骤一：导入cellpose训练好的mask array
dat = np.load('D:/11.CellSeg_ColonV2_0613/3.cellpose/DAPI_seg.npy', allow_pickle=True).item()
masks = dat['masks']
print('The image has %d cells' % np.unique(masks).max())
# The image has 82,966 cells

### 步骤二：找每个mask的质心，生成centroid.csv
from scipy import ndimage

result = ndimage.center_of_mass(masks, masks, np.unique(masks)) # 原始数组；标记数组；需要生成质心的细胞编号
resultList = list(result) #为了下一步转换成data frame
df = pd.DataFrame(
    resultList, columns=['row', 'column'])
df.drop(0, inplace=True, axis=0) #删除第一行背景缺失值
print(df)
df.to_csv("Centroid.csv")

### 步骤三：mask dilation
### dilation前
plt.figure(figsize=(20,20))
plt.style.use('dark_background')
plt.imshow(masks)
plt.savefig("masks.png",dpi=300)
### dilation后
expanded = expand_labels(masks, distance=20) #mask外扩20像素
expanded = expanded.astype('uint32')

plt.figure(figsize=(20,20))
plt.style.use('dark_background')
plt.imshow(expanded)
plt.savefig("masks_dilation.png",dpi=300)
### dilation masks叠加dapi
filepath="D:/11.CellSeg_ColonV2_0613/1.raw/DAPI.tif"
image = tifffile.imread(filepath) #dapi图
rgb_label_image = skimage.color.label2rgb(expanded, bg_label=0) #masks数值转rgb

plt.style.use('dark_background')
plt.figure(figsize=(60, 60), dpi=300)
background = plt.imshow(image)
imgplot = plt.imshow(rgb_label_image,alpha=0.6)
plt.axis('off')
plt.savefig("mask_rgb.tiff", dpi=400, bbox_inches="tight")

'''
注意事项：
plt.imshow(masks)时，Y轴是反过来的，X轴是正常的

mask array索引说明：
mask[300:600,50:125] 300-600行 50-125列

示例数据（方便解释）：
coords = np.random.rand(40, 2) *1024
mask = np.zeros((1024,1024))
mask[300:600,50:125] = 1
mask[700:800,400:650] = 2

plt.imshow(mask)
plt.scatter(coords[:,0],coords[:,1],color='red')

'''

### ISS spots
spots = pd.read_csv('D:/11.CellSeg_ColonV2_0613/Barcode_data.csv')

plt.style.use('dark_background')
plt.figure(figsize=(20, 20), dpi=300)
plt.imshow(expanded)
plt.scatter(spots.y,spots.x,color='red',s=3,alpha=0.2)
plt.savefig("spot_overlay.png",dpi=300)

'''
mask array里的row=散点图里的y
mask array的col=散点图里的x
所以画散点图时需转换为plt.scatter(spots.col,spots.row)

示例数据（方便解释）：
mask = np.zeros((1024,1024))
mask[300:600,50:125] = 1
mask[700:800,400:650] = 2

plt.imshow(mask)
plt.scatter([131,211,355,931,433],[222,313,555,561,882],color='red')
'''

### assign spots to cells
coords = spots.iloc[:,1:3] #提取信号点坐标
coords = coords[['x','y']] #列重命名
coords_int = np.round(coords).astype(int) #取整
coords_int = coords_int.to_numpy() #数据框转换为数组
values_at_coords = expanded[tuple(coords_int.T)].astype(int) #每个信号点的细胞编号
points_per_value = np.bincount(values_at_coords) #每个细胞的全基因总count数

### add a column named 'cell' to spots data frame
spots = spots.assign(cell=values_at_coords)
spots.to_csv("spots_summary.csv")
matrix = pd.crosstab(spots.gene.tolist(), spots.cell, rownames = ['gene'], colnames = ['cell']) #pd.crosstab与R里的table等价
matrix.to_csv("Matrix.csv")

### 学习笔记
'''
np.bincount学习：
np.bincount([1,1,1,7,2,3,4,5,5,5])

掩码数组学习：
coords = np.random.rand(40, 2) *1024
mask = np.zeros((1024,1024))
mask[300:600,50:125] = 1 在掩码数组的特定区域添加编号
mask[700:800,400:650] = 2
mask[723,591] 提取掩码数组里的某个值

'''

#############################
### 精细分析：scanpy独立流程
#############################
import anndata as ad

matrix.drop(0, inplace=True, axis=1) #删除第一列背景信号点

### Anndata三要素：原始矩阵
mx = matrix.T.to_numpy() #原始表达矩阵，行为细胞，列为基因
### Anndata三要素：obs
celltype = pd.read_csv("D:/11.CellSeg_ColonV2_0613/3.cellpose/2.概率模型_0722/Giotto/Giotto_Filter0.csv") #做好的细胞类型注释
celltype.index = celltype['cell'] #修改index为细胞编号
observation = pd.DataFrame({'cell':matrix.T.index,
                            'row':df.iloc[matrix.T.index-1,0],
                            'column':df.iloc[matrix.T.index-1,1],
                            'celltype':celltype['celltype']}) #细胞编号
### Anndata三要素：var
variable = pd.DataFrame(matrix.T.columns) #基因名称

adata = ad.AnnData(mx, obs=observation,var=variable) 
adata.X #原始矩阵
adata.obs #细胞
adata.var #基因

adata.write("Colon.h5ad") #AnnData











