library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(patchwork)
library(harmony)
library(rliger)
library(reshape2)
library(RColorBrewer)
library(dplyr)
library(readr)
library(tidyverse)
library(Matrix)
library(Seurat)
library(ggplot2)
library(cowplot)
library(celldex)
library(SingleR)
library(data.table)
library(scGate)
library(ggpubr)

umi<-read.delim("GSE165897_UMIcounts_HGSOC.tsv",row.names=1)
OV<- CreateSeuratObject(counts = umi, project = "OV")
OV[["percent.mt"]] <- PercentageFeatureSet(OV, pattern = "^MT-")
VlnPlot(OV, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


# Removal of poor quality cells
 OV<- subset(OV, subset = percent.mt < 15 & nFeature_RNA > 500 & nFeature_RNA < 7000)

# meta.data
meta<-fread("meta.data_update.txt")
OV@meta.data[["Treatment"]]=paste0(meta$Treatment)
OV@meta.data[["group"]]=paste0(meta$group)

OV[["RNA"]] <- split(OV[["RNA"]], f = OV$group)
OV<- NormalizeData(OV)
OV<- FindVariableFeatures(OV)
OV<- ScaleData(OV)
OV<- RunPCA(OV)

# Visualization the results of a standard analysis without integration

OV<- FindNeighbors(OV, dims = 1:30, reduction = "pca")
OV<- FindClusters(OV, resolution = 2, cluster.name = "unintegrated_clusters")
OV<- RunUMAP(OV, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")


# Integration
options(future.globals.maxSize = 8000 * 1024^2)
OV<- IntegrateLayers(
  object = OV, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca",
  verbose = FALSE
)

# Visualization and cluster the datasets
# RPCA
OV<- FindNeighbors(OV, reduction = "integrated.rpca", dims = 1:30)
OV<- FindClusters(OV, resolution = 2, cluster.name = "rpca_clusters")

OV<- RunUMAP(OV, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca")
p1 <- DimPlot(
  OV,
  reduction = "umap.rpca",
  group.by = c("group","rpca_clusters"),
  combine = FALSE, label.size = 4,label=F
)

OV<- JoinLayers(OV)

# Cluster annotation using CellDex
sample<- as.SingleCellExperiment(OV)
IGD<- BlueprintEncodeData()
pred.IGD<- SingleR(test = sample, ref =IGD, assay.type.test=1,labels = IGD$label.fine)
OV <- AddMetaData( object=OV,metadata=make.names(sub( " [(].*$","",pred.IGD$pruned.labels )),col.name="BlueprintEncodeData" )
DimPlot(OV,reduction = "umap.rpca",group.by="BlueprintEncodeData",label=T,repel=T,split.by="Treatment")+NoLegend()

#CD8T cells subset

library("stringr")
OV@meta.data=OV@meta.data%>%mutate(T.cells=case_when(rpca_clusters==14~"CD8.pos",
rpca_clusters==3~"CD8.pos",
rpca_clusters==0~"CD8.pos",
rpca_clusters==5~"CD8.pos",
rpca_clusters==23~"CD8.pos"))
CD8T=subset(OV,subset=T.cells=="CD8.pos")




