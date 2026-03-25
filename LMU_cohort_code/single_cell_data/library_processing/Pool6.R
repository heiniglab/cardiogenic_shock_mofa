---
title: "fold_change"
author: "Visha"
date: "13/10/2023"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```
source("/mnt/raidtmp/Alejandro/functions.R")
## R Markdown

### implementing SoupX to reduce ambient RNA content and reading in input matrices
```{r reading_in_data, echo=TRUE}
#.rs.restartR()
library(SoupX)
library(dplyr)
library(patchwork)
library(sctransform)
library(ggplot2)
library(dittoSeq)
library(RColorBrewer)
library(data.table)
library(magrittr)
library(readxl)
library(reshape2)
library(enrichR)
library(Seurat)
library(clusterProfiler)
library(stringr)
library(EnhancedVolcano)
library(cowplot)


foldername = dirname(Sys.glob("/mnt/raidtmp/Alejandro/shock_data/home/jovyan/bhavishya-cardiogenic-shock-volume/Pool6/outs/raw_feature_bc_matrix/*.mtx.gz"))
h5file = Read10X(foldername,unique.features = TRUE)
samplename = "Pool6"
allfiles.raw = list()
allABs.raw = list()

allfiles.raw[[samplename]] = h5file$`Gene Expression`
allABs.raw[[samplename]] = h5file$`Antibody Capture`


joint.bcs <- intersect(colnames(allfiles.raw$Pool6), colnames(allABs.raw$Pool6))
allfiles.raw$Pool6 <- allfiles.raw$Pool6[, joint.bcs]
allABs.raw$Pool6 <- as.matrix(allABs.raw$Pool6[, joint.bcs])

print("conditions")
rownames(allABs.raw$Pool6)

```

### filtering data and removing ribosomal genes
```{r filtering, echo=TRUE}
Pool6 <- CreateSeuratObject(counts = allfiles.raw$Pool6, project = "Pool6", min.cells = 3, min.features = 200)
str(Pool6)

### removing ribosomal genes from the data
#rps_rownames = rownames(Pool6[grep('^Rps', rownames(Pool6))])
#rpl_rownames = rownames(Pool6[grep('^Rpl', rownames(Pool6))])
#Pool6 <- Pool6[!(row.names(Pool6) %in% rps_rownames),]
#Pool6 <- Pool6[!(row.names(Pool6) %in% rpl_rownames),]
Pool6[["percent.mt"]] <- PercentageFeatureSet(Pool6, pattern = "^MT-")

## seurat object before filtering
print("seurat object before filtering")
p = VlnPlot(Pool6, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
save_plot(p, "/ictstr01/home/icb/bhavishya.nelakuditi/pre_processing/seurat_processig/Pool6/before_filtered_qc", fig.width=10, fig.height=6)
p=FeatureScatter(Pool6, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = 'lm')
save_plot(p, "/ictstr01/home/icb/bhavishya.nelakuditi/pre_processing/seurat_processig/Pool6/before_filtered_qc2", fig.width=10, fig.height=6)

### filtering data
Pool6 <- subset(Pool6, subset = nFeature_RNA > 100 & nFeature_RNA < 8000 & nCount_RNA > 100 & nCount_RNA < 45000)
Pool6 <- subset(Pool6, subset = percent.mt < 15)

## seurat object after filtering
print("seurat object after filtering")
p = VlnPlot(Pool6, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
save_plot(p, "/ictstr01/home/icb/bhavishya.nelakuditi/pre_processing/seurat_processig/Pool6/filtered_qc", fig.width=10, fig.height=6)
p = FeatureScatter(Pool6, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = 'lm')
save_plot(p, "/ictstr01/home/icb/bhavishya.nelakuditi/pre_processing/seurat_processig/Pool6/filtered_qc2", fig.width=10, fig.height=6)

```

### normalising data, finding variable features, and scaling data 
### regressing out the cell cycle genes
```{r cars2, echo=TRUE}

# 3. Normalize data ----------
Pool6 <- NormalizeData(Pool6, normalization.method = "LogNormalize", scale.factor = 10000)
#str(Pool6)


# 4. Identify highly variable features --------------
Pool6 <- FindVariableFeatures(Pool6, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(Pool6), 10)

# plot variable features with and without labels
print("highly variable genes plot")
plot1 <- VariableFeaturePlot(Pool6)
LabelPoints(plot = plot1, points = top10, repel = TRUE)
save_plot(plot1, "/ictstr01/home/icb/bhavishya.nelakuditi/pre_processing/seurat_processig/Pool6/highly_variable", fig.width=10, fig.height=6)


# 5. Scaling -------------
#s.genes <- cc.genes$s.genes
#g2m.genes <- cc.genes$g2m.genes
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
all.genes <- rownames(Pool6)
Pool6 <- CellCycleScoring(Pool6, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
Pool6 <- ScaleData(Pool6, vars.to.regress = c("percent.mt","S.Score", "G2M.Score"))

print("Pool6 seurat object after normalising, finding variable features, scaling, and regressing out cell cycle variable effects")
Pool6 = RunPCA(Pool6, features = c(s.genes, g2m.genes))
p = DimPlot(Pool6, reduction = "pca")
save_plot(p, "/ictstr01/home/icb/bhavishya.nelakuditi/pre_processing/seurat_processig/Pool6/dimplot", fig.width=10, fig.height=6)


str(Pool6)


```

### HTO demultiplexing and removing doublets and negative hashtagged cells
```{r cars3, echo=TRUE}
# Add HTO data as a new assay independent from RNA
obj.htos <- allABs.raw$Pool6[, colnames(Pool6)]
Pool6[["HTO"]] <- CreateAssayObject(counts = obj.htos)

# Normalize HTO data, here we use centered log-ratio (CLR) transformation
Pool6 <- NormalizeData(Pool6, assay = "HTO", normalization.method = "CLR")
Pool6 <- HTODemux(Pool6, assay = "HTO", positive.quantile = 0.99)

print("HTO classification table")
table(Pool6$HTO_classification.global)
print("HTO hashtags")
p = RidgePlot(Pool6, assay = "HTO", features = rownames(Pool6[["HTO"]])[1:10], ncol = 5)
save_plot(p, "/ictstr01/home/icb/bhavishya.nelakuditi/pre_processing/seurat_processig/Pool6/HTO_ridge", fig.width=30, fig.height=6)


Idents(Pool6) <- "HTO_classification.global"
Pool6HTO=Pool6[,Pool6$HTO_classification.global!="Doublet"]
Pool6HTO=Pool6HTO[,Pool6HTO$HTO_classification.global!="Negative"]


print("HTO classification table after removing doublets and negative hashtagged cells")
table(Pool6HTO$HTO_classification.global)

rownames(allABs.raw$Pool6)
### adding the conditions as a vector to the seurat object
cellList = colnames(Pool6HTO)
featVec <- vector(mode="character", length=length(cellList))
featVec = Pool6HTO$HTO_classification
featVec[featVec == "HTO-B0251"] = "HF"
featVec[featVec == "HTO-B0252"] = "KS"
featVec[featVec == "HTO-B0253"] = "KS"
featVec[featVec == "HTO-B0254"] = "KS"
featVec[featVec == "HTO-B0256"] = "KS"
featVec[featVec == "HTO-B0257"] = "KS"
featVec[featVec == "HTO-B0258"] = "KS"
featVec[featVec == "HTO-B0259"] = "KS"
featVec[featVec == "HTO-B0260"] = "KS"
Pool6HTO$CSclassification=featVec


cellList = colnames(Pool6HTO)
featVec <- vector(mode="character", length=length(cellList))
featVec = Pool6HTO$HTO_classification
featVec[featVec == "HTO-B0251"] = "HF10"
featVec[featVec == "HTO-B0252"] = "KS29.1"
featVec[featVec == "HTO-B0253"] = "KS36.1"
featVec[featVec == "HTO-B0254"] = "KS31.2"
featVec[featVec == "HTO-B0256"] = "KS20.1"
featVec[featVec == "HTO-B0257"] = "KS12.2"
featVec[featVec == "HTO-B0258"] = "KS18.1"
featVec[featVec == "HTO-B0259"] = "KS17.2"
featVec[featVec == "HTO-B0260"] = "KS33.2"
Pool6HTO$Patient=featVec



# Calculate a tSNE embedding of the HTO data
DefaultAssay(Pool6HTO) <- "HTO"
Pool6HTO <- ScaleData(Pool6HTO, features = rownames(Pool6HTO),
    verbose = FALSE)
Pool6HTO <- RunPCA(Pool6HTO, features = rownames(Pool6HTO), approx = FALSE)
Pool6HTO <- RunTSNE(Pool6HTO, dims = 1:8, perplexity = 100)
print("HTO dimplot")
p = DimPlot(Pool6HTO)
save_plot(p, "/ictstr01/home/icb/bhavishya.nelakuditi/pre_processing/seurat_processig/Pool6/HTO_dimplot", fig.width=10, fig.height=6)


saveRDS(Pool6HTO, "/ictstr01/home/icb/bhavishya.nelakuditi/pre_processing/seurat_processig/Pool6/Pool6HTO.Rds")   
saveRDS(Pool6HTO, "/ictstr01/home/icb/bhavishya.nelakuditi/pre_processing/seurat_processig/Pool6/Pool6HTO.RDS")   

> table(Pool6HTO$HTO_classification.global)

Singlet 
   9784 


