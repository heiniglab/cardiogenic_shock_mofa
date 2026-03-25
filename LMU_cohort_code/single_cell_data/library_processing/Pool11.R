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


foldername = dirname(Sys.glob("/mnt/raidtmp/Alejandro/shock_data/home/jovyan/bhavishya-cardiogenic-shock-volume/Pool11/outs/raw_feature_bc_matrix/*.mtx.gz"))
h5file = Read10X(foldername,unique.features = TRUE)
samplename = "Pool11"
allfiles.raw = list()
allABs.raw = list()

allfiles.raw[[samplename]] = h5file$`Gene Expression`
allABs.raw[[samplename]] = h5file$`Antibody Capture`


joint.bcs <- intersect(colnames(allfiles.raw$Pool11), colnames(allABs.raw$Pool11))
allfiles.raw$Pool11 <- allfiles.raw$Pool11[, joint.bcs]
allABs.raw$Pool11 <- as.matrix(allABs.raw$Pool11[, joint.bcs])

print("conditions")
rownames(allABs.raw$Pool11)

```

### filtering data and removing ribosomal genes
```{r filtering, echo=TRUE}
Pool11 <- CreateSeuratObject(counts = allfiles.raw$Pool11, project = "Pool11", min.cells = 3, min.features = 200)
str(Pool11)

### removing ribosomal genes from the data
#rps_rownames = rownames(Pool11[grep('^Rps', rownames(Pool11))])
#rpl_rownames = rownames(Pool11[grep('^Rpl', rownames(Pool11))])
#Pool11 <- Pool11[!(row.names(Pool11) %in% rps_rownames),]
#Pool11 <- Pool11[!(row.names(Pool11) %in% rpl_rownames),]
Pool11[["percent.mt"]] <- PercentageFeatureSet(Pool11, pattern = "^MT-")

## seurat object before filtering
print("seurat object before filtering")
p = VlnPlot(Pool11, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
save_plot(p, "/mnt/raidtmp/Alejandro/shock_data/Pool11/scRNAseq/before_filtered_qc", fig.width=10, fig.height=6)
p=FeatureScatter(Pool11, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = 'lm')
save_plot(p, "/mnt/raidtmp/Alejandro/shock_data/Pool11/scRNAseq/before_filtered_qc2", fig.width=10, fig.height=6)

### filtering data
Pool11 <- subset(Pool11, subset = nFeature_RNA > 100 & nFeature_RNA < 9500 & nCount_RNA > 100 & nCount_RNA < 50000)
Pool11 <- subset(Pool11, subset = percent.mt < 15)

## seurat object after filtering
print("seurat object after filtering")
p = VlnPlot(Pool11, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
save_plot(p, "/mnt/raidtmp/Alejandro/shock_data/Pool11/scRNAseq/filtered_qc", fig.width=10, fig.height=6)
p = FeatureScatter(Pool11, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = 'lm')
save_plot(p, "/mnt/raidtmp/Alejandro/shock_data/Pool11/scRNAseq/filtered_qc2", fig.width=10, fig.height=6)

```

### normalising data, finding variable features, and scaling data 
### regressing out the cell cycle genes
```{r cars2, echo=TRUE}

# 3. Normalize data ----------
Pool11 <- NormalizeData(Pool11, normalization.method = "LogNormalize", scale.factor = 10000)
#str(Pool11)


# 4. Identify highly variable features --------------
Pool11 <- FindVariableFeatures(Pool11, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(Pool11), 10)

# plot variable features with and without labels
print("highly variable genes plot")
plot1 <- VariableFeaturePlot(Pool11)
LabelPoints(plot = plot1, points = top10, repel = TRUE)
save_plot(plot1, "/mnt/raidtmp/Alejandro/shock_data/Pool11/scRNAseq/highly_variable", fig.width=10, fig.height=6)


# 5. Scaling -------------
#s.genes <- cc.genes$s.genes
#g2m.genes <- cc.genes$g2m.genes
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
all.genes <- rownames(Pool11)
Pool11 <- CellCycleScoring(Pool11, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
Pool11 <- ScaleData(Pool11, vars.to.regress = c("percent.mt","S.Score", "G2M.Score"))

print("Pool11 seurat object after normalising, finding variable features, scaling, and regressing out cell cycle variable effects")
Pool11 = RunPCA(Pool11, features = c(s.genes, g2m.genes))
p = DimPlot(Pool11, reduction = "pca")
save_plot(p, "/mnt/raidtmp/Alejandro/shock_data/Pool11/scRNAseq/dimplot", fig.width=10, fig.height=6)


str(Pool11)


```

### HTO demultiplexing and removing doublets and negative hashtagged cells
```{r cars3, echo=TRUE}
# Add HTO data as a new assay independent from RNA
obj.htos <- allABs.raw$Pool11[, colnames(Pool11)]
Pool11[["HTO"]] <- CreateAssayObject(counts = obj.htos)

# Normalize HTO data, here we use centered log-ratio (CLR) transformation
Pool11 <- NormalizeData(Pool11, assay = "HTO", normalization.method = "CLR")
Pool11 <- HTODemux(Pool11, assay = "HTO", positive.quantile = 0.99)

print("HTO classification table")
table(Pool11$HTO_classification.global)
print("HTO hashtags")
p = RidgePlot(Pool11, assay = "HTO", features = rownames(Pool11[["HTO"]])[1:10], ncol = 5)
save_plot(p, "/mnt/raidtmp/Alejandro/shock_data/Pool11/scRNAseq/HTO_ridge", fig.width=30, fig.height=6)


Idents(Pool11) <- "HTO_classification.global"
Pool11HTO=Pool11[,Pool11$HTO_classification.global!="Doublet"]
Pool11HTO=Pool11HTO[,Pool11HTO$HTO_classification.global!="Negative"]


print("HTO classification table after removing doublets and negative hashtagged cells")
table(Pool11HTO$HTO_classification.global)

rownames(allABs.raw$Pool11)
### adding the conditions as a vector to the seurat object
cellList = colnames(Pool11HTO)
featVec <- vector(mode="character", length=length(cellList))
featVec = Pool11HTO$HTO_classification
featVec[featVec == "HTO-B0252"] = "KS"
featVec[featVec == "HTO-B0253"] = "KS"
featVec[featVec == "HTO-B0254"] = "KS"
featVec[featVec == "HTO-B0255"] = "KS"
featVec[featVec == "HTO-B0256"] = "KS"
Pool11HTO$CSclassification=featVec


cellList = colnames(Pool11HTO)
featVec <- vector(mode="character", length=length(cellList))
featVec = Pool11HTO$HTO_classification
featVec[featVec == "HTO-B0252"] = "KS14.2"
featVec[featVec == "HTO-B0253"] = "KS4.1"
featVec[featVec == "HTO-B0254"] = "KS20.2"
featVec[featVec == "HTO-B0255"] = "KS5.1"
featVec[featVec == "HTO-B0256"] = "KS5.3"
Pool11HTO$Patient=featVec



# Calculate a tSNE embedding of the HTO data
DefaultAssay(Pool11HTO) <- "HTO"
Pool11HTO <- ScaleData(Pool11HTO, features = rownames(Pool11HTO),
    verbose = FALSE)
Pool11HTO <- RunPCA(Pool11HTO, features = rownames(Pool11HTO), approx = FALSE)
Pool11HTO <- RunTSNE(Pool11HTO, dims = 1:8, perplexity = 100)
print("HTO dimplot")
p = DimPlot(Pool11HTO)
save_plot(p, "/mnt/raidtmp/Alejandro/shock_data/Pool11/scRNAseq/HTO_dimplot", fig.width=10, fig.height=6)


saveRDS(Pool11HTO, "/mnt/raidtmp/Alejandro/shock_data/Pool11/scRNAseq/Pool11HTO.Rds")   
saveRDS(Pool11HTO, "/mnt/raidtmp/Alejandro/shock_data/Pool11/scRNAseq/Pool11HTO.RDS")   

>table(Pool11HTO$HTO_classification.global)

Singlet 
   8473



