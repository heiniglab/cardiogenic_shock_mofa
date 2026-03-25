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


foldername = dirname(Sys.glob("/mnt/raidtmp/Alejandro/shock_data/home/jovyan/bhavishya-cardiogenic-shock-volume/Pool3/outs/raw_feature_bc_matrix/*.mtx.gz"))
h5file = Read10X(foldername,unique.features = TRUE)
samplename = "Pool3"
allfiles.raw = list()
allABs.raw = list()

allfiles.raw[[samplename]] = h5file$`Gene Expression`
allABs.raw[[samplename]] = h5file$`Antibody Capture`


joint.bcs <- intersect(colnames(allfiles.raw$Pool3), colnames(allABs.raw$Pool3))
allfiles.raw$Pool3 <- allfiles.raw$Pool3[, joint.bcs]
allABs.raw$Pool3 <- as.matrix(allABs.raw$Pool3[, joint.bcs])

print("conditions")
rownames(allABs.raw$Pool3)

```

### filtering data and removing ribosomal genes
```{r filtering, echo=TRUE}
Pool3 <- CreateSeuratObject(counts = allfiles.raw$Pool3, project = "Pool3", min.cells = 3, min.features = 200)
str(Pool3)

### removing ribosomal genes from the data
#rps_rownames = rownames(Pool3[grep('^Rps', rownames(Pool3))])
#rpl_rownames = rownames(Pool3[grep('^Rpl', rownames(Pool3))])
#Pool3 <- Pool3[!(row.names(Pool3) %in% rps_rownames),]
#Pool3 <- Pool3[!(row.names(Pool3) %in% rpl_rownames),]
Pool3[["percent.mt"]] <- PercentageFeatureSet(Pool3, pattern = "^MT-")

## seurat object before filtering
print("seurat object before filtering")
p = VlnPlot(Pool3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
save_plot(p, "/mnt/raidtmp/Alejandro/shock_data/Pool3/scRNAseq/before_filtered_qc", fig.width=10, fig.height=6)
p=FeatureScatter(Pool3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = 'lm')
save_plot(p, "/mnt/raidtmp/Alejandro/shock_data/Pool3/scRNAseq/before_filtered_qc2", fig.width=10, fig.height=6)

### filtering data
Pool3 <- subset(Pool3, subset = nFeature_RNA > 100 & nFeature_RNA < 8000 & nCount_RNA > 100 & nCount_RNA < 45000)
Pool3 <- subset(Pool3, subset = percent.mt < 15)

## seurat object after filtering
print("seurat object after filtering")
p = VlnPlot(Pool3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
save_plot(p, "/mnt/raidtmp/Alejandro/shock_data/Pool3/scRNAseq/filtered_qc", fig.width=10, fig.height=6)
p = FeatureScatter(Pool3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = 'lm')
save_plot(p, "/mnt/raidtmp/Alejandro/shock_data/Pool3/scRNAseq/filtered_qc2", fig.width=10, fig.height=6)

```

### normalising data, finding variable features, and scaling data 
### regressing out the cell cycle genes
```{r cars2, echo=TRUE}

# 3. Normalize data ----------
Pool3 <- NormalizeData(Pool3, normalization.method = "LogNormalize", scale.factor = 10000)
#str(Pool3)


# 4. Identify highly variable features --------------
Pool3 <- FindVariableFeatures(Pool3, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(Pool3), 10)

# plot variable features with and without labels
print("highly variable genes plot")
plot1 <- VariableFeaturePlot(Pool3)
LabelPoints(plot = plot1, points = top10, repel = TRUE)
save_plot(plot1, "/mnt/raidtmp/Alejandro/shock_data/Pool3/scRNAseq/highly_variable", fig.width=10, fig.height=6)


# 5. Scaling -------------
#s.genes <- cc.genes$s.genes
#g2m.genes <- cc.genes$g2m.genes
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
all.genes <- rownames(Pool3)
Pool3 <- CellCycleScoring(Pool3, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
Pool3 <- ScaleData(Pool3, vars.to.regress = c("percent.mt","S.Score", "G2M.Score"))

print("Pool3 seurat object after normalising, finding variable features, scaling, and regressing out cell cycle variable effects")
Pool3 = RunPCA(Pool3, features = c(s.genes, g2m.genes))
p = DimPlot(Pool3, reduction = "pca")
save_plot(p, "/mnt/raidtmp/Alejandro/shock_data/Pool3/scRNAseq/dimplot", fig.width=10, fig.height=6)


str(Pool3)


```

### HTO demultiplexing and removing doublets and negative hashtagged cells
```{r cars3, echo=TRUE}
# Add HTO data as a new assay independent from RNA
obj.htos <- allABs.raw$Pool3[, colnames(Pool3)]
Pool3[["HTO"]] <- CreateAssayObject(counts = obj.htos)

# Normalize HTO data, here we use centered log-ratio (CLR) transformation
Pool3 <- NormalizeData(Pool3, assay = "HTO", normalization.method = "CLR")
Pool3 <- HTODemux(Pool3, assay = "HTO", positive.quantile = 0.99)

print("HTO classification table")
table(Pool3$HTO_classification.global)
print("HTO hashtags")
p = RidgePlot(Pool3, assay = "HTO", features = rownames(Pool3[["HTO"]])[1:10], ncol = 5)
save_plot(p, "/mnt/raidtmp/Alejandro/shock_data/Pool3/scRNAseq/HTO_ridge", fig.width=30, fig.height=6)


Idents(Pool3) <- "HTO_classification.global"
Pool3HTO=Pool3[,Pool3$HTO_classification.global!="Doublet"]
Pool3HTO=Pool3HTO[,Pool3HTO$HTO_classification.global!="Negative"]


print("HTO classification table after removing doublets and negative hashtagged cells")
table(Pool3HTO$HTO_classification.global)

rownames(allABs.raw$Pool3)
### adding the conditions as a vector to the seurat object
cellList = colnames(Pool3HTO)
featVec <- vector(mode="character", length=length(cellList))
featVec = Pool3HTO$HTO_classification
featVec[featVec == "HTO-B0251"] = "KS"
featVec[featVec == "HTO-B0254"] = "KS"
featVec[featVec == "HTO-B0255"] = "KS"
featVec[featVec == "HTO-B0256"] = "KS"
featVec[featVec == "HTO-B0257"] = "KS"
featVec[featVec == "HTO-B0258"] = "KS"
featVec[featVec == "HTO-B0259"] = "HF"
featVec[featVec == "HTO-B0260"] = "KS"
Pool3HTO$CSclassification=featVec


cellList = colnames(Pool3HTO)
featVec <- vector(mode="character", length=length(cellList))
featVec = Pool3HTO$HTO_classification
featVec[featVec == "HTO-B0251"] = "KS2.1"
featVec[featVec == "HTO-B0254"] = "KS26.1"
featVec[featVec == "HTO-B0255"] = "KS12.3"
featVec[featVec == "HTO-B0256"] = "KS34.3"
featVec[featVec == "HTO-B0257"] = "KS11.1"
featVec[featVec == "HTO-B0258"] = "KS24.1"
featVec[featVec == "HTO-B0259"] = "HF9"
featVec[featVec == "HTO-B0260"] = "KS9.1"
Pool3HTO$Patient=featVec



# Calculate a tSNE embedding of the HTO data
DefaultAssay(Pool3HTO) <- "HTO"
Pool3HTO <- ScaleData(Pool3HTO, features = rownames(Pool3HTO),
    verbose = FALSE)
Pool3HTO <- RunPCA(Pool3HTO, features = rownames(Pool3HTO), approx = FALSE)
Pool3HTO <- RunTSNE(Pool3HTO, dims = 1:8, perplexity = 100)
print("HTO dimplot")
p = DimPlot(Pool3HTO)
save_plot(p, "/mnt/raidtmp/Alejandro/shock_data/Pool3/scRNAseq/HTO_dimplot", fig.width=10, fig.height=6)


saveRDS(Pool3HTO, "/mnt/raidtmp/Alejandro/shock_data/Pool3/scRNAseq/Pool3HTO.RDS")

Pool3HTO = readRDS("/mnt/raidtmp/Alejandro/shock_data/Pool3/scRNAseq/Pool3HTO.RDS")

> table(Pool3HTO$HTO_classification.global)

Singlet 
   5568 
> 