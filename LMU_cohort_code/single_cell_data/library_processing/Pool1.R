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


foldername = dirname(Sys.glob("/mnt/raidtmp/Alejandro/shock_data/home/jovyan/bhavishya-cardiogenic-shock-volume/Pool1/outs/raw_feature_bc_matrix/*.mtx.gz"))
h5file = Read10X(foldername,unique.features = TRUE)
samplename = "Pool1"
allfiles.raw = list()
allABs.raw = list()

allfiles.raw[[samplename]] = h5file$`Gene Expression`
allABs.raw[[samplename]] = h5file$`Antibody Capture`


joint.bcs <- intersect(colnames(allfiles.raw$Pool1), colnames(allABs.raw$Pool1))
allfiles.raw$Pool1 <- allfiles.raw$Pool1[, joint.bcs]
allABs.raw$Pool1 <- as.matrix(allABs.raw$Pool1[, joint.bcs])

print("conditions")
rownames(allABs.raw$Pool1)

```

### filtering data and removing ribosomal genes
```{r filtering, echo=TRUE}
Pool1 <- CreateSeuratObject(counts = allfiles.raw$Pool1, project = "Pool1", min.cells = 3, min.features = 200)
str(Pool1)

### removing ribosomal genes from the data
#rps_rownames = rownames(Pool1[grep('^Rps', rownames(Pool1))])
#rpl_rownames = rownames(Pool1[grep('^Rpl', rownames(Pool1))])
#Pool1 <- Pool1[!(row.names(Pool1) %in% rps_rownames),]
#Pool1 <- Pool1[!(row.names(Pool1) %in% rpl_rownames),]
Pool1[["percent.mt"]] <- PercentageFeatureSet(Pool1, pattern = "^MT-")

## seurat object before filtering
print("seurat object before filtering")
p = VlnPlot(Pool1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
save_plot(p, "/mnt/raidtmp/Alejandro/shock_data/Pool1/scRNAseq/before_filtered_qc", fig.width=10, fig.height=6)
p=FeatureScatter(Pool1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = 'lm')
save_plot(p, "/mnt/raidtmp/Alejandro/shock_data/Pool1/scRNAseq/before_filtered_qc2", fig.width=10, fig.height=6)

### filtering data
Pool1 <- subset(Pool1, subset = nFeature_RNA > 100 & nFeature_RNA < 6000 & nCount_RNA > 100 & nCount_RNA < 15000)
Pool1 <- subset(Pool1, subset = percent.mt < 15)

## seurat object after filtering
print("seurat object after filtering")
p = VlnPlot(Pool1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
save_plot(p, "/mnt/raidtmp/Alejandro/shock_data/Pool1/scRNAseq/filtered_qc", fig.width=10, fig.height=6)
p = FeatureScatter(Pool1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = 'lm')
save_plot(p, "/mnt/raidtmp/Alejandro/shock_data/Pool1/scRNAseq/filtered_qc2", fig.width=10, fig.height=6)

```

### normalising data, finding variable features, and scaling data 
### regressing out the cell cycle genes
```{r cars2, echo=TRUE}

# 3. Normalize data ----------
Pool1 <- NormalizeData(Pool1, normalization.method = "LogNormalize", scale.factor = 10000)
#str(Pool1)


# 4. Identify highly variable features --------------
Pool1 <- FindVariableFeatures(Pool1, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(Pool1), 10)

# plot variable features with and without labels
print("highly variable genes plot")
plot1 <- VariableFeaturePlot(Pool1)
LabelPoints(plot = plot1, points = top10, repel = TRUE)
save_plot(plot1, "/mnt/raidtmp/Alejandro/shock_data/Pool1/scRNAseq/highly_variable", fig.width=10, fig.height=6)


# 5. Scaling -------------
#s.genes <- cc.genes$s.genes
#g2m.genes <- cc.genes$g2m.genes
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
all.genes <- rownames(Pool1)
Pool1 <- CellCycleScoring(Pool1, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
Pool1 <- ScaleData(Pool1, vars.to.regress = c("percent.mt","S.Score", "G2M.Score"))

print("Pool1 seurat object after normalising, finding variable features, scaling, and regressing out cell cycle variable effects")
Pool1 = RunPCA(Pool1, features = c(s.genes, g2m.genes))
p = DimPlot(Pool1, reduction = "pca")
save_plot(p, "/mnt/raidtmp/Alejandro/shock_data/Pool1/scRNAseq/dimplot", fig.width=10, fig.height=6)


str(Pool1)


```

### HTO demultiplexing and removing doublets and negative hashtagged cells
```{r cars3, echo=TRUE}
# Add HTO data as a new assay independent from RNA
obj.htos <- allABs.raw$Pool1[, colnames(Pool1)]
Pool1[["HTO"]] <- CreateAssayObject(counts = obj.htos)

# Normalize HTO data, here we use centered log-ratio (CLR) transformation
Pool1 <- NormalizeData(Pool1, assay = "HTO", normalization.method = "CLR")
Pool1 <- HTODemux(Pool1, assay = "HTO", positive.quantile = 0.99)

print("HTO classification table")
table(Pool1$HTO_classification.global)
print("HTO hashtags")
p = RidgePlot(Pool1, assay = "HTO", features = rownames(Pool1[["HTO"]])[1:10], ncol = 5)
save_plot(p, "/mnt/raidtmp/Alejandro/shock_data/Pool1/scRNAseq/HTO_ridge", fig.width=30, fig.height=6)


Idents(Pool1) <- "HTO_classification.global"
Pool1HTO=Pool1[,Pool1$HTO_classification.global!="Doublet"]
Pool1HTO=Pool1HTO[,Pool1HTO$HTO_classification.global!="Negative"]


print("HTO classification table after removing doublets and negative hashtagged cells")
table(Pool1HTO$HTO_classification.global)

rownames(allABs.raw$Pool1)
### adding the conditions as a vector to the seurat object
cellList = colnames(Pool1HTO)
featVec <- vector(mode="character", length=length(cellList))
featVec = Pool1HTO$HTO_classification
featVec[featVec == "HTO-B0251"] = "KS"
featVec[featVec == "HTO-B0252"] = "HF"
featVec[featVec == "HTO-B0253"] = "KS"
featVec[featVec == "HTO-B0254"] = "KS"
featVec[featVec == "HTO-B0255"] = "KS"
featVec[featVec == "HTO-B0256"] = "KS"
featVec[featVec == "HTO-B0257"] = "KS"
featVec[featVec == "HTO-B0258"] = "KS"
featVec[featVec == "HTO-B0259"] = "KS"
featVec[featVec == "HTO-B0260"] = "KS"
Pool1HTO$CSclassification=featVec


cellList = colnames(Pool1HTO)
featVec <- vector(mode="character", length=length(cellList))
featVec = Pool1HTO$HTO_classification
featVec[featVec == "HTO-B0251"] = "KS22.1"
featVec[featVec == "HTO-B0252"] = "HF15"
featVec[featVec == "HTO-B0253"] = "KS32.1"
featVec[featVec == "HTO-B0254"] = "KS27.1"
featVec[featVec == "HTO-B0255"] = "KS16.2"
featVec[featVec == "HTO-B0256"] = "KS24.3"
featVec[featVec == "HTO-B0257"] = "KS4.2"
featVec[featVec == "HTO-B0258"] = "KS35.1"
featVec[featVec == "HTO-B0259"] = "KS35.2"
featVec[featVec == "HTO-B0260"] = "KS5.2"
Pool1HTO$Patient=featVec



# Calculate a tSNE embedding of the HTO data
DefaultAssay(Pool1HTO) <- "HTO"
Pool1HTO <- ScaleData(Pool1HTO, features = rownames(Pool1HTO),
    verbose = FALSE)
Pool1HTO <- RunPCA(Pool1HTO, features = rownames(Pool1HTO), approx = FALSE)
Pool1HTO <- RunTSNE(Pool1HTO, dims = 1:8, perplexity = 100)
print("HTO dimplot")
p = DimPlot(Pool1HTO)
save_plot(p, "/mnt/raidtmp/Alejandro/shock_data/Pool1/scRNAseq/HTO_dimplot", fig.width=10, fig.height=6)



saveRDS(Pool1HTO, "/mnt/raidtmp/Alejandro/shock_data/Pool1/scRNAseq/Pool1HTO.Rds")   
saveRDS(Pool1HTO, "/mnt/raidtmp/Alejandro/shock_data/Pool1/scRNAseq/Pool1HTO.RDS")   

> table(Pool1HTO$HTO_classification.global)

Singlet 
  12072 