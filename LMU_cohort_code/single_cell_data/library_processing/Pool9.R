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


foldername = dirname(Sys.glob("/mnt/raidtmp/Alejandro/shock_data/home/jovyan/bhavishya-cardiogenic-shock-volume/Pool9/outs/raw_feature_bc_matrix/*.mtx.gz"))
h5file = Read10X(foldername,unique.features = TRUE)
samplename = "Pool9"
allfiles.raw = list()
allABs.raw = list()

allfiles.raw[[samplename]] = h5file$`Gene Expression`
allABs.raw[[samplename]] = h5file$`Antibody Capture`


joint.bcs <- intersect(colnames(allfiles.raw$Pool9), colnames(allABs.raw$Pool9))
allfiles.raw$Pool9 <- allfiles.raw$Pool9[, joint.bcs]
allABs.raw$Pool9 <- as.matrix(allABs.raw$Pool9[, joint.bcs])

print("conditions")
rownames(allABs.raw$Pool9)

```

### filtering data and removing ribosomal genes
```{r filtering, echo=TRUE}
Pool9 <- CreateSeuratObject(counts = allfiles.raw$Pool9, project = "Pool9", min.cells = 3, min.features = 200)
str(Pool9)

### removing ribosomal genes from the data
#rps_rownames = rownames(Pool9[grep('^Rps', rownames(Pool9))])
#rpl_rownames = rownames(Pool9[grep('^Rpl', rownames(Pool9))])
#Pool9 <- Pool9[!(row.names(Pool9) %in% rps_rownames),]
#Pool9 <- Pool9[!(row.names(Pool9) %in% rpl_rownames),]
Pool9[["percent.mt"]] <- PercentageFeatureSet(Pool9, pattern = "^MT-")

## seurat object before filtering
print("seurat object before filtering")
p = VlnPlot(Pool9, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
save_plot(p, "/mnt/raidtmp/Alejandro/shock_data/Pool9/scRNAseq/before_filtered_qc", fig.width=10, fig.height=6)
p=FeatureScatter(Pool9, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = 'lm')
save_plot(p, "/mnt/raidtmp/Alejandro/shock_data/Pool9/scRNAseq/before_filtered_qc2", fig.width=10, fig.height=6)

### filtering data
Pool9 <- subset(Pool9, subset = nFeature_RNA > 100 & nFeature_RNA < 9000 & nCount_RNA > 100 & nCount_RNA < 45000)
Pool9 <- subset(Pool9, subset = percent.mt < 15)

## seurat object after filtering
print("seurat object after filtering")
p = VlnPlot(Pool9, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
save_plot(p, "/mnt/raidtmp/Alejandro/shock_data/Pool9/scRNAseq/filtered_qc", fig.width=10, fig.height=6)
p = FeatureScatter(Pool9, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = 'lm')
save_plot(p, "/mnt/raidtmp/Alejandro/shock_data/Pool9/scRNAseq/filtered_qc2", fig.width=10, fig.height=6)

```

### normalising data, finding variable features, and scaling data 
### regressing out the cell cycle genes
```{r cars2, echo=TRUE}

# 3. Normalize data ----------
Pool9 <- NormalizeData(Pool9, normalization.method = "LogNormalize", scale.factor = 10000)
#str(Pool9)


# 4. Identify highly variable features --------------
Pool9 <- FindVariableFeatures(Pool9, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(Pool9), 10)

# plot variable features with and without labels
print("highly variable genes plot")
plot1 <- VariableFeaturePlot(Pool9)
LabelPoints(plot = plot1, points = top10, repel = TRUE)
save_plot(plot1, "/mnt/raidtmp/Alejandro/shock_data/Pool9/scRNAseq/highly_variable", fig.width=10, fig.height=6)


# 5. Scaling -------------
#s.genes <- cc.genes$s.genes
#g2m.genes <- cc.genes$g2m.genes
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
all.genes <- rownames(Pool9)
Pool9 <- CellCycleScoring(Pool9, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
Pool9 <- ScaleData(Pool9, vars.to.regress = c("percent.mt","S.Score", "G2M.Score"))

print("Pool9 seurat object after normalising, finding variable features, scaling, and regressing out cell cycle variable effects")
Pool9 = RunPCA(Pool9, features = c(s.genes, g2m.genes))
p = DimPlot(Pool9, reduction = "pca")
save_plot(p, "/mnt/raidtmp/Alejandro/shock_data/Pool9/scRNAseq/dimplot", fig.width=10, fig.height=6)


str(Pool9)


```

### HTO demultiplexing and removing doublets and negative hashtagged cells
```{r cars3, echo=TRUE}
# Add HTO data as a new assay independent from RNA
obj.htos <- allABs.raw$Pool9[, colnames(Pool9)]
Pool9[["HTO"]] <- CreateAssayObject(counts = obj.htos)

# Normalize HTO data, here we use centered log-ratio (CLR) transformation
Pool9 <- NormalizeData(Pool9, assay = "HTO", normalization.method = "CLR")
Pool9 <- HTODemux(Pool9, assay = "HTO", positive.quantile = 0.99)

print("HTO classification table")
table(Pool9$HTO_classification.global)
print("HTO hashtags")
p = RidgePlot(Pool9, assay = "HTO", features = rownames(Pool9[["HTO"]])[1:10], ncol = 5)
save_plot(p, "/mnt/raidtmp/Alejandro/shock_data/Pool9/scRNAseq/HTO_ridge", fig.width=30, fig.height=6)


Idents(Pool9) <- "HTO_classification.global"
Pool9HTO=Pool9[,Pool9$HTO_classification.global!="Doublet"]
Pool9HTO=Pool9HTO[,Pool9HTO$HTO_classification.global!="Negative"]


print("HTO classification table after removing doublets and negative hashtagged cells")
table(Pool9HTO$HTO_classification.global)

rownames(allABs.raw$Pool9)
### adding the conditions as a vector to the seurat object
cellList = colnames(Pool9HTO)
featVec <- vector(mode="character", length=length(cellList))
featVec = Pool9HTO$HTO_classification
featVec[featVec == "HTO-B0251"] = "KS"
featVec[featVec == "HTO-B0252"] = "KS"
featVec[featVec == "HTO-B0254"] = "HF"
featVec[featVec == "HTO-B0255"] = "HF"
featVec[featVec == "HTO-B0256"] = "KS"
featVec[featVec == "HTO-B0257"] = "KS"
featVec[featVec == "HTO-B0258"] = "KS"
featVec[featVec == "HTO-B0259"] = "KS"
featVec[featVec == "HTO-B0260"] = "KS"
Pool9HTO$CSclassification=featVec


cellList = colnames(Pool9HTO)
featVec <- vector(mode="character", length=length(cellList))
featVec = Pool9HTO$HTO_classification
featVec[featVec == "HTO-B0251"] = "KS34.1"
featVec[featVec == "HTO-B0252"] = "KS13.3"
featVec[featVec == "HTO-B0254"] = "HF7"
featVec[featVec == "HTO-B0255"] = "HF16"
featVec[featVec == "HTO-B0256"] = "KS11.2"
featVec[featVec == "HTO-B0257"] = "KS10.1"
featVec[featVec == "HTO-B0258"] = "KS10.2"
featVec[featVec == "HTO-B0259"] = "KS14.1"
featVec[featVec == "HTO-B0260"] = "KS13.1"
Pool9HTO$Patient=featVec



# Calculate a tSNE embedding of the HTO data
DefaultAssay(Pool9HTO) <- "HTO"
Pool9HTO <- ScaleData(Pool9HTO, features = rownames(Pool9HTO),
    verbose = FALSE)
Pool9HTO <- RunPCA(Pool9HTO, features = rownames(Pool9HTO), approx = FALSE)
Pool9HTO <- RunTSNE(Pool9HTO, dims = 1:8, perplexity = 100)
print("HTO dimplot")
p = DimPlot(Pool9HTO)
save_plot(p, "/mnt/raidtmp/Alejandro/shock_data/Pool9/scRNAseq/HTO_dimplot", fig.width=10, fig.height=6)


saveRDS(Pool9HTO, "/mnt/raidtmp/Alejandro/shock_data/Pool9/scRNAseq/Pool9HTO.Rds")   
saveRDS(Pool9HTO, "/mnt/raidtmp/Alejandro/shock_data/Pool9/scRNAseq/Pool9HTO.RDS")   

> table(Pool9HTO$HTO_classification.global)

Singlet 
  12639 


