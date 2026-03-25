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


foldername = dirname(Sys.glob("/mnt/raidtmp/Alejandro/shock_data/home/jovyan/bhavishya-cardiogenic-shock-volume/Pool8/outs/raw_feature_bc_matrix/*.mtx.gz"))
h5file = Read10X(foldername,unique.features = TRUE)
samplename = "Pool8"
allfiles.raw = list()
allABs.raw = list()

allfiles.raw[[samplename]] = h5file$`Gene Expression`
allABs.raw[[samplename]] = h5file$`Antibody Capture`


joint.bcs <- intersect(colnames(allfiles.raw$Pool8), colnames(allABs.raw$Pool8))
allfiles.raw$Pool8 <- allfiles.raw$Pool8[, joint.bcs]
allABs.raw$Pool8 <- as.matrix(allABs.raw$Pool8[, joint.bcs])

print("conditions")
rownames(allABs.raw$Pool8)

```

### filtering data and removing ribosomal genes
```{r filtering, echo=TRUE}
Pool8 <- CreateSeuratObject(counts = allfiles.raw$Pool8, project = "Pool8", min.cells = 3, min.features = 200)
str(Pool8)

### removing ribosomal genes from the data
#rps_rownames = rownames(Pool8[grep('^Rps', rownames(Pool8))])
#rpl_rownames = rownames(Pool8[grep('^Rpl', rownames(Pool8))])
#Pool8 <- Pool8[!(row.names(Pool8) %in% rps_rownames),]
#Pool8 <- Pool8[!(row.names(Pool8) %in% rpl_rownames),]
Pool8[["percent.mt"]] <- PercentageFeatureSet(Pool8, pattern = "^MT-")

## seurat object before filtering
print("seurat object before filtering")
p = VlnPlot(Pool8, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
save_plot(p, "/mnt/raidtmp/Alejandro/shock_data/Pool8/scRNAseq/before_filtered_qc", fig.width=10, fig.height=6)
p=FeatureScatter(Pool8, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = 'lm')
save_plot(p, "/mnt/raidtmp/Alejandro/shock_data/Pool8/scRNAseq/before_filtered_qc2", fig.width=10, fig.height=6)

### filtering data
Pool8 <- subset(Pool8, subset = nFeature_RNA > 100 & nFeature_RNA < 9000 & nCount_RNA > 100 & nCount_RNA < 45000)
Pool8 <- subset(Pool8, subset = percent.mt < 15)

## seurat object after filtering
print("seurat object after filtering")
p = VlnPlot(Pool8, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
save_plot(p, "/mnt/raidtmp/Alejandro/shock_data/Pool8/scRNAseq/filtered_qc", fig.width=10, fig.height=6)
p = FeatureScatter(Pool8, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = 'lm')
save_plot(p, "/mnt/raidtmp/Alejandro/shock_data/Pool8/scRNAseq/filtered_qc2", fig.width=10, fig.height=6)

```

### normalising data, finding variable features, and scaling data 
### regressing out the cell cycle genes
```{r cars2, echo=TRUE}

# 3. Normalize data ----------
Pool8 <- NormalizeData(Pool8, normalization.method = "LogNormalize", scale.factor = 10000)
#str(Pool8)


# 4. Identify highly variable features --------------
Pool8 <- FindVariableFeatures(Pool8, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(Pool8), 10)

# plot variable features with and without labels
print("highly variable genes plot")
plot1 <- VariableFeaturePlot(Pool8)
LabelPoints(plot = plot1, points = top10, repel = TRUE)
save_plot(plot1, "/mnt/raidtmp/Alejandro/shock_data/Pool8/scRNAseq/highly_variable", fig.width=10, fig.height=6)


# 5. Scaling -------------
#s.genes <- cc.genes$s.genes
#g2m.genes <- cc.genes$g2m.genes
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
all.genes <- rownames(Pool8)
Pool8 <- CellCycleScoring(Pool8, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
Pool8 <- ScaleData(Pool8, vars.to.regress = c("percent.mt","S.Score", "G2M.Score"))

print("Pool8 seurat object after normalising, finding variable features, scaling, and regressing out cell cycle variable effects")
Pool8 = RunPCA(Pool8, features = c(s.genes, g2m.genes))
p = DimPlot(Pool8, reduction = "pca")
save_plot(p, "/mnt/raidtmp/Alejandro/shock_data/Pool8/scRNAseq/dimplot", fig.width=10, fig.height=6)


str(Pool8)


```

### HTO demultiplexing and removing doublets and negative hashtagged cells
```{r cars3, echo=TRUE}
# Add HTO data as a new assay independent from RNA
obj.htos <- allABs.raw$Pool8[, colnames(Pool8)]
Pool8[["HTO"]] <- CreateAssayObject(counts = obj.htos)

# Normalize HTO data, here we use centered log-ratio (CLR) transformation
Pool8 <- NormalizeData(Pool8, assay = "HTO", normalization.method = "CLR")
Pool8 <- HTODemux(Pool8, assay = "HTO", positive.quantile = 0.99)

print("HTO classification table")
table(Pool8$HTO_classification.global)
print("HTO hashtags")
p = RidgePlot(Pool8, assay = "HTO", features = rownames(Pool8[["HTO"]])[1:10], ncol = 5)
save_plot(p, "/mnt/raidtmp/Alejandro/shock_data/Pool8/scRNAseq/HTO_ridge", fig.width=30, fig.height=6)


Idents(Pool8) <- "HTO_classification.global"
Pool8HTO=Pool8[,Pool8$HTO_classification.global!="Doublet"]
Pool8HTO=Pool8HTO[,Pool8HTO$HTO_classification.global!="Negative"]


print("HTO classification table after removing doublets and negative hashtagged cells")
table(Pool8HTO$HTO_classification.global)

rownames(allABs.raw$Pool8)
### adding the conditions as a vector to the seurat object
cellList = colnames(Pool8HTO)
featVec <- vector(mode="character", length=length(cellList))
featVec = Pool8HTO$HTO_classification
featVec[featVec == "HTO-B0252"] = "HF"
featVec[featVec == "HTO-B0254"] = "KS"
featVec[featVec == "HTO-B0255"] = "KS"
featVec[featVec == "HTO-B0256"] = "HF"
featVec[featVec == "HTO-B0257"] = "KS"
featVec[featVec == "HTO-B0258"] = "KS"
featVec[featVec == "HTO-B0259"] = "KS"
featVec[featVec == "HTO-B0260"] = "KS"
Pool8HTO$CSclassification=featVec


cellList = colnames(Pool8HTO)
featVec <- vector(mode="character", length=length(cellList))
featVec = Pool8HTO$HTO_classification
featVec[featVec == "HTO-B0252"] = "HF12"
featVec[featVec == "HTO-B0254"] = "KS34.2"
featVec[featVec == "HTO-B0255"] = "KS29.2"
featVec[featVec == "HTO-B0256"] = "HF2"
featVec[featVec == "HTO-B0257"] = "KS26.3"
featVec[featVec == "HTO-B0258"] = "KS31.3"
featVec[featVec == "HTO-B0259"] = "KS19.1"
featVec[featVec == "HTO-B0260"] = "KS3.1"
Pool8HTO$Patient=featVec



# Calculate a tSNE embedding of the HTO data
DefaultAssay(Pool8HTO) <- "HTO"
Pool8HTO <- ScaleData(Pool8HTO, features = rownames(Pool8HTO),
    verbose = FALSE)
Pool8HTO <- RunPCA(Pool8HTO, features = rownames(Pool8HTO), approx = FALSE)
Pool8HTO <- RunTSNE(Pool8HTO, dims = 1:8, perplexity = 100)
print("HTO dimplot")
p = DimPlot(Pool8HTO)
save_plot(p, "/mnt/raidtmp/Alejandro/shock_data/Pool8/scRNAseq/HTO_dimplot", fig.width=10, fig.height=6)


saveRDS(Pool8HTO, "/mnt/raidtmp/Alejandro/shock_data/Pool8/scRNAseq/Pool8HTO.Rds")   
saveRDS(Pool8HTO, "/mnt/raidtmp/Alejandro/shock_data/Pool8/scRNAseq/Pool8HTO.RDS")   

> table(Pool8HTO$HTO_classification.global)

Singlet 
   8461 



