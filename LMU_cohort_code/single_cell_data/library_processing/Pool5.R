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


foldername = dirname(Sys.glob("/mnt/raidtmp/Alejandro/shock_data/home/jovyan/bhavishya-cardiogenic-shock-volume/Pool5/outs/raw_feature_bc_matrix/*.mtx.gz"))
h5file = Read10X(foldername,unique.features = TRUE)
samplename = "Pool5"
allfiles.raw = list()
allABs.raw = list()

allfiles.raw[[samplename]] = h5file$`Gene Expression`
allABs.raw[[samplename]] = h5file$`Antibody Capture`


joint.bcs <- intersect(colnames(allfiles.raw$Pool5), colnames(allABs.raw$Pool5))
allfiles.raw$Pool5 <- allfiles.raw$Pool5[, joint.bcs]
allABs.raw$Pool5 <- as.matrix(allABs.raw$Pool5[, joint.bcs])

print("conditions")
rownames(allABs.raw$Pool5)

```

### filtering data and removing ribosomal genes
```{r filtering, echo=TRUE}
Pool5 <- CreateSeuratObject(counts = allfiles.raw$Pool5, project = "Pool5", min.cells = 3, min.features = 200)
str(Pool5)

### removing ribosomal genes from the data
#rps_rownames = rownames(Pool5[grep('^Rps', rownames(Pool5))])
#rpl_rownames = rownames(Pool5[grep('^Rpl', rownames(Pool5))])
#Pool5 <- Pool5[!(row.names(Pool5) %in% rps_rownames),]
#Pool5 <- Pool5[!(row.names(Pool5) %in% rpl_rownames),]
Pool5[["percent.mt"]] <- PercentageFeatureSet(Pool5, pattern = "^MT-")

## seurat object before filtering
print("seurat object before filtering")
p = VlnPlot(Pool5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
save_plot(p, "/mnt/raidtmp/Alejandro/shock_data/Pool5/scRNAseq/before_filtered_qc", fig.width=10, fig.height=6)
p=FeatureScatter(Pool5, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = 'lm')
save_plot(p, "/mnt/raidtmp/Alejandro/shock_data/Pool5/scRNAseq/before_filtered_qc2", fig.width=10, fig.height=6)

### filtering data
Pool5 <- subset(Pool5, subset = nFeature_RNA > 100 & nFeature_RNA < 8000 & nCount_RNA > 100 & nCount_RNA < 45000)
Pool5 <- subset(Pool5, subset = percent.mt < 15)

## seurat object after filtering
print("seurat object after filtering")
p = VlnPlot(Pool5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
save_plot(p, "/mnt/raidtmp/Alejandro/shock_data/Pool5/scRNAseq/filtered_qc", fig.width=10, fig.height=6)
p = FeatureScatter(Pool5, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = 'lm')
save_plot(p, "/mnt/raidtmp/Alejandro/shock_data/Pool5/scRNAseq/filtered_qc2", fig.width=10, fig.height=6)

```

### normalising data, finding variable features, and scaling data 
### regressing out the cell cycle genes
```{r cars2, echo=TRUE}

# 3. Normalize data ----------
Pool5 <- NormalizeData(Pool5, normalization.method = "LogNormalize", scale.factor = 10000)
#str(Pool5)


# 4. Identify highly variable features --------------
Pool5 <- FindVariableFeatures(Pool5, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(Pool5), 10)

# plot variable features with and without labels
print("highly variable genes plot")
plot1 <- VariableFeaturePlot(Pool5)
LabelPoints(plot = plot1, points = top10, repel = TRUE)
save_plot(plot1, "/mnt/raidtmp/Alejandro/shock_data/Pool5/scRNAseq/highly_variable", fig.width=10, fig.height=6)


# 5. Scaling -------------
#s.genes <- cc.genes$s.genes
#g2m.genes <- cc.genes$g2m.genes
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
all.genes <- rownames(Pool5)
Pool5 <- CellCycleScoring(Pool5, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
Pool5 <- ScaleData(Pool5, vars.to.regress = c("percent.mt","S.Score", "G2M.Score"))

print("Pool5 seurat object after normalising, finding variable features, scaling, and regressing out cell cycle variable effects")
Pool5 = RunPCA(Pool5, features = c(s.genes, g2m.genes))
p = DimPlot(Pool5, reduction = "pca")
save_plot(p, "/mnt/raidtmp/Alejandro/shock_data/Pool5/scRNAseq/dimplot", fig.width=10, fig.height=6)


str(Pool5)


```

### HTO demultiplexing and removing doublets and negative hashtagged cells
```{r cars3, echo=TRUE}
# Add HTO data as a new assay independent from RNA
obj.htos <- allABs.raw$Pool5[, colnames(Pool5)]
Pool5[["HTO"]] <- CreateAssayObject(counts = obj.htos)

# Normalize HTO data, here we use centered log-ratio (CLR) transformation
Pool5 <- NormalizeData(Pool5, assay = "HTO", normalization.method = "CLR")
Pool5 <- HTODemux(Pool5, assay = "HTO", positive.quantile = 0.99)

print("HTO classification table")
table(Pool5$HTO_classification.global)
print("HTO hashtags")
p = RidgePlot(Pool5, assay = "HTO", features = rownames(Pool5[["HTO"]])[1:10], ncol = 5)
save_plot(p, "/mnt/raidtmp/Alejandro/shock_data/Pool5/scRNAseq/HTO_ridge", fig.width=30, fig.height=6)


Idents(Pool5) <- "HTO_classification.global"
Pool5HTO=Pool5[,Pool5$HTO_classification.global!="Doublet"]
Pool5HTO=Pool5HTO[,Pool5HTO$HTO_classification.global!="Negative"]


print("HTO classification table after removing doublets and negative hashtagged cells")
table(Pool5HTO$HTO_classification.global)

rownames(allABs.raw$Pool5)
### adding the conditions as a vector to the seurat object
cellList = colnames(Pool5HTO)
featVec <- vector(mode="character", length=length(cellList))
featVec = Pool5HTO$HTO_classification
featVec[featVec == "HTO-B0251"] = "KS"
featVec[featVec == "HTO-B0252"] = "KS"
featVec[featVec == "HTO-B0253"] = "KS"
featVec[featVec == "HTO-B0254"] = "KS"
featVec[featVec == "HTO-B0256"] = "KS"
featVec[featVec == "HTO-B0257"] = "KS"
featVec[featVec == "HTO-B0258"] = "HF"
featVec[featVec == "HTO-B0259"] = "HF"
featVec[featVec == "HTO-B0260"] = "KS"
Pool5HTO$CSclassification=featVec


cellList = colnames(Pool5HTO)
featVec <- vector(mode="character", length=length(cellList))
featVec = Pool5HTO$HTO_classification
featVec[featVec == "HTO-B0251"] = "KS6.1"
featVec[featVec == "HTO-B0252"] = "KS21.2"
featVec[featVec == "HTO-B0253"] = "KS36.2"
featVec[featVec == "HTO-B0254"] = "KS2.3"
featVec[featVec == "HTO-B0256"] = "KS26.2"
featVec[featVec == "HTO-B0257"] = "KS25.3"
featVec[featVec == "HTO-B0258"] = "HF13"
featVec[featVec == "HTO-B0259"] = "HF6"
featVec[featVec == "HTO-B0260"] = "KS17.3"
Pool5HTO$Patient=featVec



# Calculate a tSNE embedding of the HTO data
DefaultAssay(Pool5HTO) <- "HTO"
Pool5HTO <- ScaleData(Pool5HTO, features = rownames(Pool5HTO),
    verbose = FALSE)
Pool5HTO <- RunPCA(Pool5HTO, features = rownames(Pool5HTO), approx = FALSE)
Pool5HTO <- RunTSNE(Pool5HTO, dims = 1:8, perplexity = 100)
print("HTO dimplot")
p = DimPlot(Pool5HTO)
save_plot(p, "/mnt/raidtmp/Alejandro/shock_data/Pool5/scRNAseq/HTO_dimplot", fig.width=10, fig.height=6)


saveRDS(Pool5HTO, "/mnt/raidtmp/Alejandro/shock_data/Pool5/scRNAseq/Pool5HTO.Rds")   
saveRDS(Pool5HTO, "/mnt/raidtmp/Alejandro/shock_data/Pool5/scRNAseq/Pool5HTO.RDS")   

> table(Pool5HTO$HTO_classification.global)

Singlet 
  11096 

