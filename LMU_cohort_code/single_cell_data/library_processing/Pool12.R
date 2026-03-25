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


foldername = dirname(Sys.glob("/mnt/raidtmp/Alejandro/shock_data/home/jovyan/bhavishya-cardiogenic-shock-volume/Pool12/outs/raw_feature_bc_matrix/*.mtx.gz"))
h5file = Read10X(foldername,unique.features = TRUE)
samplename = "Pool12"
allfiles.raw = list()
allABs.raw = list()

allfiles.raw[[samplename]] = h5file$`Gene Expression`
allABs.raw[[samplename]] = h5file$`Antibody Capture`


joint.bcs <- intersect(colnames(allfiles.raw$Pool12), colnames(allABs.raw$Pool12))
allfiles.raw$Pool12 <- allfiles.raw$Pool12[, joint.bcs]
allABs.raw$Pool12 <- as.matrix(allABs.raw$Pool12[, joint.bcs])

print("conditions")
rownames(allABs.raw$Pool12)

```

### filtering data and removing ribosomal genes
```{r filtering, echo=TRUE}
Pool12 <- CreateSeuratObject(counts = allfiles.raw$Pool12, project = "Pool12", min.cells = 3, min.features = 200)
str(Pool12)

### removing ribosomal genes from the data
#rps_rownames = rownames(Pool12[grep('^Rps', rownames(Pool12))])
#rpl_rownames = rownames(Pool12[grep('^Rpl', rownames(Pool12))])
#Pool12 <- Pool12[!(row.names(Pool12) %in% rps_rownames),]
#Pool12 <- Pool12[!(row.names(Pool12) %in% rpl_rownames),]
Pool12[["percent.mt"]] <- PercentageFeatureSet(Pool12, pattern = "^MT-")

## seurat object before filtering
print("seurat object before filtering")
p = VlnPlot(Pool12, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
save_plot(p, "/mnt/raidtmp/Alejandro/shock_data/Pool12/scRNAseq/before_filtered_qc", fig.width=10, fig.height=6)
p=FeatureScatter(Pool12, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = 'lm')
save_plot(p, "/mnt/raidtmp/Alejandro/shock_data/Pool12/scRNAseq/before_filtered_qc2", fig.width=10, fig.height=6)

### filtering data
Pool12 <- subset(Pool12, subset = nFeature_RNA > 100 & nFeature_RNA < 9500 & nCount_RNA > 100 & nCount_RNA < 45000)
Pool12 <- subset(Pool12, subset = percent.mt < 15)

## seurat object after filtering
print("seurat object after filtering")
p = VlnPlot(Pool12, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
save_plot(p, "/mnt/raidtmp/Alejandro/shock_data/Pool12/scRNAseq/filtered_qc", fig.width=10, fig.height=6)
p = FeatureScatter(Pool12, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = 'lm')
save_plot(p, "/mnt/raidtmp/Alejandro/shock_data/Pool12/scRNAseq/filtered_qc2", fig.width=10, fig.height=6)

```

### normalising data, finding variable features, and scaling data 
### regressing out the cell cycle genes
```{r cars2, echo=TRUE}

# 3. Normalize data ----------
Pool12 <- NormalizeData(Pool12, normalization.method = "LogNormalize", scale.factor = 10000)
#str(Pool12)


# 4. Identify highly variable features --------------
Pool12 <- FindVariableFeatures(Pool12, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(Pool12), 10)

# plot variable features with and without labels
print("highly variable genes plot")
plot1 <- VariableFeaturePlot(Pool12)
LabelPoints(plot = plot1, points = top10, repel = TRUE)
save_plot(plot1, "/mnt/raidtmp/Alejandro/shock_data/Pool12/scRNAseq/highly_variable", fig.width=10, fig.height=6)


# 5. Scaling -------------
#s.genes <- cc.genes$s.genes
#g2m.genes <- cc.genes$g2m.genes
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
all.genes <- rownames(Pool12)
Pool12 <- CellCycleScoring(Pool12, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
Pool12 <- ScaleData(Pool12, vars.to.regress = c("percent.mt","S.Score", "G2M.Score"))

print("Pool12 seurat object after normalising, finding variable features, scaling, and regressing out cell cycle variable effects")
Pool12 = RunPCA(Pool12, features = c(s.genes, g2m.genes))
p = DimPlot(Pool12, reduction = "pca")
save_plot(p, "/mnt/raidtmp/Alejandro/shock_data/Pool12/scRNAseq/dimplot", fig.width=10, fig.height=6)


str(Pool12)


```

### HTO demultiplexing and removing doublets and negative hashtagged cells
```{r cars3, echo=TRUE}
# Add HTO data as a new assay independent from RNA
obj.htos <- allABs.raw$Pool12[, colnames(Pool12)]
Pool12[["HTO"]] <- CreateAssayObject(counts = obj.htos)

# Normalize HTO data, here we use centered log-ratio (CLR) transformation
Pool12 <- NormalizeData(Pool12, assay = "HTO", normalization.method = "CLR")
Pool12 <- HTODemux(Pool12, assay = "HTO", positive.quantile = 0.99)

print("HTO classification table")
table(Pool12$HTO_classification.global)
print("HTO hashtags")
p = RidgePlot(Pool12, assay = "HTO", features = rownames(Pool12[["HTO"]])[1:10], ncol = 5)
save_plot(p, "/mnt/raidtmp/Alejandro/shock_data/Pool12/scRNAseq/HTO_ridge", fig.width=30, fig.height=6)


Idents(Pool12) <- "HTO_classification.global"
Pool12HTO=Pool12[,Pool12$HTO_classification.global!="Doublet"]
Pool12HTO=Pool12HTO[,Pool12HTO$HTO_classification.global!="Negative"]


print("HTO classification table after removing doublets and negative hashtagged cells")
table(Pool12HTO$HTO_classification.global)

rownames(allABs.raw$Pool12)
### adding the conditions as a vector to the seurat object
cellList = colnames(Pool12HTO)
featVec <- vector(mode="character", length=length(cellList))
featVec = Pool12HTO$HTO_classification
featVec[featVec == "HTO-B0251"] = "KS"
featVec[featVec == "HTO-B0252"] = "KS"
featVec[featVec == "HTO-B0254"] = "KS"
featVec[featVec == "HTO-B0255"] = "KS"
featVec[featVec == "HTO-B0256"] = "KS"
featVec[featVec == "HTO-B0257"] = "KS"
Pool12HTO$CSclassification=featVec


cellList = colnames(Pool12HTO)
featVec <- vector(mode="character", length=length(cellList))
featVec = Pool12HTO$HTO_classification
featVec[featVec == "HTO-B0251"] = "KS8.1"
featVec[featVec == "HTO-B0252"] = "KS19.3"
featVec[featVec == "HTO-B0254"] = "KS6.3"
featVec[featVec == "HTO-B0255"] = "KS7.1"
featVec[featVec == "HTO-B0256"] = "KS6.2"
featVec[featVec == "HTO-B0257"] = "KS31.1"
Pool12HTO$Patient=featVec



# Calculate a tSNE embedding of the HTO data
DefaultAssay(Pool12HTO) <- "HTO"
Pool12HTO <- ScaleData(Pool12HTO, features = rownames(Pool12HTO),
    verbose = FALSE)
Pool12HTO <- RunPCA(Pool12HTO, features = rownames(Pool12HTO), approx = FALSE)
Pool12HTO <- RunTSNE(Pool12HTO, dims = 1:8, perplexity = 100)
print("HTO dimplot")
p = DimPlot(Pool12HTO)
save_plot(p, "/mnt/raidtmp/Alejandro/shock_data/Pool12/scRNAseq/HTO_dimplot", fig.width=10, fig.height=6)


saveRDS(Pool12HTO, "/mnt/raidtmp/Alejandro/shock_data/Pool12/scRNAseq/Pool12HTO.Rds")   
saveRDS(Pool12HTO, "/mnt/raidtmp/Alejandro/shock_data/Pool12/scRNAseq/Pool12HTO.RDS")   

> table(Pool12HTO$HTO_classification.global)

Singlet 
   13523



