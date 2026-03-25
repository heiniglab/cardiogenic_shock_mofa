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


foldername = dirname(Sys.glob("/mnt/raidtmp/Alejandro/shock_data/home/jovyan/bhavishya-cardiogenic-shock-volume/Pool13/outs/raw_feature_bc_matrix/*.mtx.gz"))
h5file = Read10X(foldername,unique.features = TRUE)
samplename = "Pool13"
allfiles.raw = list()
allABs.raw = list()

allfiles.raw[[samplename]] = h5file$`Gene Expression`
allABs.raw[[samplename]] = h5file$`Antibody Capture`


joint.bcs <- intersect(colnames(allfiles.raw$Pool13), colnames(allABs.raw$Pool13))
allfiles.raw$Pool13 <- allfiles.raw$Pool13[, joint.bcs]
allABs.raw$Pool13 <- as.matrix(allABs.raw$Pool13[, joint.bcs])

print("conditions")
rownames(allABs.raw$Pool13)

```

### filtering data and removing ribosomal genes
```{r filtering, echo=TRUE}
Pool13 <- CreateSeuratObject(counts = allfiles.raw$Pool13, project = "Pool13", min.cells = 3, min.features = 200)
str(Pool13)

### removing ribosomal genes from the data
#rps_rownames = rownames(Pool13[grep('^Rps', rownames(Pool13))])
#rpl_rownames = rownames(Pool13[grep('^Rpl', rownames(Pool13))])
#Pool13 <- Pool13[!(row.names(Pool13) %in% rps_rownames),]
#Pool13 <- Pool13[!(row.names(Pool13) %in% rpl_rownames),]
Pool13[["percent.mt"]] <- PercentageFeatureSet(Pool13, pattern = "^MT-")

## seurat object before filtering
print("seurat object before filtering")
p = VlnPlot(Pool13, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
save_plot(p, "/mnt/raidtmp/Alejandro/shock_data/Pool13/scRNAseq/before_filtered_qc", fig.width=10, fig.height=6)
p=FeatureScatter(Pool13, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = 'lm')
save_plot(p, "/mnt/raidtmp/Alejandro/shock_data/Pool13/scRNAseq/before_filtered_qc2", fig.width=10, fig.height=6)

### filtering data
Pool13 <- subset(Pool13, subset = nFeature_RNA > 100 & nFeature_RNA < 9000 & nCount_RNA > 100 & nCount_RNA < 45000)
Pool13 <- subset(Pool13, subset = percent.mt < 15)

## seurat object after filtering
print("seurat object after filtering")
p = VlnPlot(Pool13, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
save_plot(p, "/mnt/raidtmp/Alejandro/shock_data/Pool13/scRNAseq/filtered_qc", fig.width=10, fig.height=6)
p = FeatureScatter(Pool13, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = 'lm')
save_plot(p, "/mnt/raidtmp/Alejandro/shock_data/Pool13/scRNAseq/filtered_qc2", fig.width=10, fig.height=6)

```

### normalising data, finding variable features, and scaling data 
### regressing out the cell cycle genes
```{r cars2, echo=TRUE}

# 3. Normalize data ----------
Pool13 <- NormalizeData(Pool13, normalization.method = "LogNormalize", scale.factor = 10000)
#str(Pool13)


# 4. Identify highly variable features --------------
Pool13 <- FindVariableFeatures(Pool13, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(Pool13), 10)

# plot variable features with and without labels
print("highly variable genes plot")
plot1 <- VariableFeaturePlot(Pool13)
LabelPoints(plot = plot1, points = top10, repel = TRUE)
save_plot(plot1, "/mnt/raidtmp/Alejandro/shock_data/Pool13/scRNAseq/highly_variable", fig.width=10, fig.height=6)


# 5. Scaling -------------
#s.genes <- cc.genes$s.genes
#g2m.genes <- cc.genes$g2m.genes
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
all.genes <- rownames(Pool13)
Pool13 <- CellCycleScoring(Pool13, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
Pool13 <- ScaleData(Pool13, vars.to.regress = c("percent.mt","S.Score", "G2M.Score"))

print("Pool13 seurat object after normalising, finding variable features, scaling, and regressing out cell cycle variable effects")
Pool13 = RunPCA(Pool13, features = c(s.genes, g2m.genes))
p = DimPlot(Pool13, reduction = "pca")
save_plot(p, "/mnt/raidtmp/Alejandro/shock_data/Pool13/scRNAseq/dimplot", fig.width=10, fig.height=6)


str(Pool13)


```

### HTO demultiplexing and removing doublets and negative hashtagged cells
```{r cars3, echo=TRUE}
# Add HTO data as a new assay independent from RNA
obj.htos <- allABs.raw$Pool13[, colnames(Pool13)]
Pool13[["HTO"]] <- CreateAssayObject(counts = obj.htos)

# Normalize HTO data, here we use centered log-ratio (CLR) transformation
Pool13 <- NormalizeData(Pool13, assay = "HTO", normalization.method = "CLR")
Pool13 <- HTODemux(Pool13, assay = "HTO", positive.quantile = 0.99)

print("HTO classification table")
table(Pool13$HTO_classification.global)
print("HTO hashtags")
p = RidgePlot(Pool13, assay = "HTO", features = rownames(Pool13[["HTO"]])[1:10], ncol = 5)
save_plot(p, "/mnt/raidtmp/Alejandro/shock_data/Pool13/scRNAseq/HTO_ridge", fig.width=30, fig.height=6)


Idents(Pool13) <- "HTO_classification.global"
Pool13HTO=Pool13[,Pool13$HTO_classification.global!="Doublet"]
Pool13HTO=Pool13HTO[,Pool13HTO$HTO_classification.global!="Negative"]


print("HTO classification table after removing doublets and negative hashtagged cells")
table(Pool13HTO$HTO_classification.global)

rownames(allABs.raw$Pool13)
### adding the conditions as a vector to the seurat object
cellList = colnames(Pool13HTO)
featVec <- vector(mode="character", length=length(cellList))
featVec = Pool13HTO$HTO_classification
featVec[featVec == "HTO-B0251"] = "HF"
featVec[featVec == "HTO-B0252"] = "KS"
featVec[featVec == "HTO-B0253"] = "KS"
featVec[featVec == "HTO-B0254"] = "KS"
featVec[featVec == "HTO-B0255"] = "KS"
featVec[featVec == "HTO-B0256"] = "KS"
featVec[featVec == "HTO-B0257"] = "HF"
Pool13HTO$CSclassification=featVec


cellList = colnames(Pool13HTO)
featVec <- vector(mode="character", length=length(cellList))
featVec = Pool13HTO$HTO_classification
featVec[featVec == "HTO-B0251"] = "HF3"
featVec[featVec == "HTO-B0252"] = "KS23.1"
featVec[featVec == "HTO-B0253"] = "KS36.3"
featVec[featVec == "HTO-B0254"] = "KS16.1"
featVec[featVec == "HTO-B0255"] = "KS8.2"
featVec[featVec == "HTO-B0256"] = "KS15.1"
featVec[featVec == "HTO-B0257"] = "HF5"
Pool13HTO$Patient=featVec



# Calculate a tSNE embedding of the HTO data
DefaultAssay(Pool13HTO) <- "HTO"
Pool13HTO <- ScaleData(Pool13HTO, features = rownames(Pool13HTO),
    verbose = FALSE)
Pool13HTO <- RunPCA(Pool13HTO, features = rownames(Pool13HTO), approx = FALSE)
Pool13HTO <- RunTSNE(Pool13HTO, dims = 1:8, perplexity = 100)
print("HTO dimplot")
p = DimPlot(Pool13HTO)
save_plot(p, "/mnt/raidtmp/Alejandro/shock_data/Pool13/scRNAseq/HTO_dimplot", fig.width=10, fig.height=6)


saveRDS(Pool13HTO, "/mnt/raidtmp/Alejandro/shock_data/Pool13/scRNAseq/Pool13HTO.Rds")   
saveRDS(Pool13HTO, "/mnt/raidtmp/Alejandro/shock_data/Pool13/scRNAseq/Pool13HTO.RDS")   

> table(Pool13HTO$HTO_classification.global)

Singlet 
   15682 



