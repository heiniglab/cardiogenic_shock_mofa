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


foldername = dirname(Sys.glob("/mnt/raidtmp/Alejandro/shock_data/home/jovyan/bhavishya-cardiogenic-shock-volume/Pool7/outs/raw_feature_bc_matrix/*.mtx.gz"))
h5file = Read10X(foldername,unique.features = TRUE)
samplename = "Pool7"
allfiles.raw = list()
allABs.raw = list()

allfiles.raw[[samplename]] = h5file$`Gene Expression`
allABs.raw[[samplename]] = h5file$`Antibody Capture`


joint.bcs <- intersect(colnames(allfiles.raw$Pool7), colnames(allABs.raw$Pool7))
allfiles.raw$Pool7 <- allfiles.raw$Pool7[, joint.bcs]
allABs.raw$Pool7 <- as.matrix(allABs.raw$Pool7[, joint.bcs])

print("conditions")
rownames(allABs.raw$Pool7)

```

### filtering data and removing ribosomal genes
```{r filtering, echo=TRUE}
Pool7 <- CreateSeuratObject(counts = allfiles.raw$Pool7, project = "Pool7", min.cells = 3, min.features = 200)
str(Pool7)

### removing ribosomal genes from the data
#rps_rownames = rownames(Pool7[grep('^Rps', rownames(Pool7))])
#rpl_rownames = rownames(Pool7[grep('^Rpl', rownames(Pool7))])
#Pool7 <- Pool7[!(row.names(Pool7) %in% rps_rownames),]
#Pool7 <- Pool7[!(row.names(Pool7) %in% rpl_rownames),]
Pool7[["percent.mt"]] <- PercentageFeatureSet(Pool7, pattern = "^MT-")

## seurat object before filtering
print("seurat object before filtering")
p = VlnPlot(Pool7, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
save_plot(p, "/mnt/raidtmp/Alejandro/shock_data/Pool7/scRNAseq/before_filtered_qc", fig.width=10, fig.height=6)
p=FeatureScatter(Pool7, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = 'lm')
save_plot(p, "/mnt/raidtmp/Alejandro/shock_data/Pool7/scRNAseq/before_filtered_qc2", fig.width=10, fig.height=6)

### filtering data
Pool7 <- subset(Pool7, subset = nFeature_RNA > 100 & nFeature_RNA < 9000 & nCount_RNA > 100 & nCount_RNA < 45000)
Pool7 <- subset(Pool7, subset = percent.mt < 15)

## seurat object after filtering
print("seurat object after filtering")
p = VlnPlot(Pool7, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
save_plot(p, "/mnt/raidtmp/Alejandro/shock_data/Pool7/scRNAseq/filtered_qc", fig.width=10, fig.height=6)
p = FeatureScatter(Pool7, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = 'lm')
save_plot(p, "/mnt/raidtmp/Alejandro/shock_data/Pool7/scRNAseq/filtered_qc2", fig.width=10, fig.height=6)

```

### normalising data, finding variable features, and scaling data 
### regressing out the cell cycle genes
```{r cars2, echo=TRUE}

# 3. Normalize data ----------
Pool7 <- NormalizeData(Pool7, normalization.method = "LogNormalize", scale.factor = 10000)
#str(Pool7)


# 4. Identify highly variable features --------------
Pool7 <- FindVariableFeatures(Pool7, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(Pool7), 10)

# plot variable features with and without labels
print("highly variable genes plot")
plot1 <- VariableFeaturePlot(Pool7)
LabelPoints(plot = plot1, points = top10, repel = TRUE)
save_plot(plot1, "/mnt/raidtmp/Alejandro/shock_data/Pool7/scRNAseq/highly_variable", fig.width=10, fig.height=6)


# 5. Scaling -------------
#s.genes <- cc.genes$s.genes
#g2m.genes <- cc.genes$g2m.genes
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
all.genes <- rownames(Pool7)
Pool7 <- CellCycleScoring(Pool7, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
Pool7 <- ScaleData(Pool7, vars.to.regress = c("percent.mt","S.Score", "G2M.Score"))

print("Pool7 seurat object after normalising, finding variable features, scaling, and regressing out cell cycle variable effects")
Pool7 = RunPCA(Pool7, features = c(s.genes, g2m.genes))
p = DimPlot(Pool7, reduction = "pca")
save_plot(p, "/mnt/raidtmp/Alejandro/shock_data/Pool7/scRNAseq/dimplot", fig.width=10, fig.height=6)


str(Pool7)


```

### HTO demultiplexing and removing doublets and negative hashtagged cells
```{r cars3, echo=TRUE}
# Add HTO data as a new assay independent from RNA
obj.htos <- allABs.raw$Pool7[, colnames(Pool7)]
Pool7[["HTO"]] <- CreateAssayObject(counts = obj.htos)

# Normalize HTO data, here we use centered log-ratio (CLR) transformation
Pool7 <- NormalizeData(Pool7, assay = "HTO", normalization.method = "CLR")
Pool7 <- HTODemux(Pool7, assay = "HTO", positive.quantile = 0.99)

print("HTO classification table")
table(Pool7$HTO_classification.global)
print("HTO hashtags")
p = RidgePlot(Pool7, assay = "HTO", features = rownames(Pool7[["HTO"]])[1:10], ncol = 5)
save_plot(p, "/mnt/raidtmp/Alejandro/shock_data/Pool7/scRNAseq/HTO_ridge", fig.width=30, fig.height=6)


Idents(Pool7) <- "HTO_classification.global"
Pool7HTO=Pool7[,Pool7$HTO_classification.global!="Doublet"]
Pool7HTO=Pool7HTO[,Pool7HTO$HTO_classification.global!="Negative"]


print("HTO classification table after removing doublets and negative hashtagged cells")
table(Pool7HTO$HTO_classification.global)

rownames(allABs.raw$Pool7)
### adding the conditions as a vector to the seurat object
cellList = colnames(Pool7HTO)
featVec <- vector(mode="character", length=length(cellList))
featVec = Pool7HTO$HTO_classification
featVec[featVec == "HTO-B0251"] = "KS"
featVec[featVec == "HTO-B0252"] = "KS"
featVec[featVec == "HTO-B0254"] = "HF"
featVec[featVec == "HTO-B0255"] = "KS"
featVec[featVec == "HTO-B0256"] = "HF"
featVec[featVec == "HTO-B0257"] = "KS"
featVec[featVec == "HTO-B0258"] = "KS"
featVec[featVec == "HTO-B0259"] = "KS"
featVec[featVec == "HTO-B0260"] = "KS"
Pool7HTO$CSclassification=featVec


cellList = colnames(Pool7HTO)
featVec <- vector(mode="character", length=length(cellList))
featVec = Pool7HTO$HTO_classification
featVec[featVec == "HTO-B0251"] = "KS9.2"
featVec[featVec == "HTO-B0252"] = "KS16.3"
featVec[featVec == "HTO-B0254"] = "HF11"
featVec[featVec == "HTO-B0255"] = "KS25.2"
featVec[featVec == "HTO-B0256"] = "HF4"
featVec[featVec == "HTO-B0257"] = "KS19.2"
featVec[featVec == "HTO-B0258"] = "KS33.1"
featVec[featVec == "HTO-B0259"] = "KS17.1"
featVec[featVec == "HTO-B0260"] = "KS3.2"
Pool7HTO$Patient=featVec



# Calculate a tSNE embedding of the HTO data
DefaultAssay(Pool7HTO) <- "HTO"
Pool7HTO <- ScaleData(Pool7HTO, features = rownames(Pool7HTO),
    verbose = FALSE)
Pool7HTO <- RunPCA(Pool7HTO, features = rownames(Pool7HTO), approx = FALSE)
Pool7HTO <- RunTSNE(Pool7HTO, dims = 1:8, perplexity = 100)
print("HTO dimplot")
p = DimPlot(Pool7HTO)
save_plot(p, "/mnt/raidtmp/Alejandro/shock_data/Pool7/scRNAseq/HTO_dimplot", fig.width=10, fig.height=6)


saveRDS(Pool7HTO, "/mnt/raidtmp/Alejandro/shock_data/Pool7/scRNAseq/Pool7HTO.Rds")   
saveRDS(Pool7HTO, "/mnt/raidtmp/Alejandro/shock_data/Pool7/scRNAseq/Pool7HTO.RDS")   

> table(Pool7HTO$HTO_classification.global)

Singlet 
  12204 



