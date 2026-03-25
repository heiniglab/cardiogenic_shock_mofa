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


foldername = dirname(Sys.glob("/mnt/raidtmp/Alejandro/shock_data/home/jovyan/bhavishya-cardiogenic-shock-volume/Pool2/outs/raw_feature_bc_matrix/*.mtx.gz"))
h5file = Read10X(foldername,unique.features = TRUE)
samplename = "Pool2"
allfiles.raw = list()
allABs.raw = list()

allfiles.raw[[samplename]] = h5file$`Gene Expression`
allABs.raw[[samplename]] = h5file$`Antibody Capture`


joint.bcs <- intersect(colnames(allfiles.raw$Pool2), colnames(allABs.raw$Pool2))
allfiles.raw$Pool2 <- allfiles.raw$Pool2[, joint.bcs]
allABs.raw$Pool2 <- as.matrix(allABs.raw$Pool2[, joint.bcs])

print("conditions")
rownames(allABs.raw$Pool2)

```

### filtering data and removing ribosomal genes
```{r filtering, echo=TRUE}
Pool2 <- CreateSeuratObject(counts = allfiles.raw$Pool2, project = "Pool2", min.cells = 3, min.features = 200)
str(Pool2)

### removing ribosomal genes from the data
#rps_rownames = rownames(Pool2[grep('^Rps', rownames(Pool2))])
#rpl_rownames = rownames(Pool2[grep('^Rpl', rownames(Pool2))])
#Pool2 <- Pool2[!(row.names(Pool2) %in% rps_rownames),]
#Pool2 <- Pool2[!(row.names(Pool2) %in% rpl_rownames),]
Pool2[["percent.mt"]] <- PercentageFeatureSet(Pool2, pattern = "^MT-")

## seurat object before filtering
print("seurat object before filtering")
p = VlnPlot(Pool2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
save_plot(p, "/mnt/raidtmp/Alejandro/shock_data/Pool2/scRNAseq/before_filtered_qc", fig.width=10, fig.height=6)
p=FeatureScatter(Pool2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = 'lm')
save_plot(p, "/mnt/raidtmp/Alejandro/shock_data/Pool2/scRNAseq/before_filtered_qc2", fig.width=10, fig.height=6)

### filtering data
Pool2 <- subset(Pool2, subset = nFeature_RNA > 100 & nFeature_RNA < 8000 & nCount_RNA > 100 & nCount_RNA < 45000)
Pool2 <- subset(Pool2, subset = percent.mt < 15)

## seurat object after filtering
print("seurat object after filtering")
p = VlnPlot(Pool2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
save_plot(p, "/mnt/raidtmp/Alejandro/shock_data/Pool2/scRNAseq/filtered_qc", fig.width=10, fig.height=6)
p = FeatureScatter(Pool2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = 'lm')
save_plot(p, "/mnt/raidtmp/Alejandro/shock_data/Pool2/scRNAseq/filtered_qc2", fig.width=10, fig.height=6)

```

### normalising data, finding variable features, and scaling data 
### regressing out the cell cycle genes
```{r cars2, echo=TRUE}

# 3. Normalize data ----------
Pool2 <- NormalizeData(Pool2, normalization.method = "LogNormalize", scale.factor = 10000)
#str(Pool2)


# 4. Identify highly variable features --------------
Pool2 <- FindVariableFeatures(Pool2, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(Pool2), 10)

# plot variable features with and without labels
print("highly variable genes plot")
plot1 <- VariableFeaturePlot(Pool2)
LabelPoints(plot = plot1, points = top10, repel = TRUE)
save_plot(plot1, "/mnt/raidtmp/Alejandro/shock_data/Pool2/scRNAseq/highly_variable", fig.width=10, fig.height=6)


# 5. Scaling -------------
#s.genes <- cc.genes$s.genes
#g2m.genes <- cc.genes$g2m.genes
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
all.genes <- rownames(Pool2)
Pool2 <- CellCycleScoring(Pool2, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
Pool2 <- ScaleData(Pool2, vars.to.regress = c("percent.mt","S.Score", "G2M.Score"))

print("Pool2 seurat object after normalising, finding variable features, scaling, and regressing out cell cycle variable effects")
Pool2 = RunPCA(Pool2, features = c(s.genes, g2m.genes))
p = DimPlot(Pool2, reduction = "pca")
save_plot(p, "/mnt/raidtmp/Alejandro/shock_data/Pool2/scRNAseq/dimplot", fig.width=10, fig.height=6)


str(Pool2)


```

### HTO demultiplexing and removing doublets and negative hashtagged cells
```{r cars3, echo=TRUE}
# Add HTO data as a new assay independent from RNA
obj.htos <- allABs.raw$Pool2[, colnames(Pool2)]
Pool2[["HTO"]] <- CreateAssayObject(counts = obj.htos)

# Normalize HTO data, here we use centered log-ratio (CLR) transformation
Pool2 <- NormalizeData(Pool2, assay = "HTO", normalization.method = "CLR")
Pool2 <- HTODemux(Pool2, assay = "HTO", positive.quantile = 0.99)

print("HTO classification table")
table(Pool2$HTO_classification.global)
print("HTO hashtags")
p = RidgePlot(Pool2, assay = "HTO", features = rownames(Pool2[["HTO"]])[1:10], ncol = 5)
save_plot(p, "/mnt/raidtmp/Alejandro/shock_data/Pool2/scRNAseq/HTO_ridge", fig.width=30, fig.height=6)


Idents(Pool2) <- "HTO_classification.global"
Pool2HTO=Pool2[,Pool2$HTO_classification.global!="Doublet"]
Pool2HTO=Pool2HTO[,Pool2HTO$HTO_classification.global!="Negative"]


print("HTO classification table after removing doublets and negative hashtagged cells")
table(Pool2HTO$HTO_classification.global)

rownames(allABs.raw$Pool2)
### adding the conditions as a vector to the seurat object
cellList = colnames(Pool2HTO)
featVec <- vector(mode="character", length=length(cellList))
featVec = Pool2HTO$HTO_classification
featVec[featVec == "HTO-B0251"] = "KS"
featVec[featVec == "HTO-B0252"] = "KS"
featVec[featVec == "HTO-B0253"] = "KS"
featVec[featVec == "HTO-B0254"] = "KS"
featVec[featVec == "HTO-B0255"] = "KS"
featVec[featVec == "HTO-B0256"] = "KS"
featVec[featVec == "HTO-B0257"] = "KS"
featVec[featVec == "HTO-B0259"] = "KS"
featVec[featVec == "HTO-B0260"] = "KS"
Pool2HTO$CSclassification=featVec


cellList = colnames(Pool2HTO)
featVec <- vector(mode="character", length=length(cellList))
featVec = Pool2HTO$HTO_classification
featVec[featVec == "HTO-B0251"] = "KS13.2"
featVec[featVec == "HTO-B0252"] = "KS11.3"
featVec[featVec == "HTO-B0253"] = "KS24.2"
featVec[featVec == "HTO-B0254"] = "KS7.2"
featVec[featVec == "HTO-B0255"] = "KS23.2"
featVec[featVec == "HTO-B0256"] = "KS33.3"
featVec[featVec == "HTO-B0257"] = "KS1.1"
featVec[featVec == "HTO-B0259"] = "KS27.2"
featVec[featVec == "HTO-B0260"] = "KS22.2"
Pool2HTO$Patient=featVec



# Calculate a tSNE embedding of the HTO data
DefaultAssay(Pool2HTO) <- "HTO"
Pool2HTO <- ScaleData(Pool2HTO, features = rownames(Pool2HTO),
    verbose = FALSE)
Pool2HTO <- RunPCA(Pool2HTO, features = rownames(Pool2HTO), approx = FALSE)
Pool2HTO <- RunTSNE(Pool2HTO, dims = 1:8, perplexity = 100)
print("HTO dimplot")
p = DimPlot(Pool2HTO)
save_plot(p, "/mnt/raidtmp/Alejandro/shock_data/Pool2/scRNAseq/HTO_dimplot", fig.width=10, fig.height=6)




> table(Pool2HTO$HTO_classification.global)

Singlet 
   9163 


saveRDS(Pool2HTO, "/mnt/raidtmp/Alejandro/shock_data/Pool2/scRNAseq/Pool2HTO.Rds")   