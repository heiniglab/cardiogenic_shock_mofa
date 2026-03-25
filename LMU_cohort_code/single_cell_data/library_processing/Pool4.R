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


foldername = dirname(Sys.glob("/mnt/raidtmp/Alejandro/shock_data/home/jovyan/bhavishya-cardiogenic-shock-volume/Pool4/outs/raw_feature_bc_matrix/*.mtx.gz"))
h5file = Read10X(foldername,unique.features = TRUE)
samplename = "Pool4"
allfiles.raw = list()
allABs.raw = list()

allfiles.raw[[samplename]] = h5file$`Gene Expression`
allABs.raw[[samplename]] = h5file$`Antibody Capture`


joint.bcs <- intersect(colnames(allfiles.raw$Pool4), colnames(allABs.raw$Pool4))
allfiles.raw$Pool4 <- allfiles.raw$Pool4[, joint.bcs]
allABs.raw$Pool4 <- as.matrix(allABs.raw$Pool4[, joint.bcs])

print("conditions")
rownames(allABs.raw$Pool4)

```

### filtering data and removing ribosomal genes
```{r filtering, echo=TRUE}
Pool4 <- CreateSeuratObject(counts = allfiles.raw$Pool4, project = "Pool4", min.cells = 3, min.features = 200)
str(Pool4)

### removing ribosomal genes from the data
#rps_rownames = rownames(Pool4[grep('^Rps', rownames(Pool4))])
#rpl_rownames = rownames(Pool4[grep('^Rpl', rownames(Pool4))])
#Pool4 <- Pool4[!(row.names(Pool4) %in% rps_rownames),]
#Pool4 <- Pool4[!(row.names(Pool4) %in% rpl_rownames),]
Pool4[["percent.mt"]] <- PercentageFeatureSet(Pool4, pattern = "^MT-")

## seurat object before filtering
print("seurat object before filtering")
p = VlnPlot(Pool4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
save_plot(p, "/mnt/raidtmp/Alejandro/shock_data/Pool4/scRNAseq/before_filtered_qc", fig.width=10, fig.height=6)
p=FeatureScatter(Pool4, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = 'lm')
save_plot(p, "/mnt/raidtmp/Alejandro/shock_data/Pool4/scRNAseq/before_filtered_qc2", fig.width=10, fig.height=6)

### filtering data
Pool4 <- subset(Pool4, subset = nFeature_RNA > 100 & nFeature_RNA < 8000 & nCount_RNA > 100 & nCount_RNA < 45000)
Pool4 <- subset(Pool4, subset = percent.mt < 15)

## seurat object after filtering
print("seurat object after filtering")
p = VlnPlot(Pool4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
save_plot(p, "/mnt/raidtmp/Alejandro/shock_data/Pool4/scRNAseq/filtered_qc", fig.width=10, fig.height=6)
p = FeatureScatter(Pool4, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = 'lm')
save_plot(p, "/mnt/raidtmp/Alejandro/shock_data/Pool4/scRNAseq/filtered_qc2", fig.width=10, fig.height=6)

```

### normalising data, finding variable features, and scaling data 
### regressing out the cell cycle genes
```{r cars2, echo=TRUE}

# 3. Normalize data ----------
Pool4 <- NormalizeData(Pool4, normalization.method = "LogNormalize", scale.factor = 10000)
#str(Pool4)


# 4. Identify highly variable features --------------
Pool4 <- FindVariableFeatures(Pool4, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(Pool4), 10)

# plot variable features with and without labels
print("highly variable genes plot")
plot1 <- VariableFeaturePlot(Pool4)
LabelPoints(plot = plot1, points = top10, repel = TRUE)
save_plot(plot1, "/mnt/raidtmp/Alejandro/shock_data/Pool4/scRNAseq/highly_variable", fig.width=10, fig.height=6)


# 5. Scaling -------------
#s.genes <- cc.genes$s.genes
#g2m.genes <- cc.genes$g2m.genes
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
all.genes <- rownames(Pool4)
Pool4 <- CellCycleScoring(Pool4, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
Pool4 <- ScaleData(Pool4, vars.to.regress = c("percent.mt","S.Score", "G2M.Score"))

print("Pool4 seurat object after normalising, finding variable features, scaling, and regressing out cell cycle variable effects")
Pool4 = RunPCA(Pool4, features = c(s.genes, g2m.genes))
p = DimPlot(Pool4, reduction = "pca")
save_plot(p, "/mnt/raidtmp/Alejandro/shock_data/Pool4/scRNAseq/dimplot", fig.width=10, fig.height=6)


str(Pool4)


```

### HTO demultiplexing and removing doublets and negative hashtagged cells
```{r cars3, echo=TRUE}
# Add HTO data as a new assay independent from RNA
obj.htos <- allABs.raw$Pool4[, colnames(Pool4)]
Pool4[["HTO"]] <- CreateAssayObject(counts = obj.htos)

# Normalize HTO data, here we use centered log-ratio (CLR) transformation
Pool4 <- NormalizeData(Pool4, assay = "HTO", normalization.method = "CLR")
Pool4 <- HTODemux(Pool4, assay = "HTO", positive.quantile = 0.99)

print("HTO classification table")
table(Pool4$HTO_classification.global)
print("HTO hashtags")
p = RidgePlot(Pool4, assay = "HTO", features = rownames(Pool4[["HTO"]])[1:10], ncol = 5)
save_plot(p, "/mnt/raidtmp/Alejandro/shock_data/Pool4/scRNAseq/HTO_ridge", fig.width=30, fig.height=6)


Idents(Pool4) <- "HTO_classification.global"
Pool4HTO=Pool4[,Pool4$HTO_classification.global!="Doublet"]
Pool4HTO=Pool4HTO[,Pool4HTO$HTO_classification.global!="Negative"]


print("HTO classification table after removing doublets and negative hashtagged cells")
table(Pool4HTO$HTO_classification.global)

rownames(allABs.raw$Pool4)
### adding the conditions as a vector to the seurat object
cellList = colnames(Pool4HTO)
featVec <- vector(mode="character", length=length(cellList))
featVec = Pool4HTO$HTO_classification
featVec[featVec == "HTO-B0251"] = "KS"
featVec[featVec == "HTO-B0252"] = "KS"
featVec[featVec == "HTO-B0254"] = "KS"
featVec[featVec == "HTO-B0255"] = "KS"
featVec[featVec == "HTO-B0257"] = "HF"
featVec[featVec == "HTO-B0258"] = "HF"
featVec[featVec == "HTO-B0259"] = "KS"
featVec[featVec == "HTO-B0260"] = "KS"
Pool4HTO$CSclassification=featVec


cellList = colnames(Pool4HTO)
featVec <- vector(mode="character", length=length(cellList))
featVec = Pool4HTO$HTO_classification
featVec[featVec == "HTO-B0251"] = "KS29.3"
featVec[featVec == "HTO-B0252"] = "KS35.3"
featVec[featVec == "HTO-B0254"] = "KS12.1"
featVec[featVec == "HTO-B0255"] = "KS25.1"
featVec[featVec == "HTO-B0257"] = "HF17"
featVec[featVec == "HTO-B0258"] = "HF14"
featVec[featVec == "HTO-B0259"] = "KS21.1"
featVec[featVec == "HTO-B0260"] = "KS2.2"
Pool4HTO$Patient=featVec



# Calculate a tSNE embedding of the HTO data
DefaultAssay(Pool4HTO) <- "HTO"
Pool4HTO <- ScaleData(Pool4HTO, features = rownames(Pool4HTO),
    verbose = FALSE)
Pool4HTO <- RunPCA(Pool4HTO, features = rownames(Pool4HTO), approx = FALSE)
Pool4HTO <- RunTSNE(Pool4HTO, dims = 1:8, perplexity = 100)
print("HTO dimplot")
p = DimPlot(Pool4HTO)
save_plot(p, "/mnt/raidtmp/Alejandro/shock_data/Pool4/scRNAseq/HTO_dimplot", fig.width=10, fig.height=6)




> table(Pool4HTO$HTO_classification.global)

Singlet 
   6604 

saveRDS(Pool4HTO, "/mnt/raidtmp/Alejandro/shock_data/Pool4/scRNAseq/Pool4HTO.Rds")   
saveRDS(Pool4HTO, "/mnt/raidtmp/Alejandro/shock_data/Pool4/scRNAseq/Pool4HTO.RDS")   