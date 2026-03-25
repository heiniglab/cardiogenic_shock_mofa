```


Pool1HTO= readRDS("shock_data/Pool1/scRNAseq/Pool1HTO.Rds")
Pool2HTO= readRDS("shock_data/Pool2/scRNAseq/Pool2HTO.Rds")
Pool3HTO= readRDS("shock_data/Pool3/scRNAseq/Pool3HTO.RDS")
Pool4HTO= readRDS("shock_data/Pool4/scRNAseq/Pool4HTO.Rds")
Pool5HTO= readRDS("shock_data/Pool5/scRNAseq/Pool5HTO.Rds")
Pool6HTO= readRDS("shock_data/Pool6/scRNAseq/Pool6HTO.Rds")
Pool7HTO= readRDS("shock_data/Pool7/scRNAseq/Pool7HTO.Rds")
Pool8HTO= readRDS("shock_data/Pool8/scRNAseq/Pool8HTO.Rds")
Pool9HTO= readRDS("shock_data/Pool9/scRNAseq/Pool9HTO.Rds")
Pool11HTO= readRDS("shock_data/Pool11/scRNAseq/Pool11HTO.Rds")
Pool12HTO= readRDS("shock_data/Pool12/scRNAseq/Pool12HTO.Rds")
Pool13HTO= readRDS("shock_data/Pool13/scRNAseq/Pool13HTO.Rds")









### finding cluster specific genes
```{r cars4, echo=TRUE}
DefaultAssay(Pool1HTO) <- "RNA"
DefaultAssay(Pool2HTO) <- "RNA"
DefaultAssay(Pool3HTO) <- "RNA"
DefaultAssay(Pool4HTO) <- "RNA"
DefaultAssay(Pool5HTO) <- "RNA"
DefaultAssay(Pool6HTO) <- "RNA"
DefaultAssay(Pool7HTO) <- "RNA"
DefaultAssay(Pool8HTO) <- "RNA"
DefaultAssay(Pool9HTO) <- "RNA"
DefaultAssay(Pool11HTO) <- "RNA"
DefaultAssay(Pool12HTO) <- "RNA"
DefaultAssay(Pool13HTO) <- "RNA"
# Select the top 2000 most variable features
#Pool1HTO <- FindVariableFeatures(Pool1HTO, selection.method = "vst", nfeatures = 2000)

# Scaling RNA data,
#Pool1HTO <- ScaleData(Pool1HTO, features = rownames(Pool1HTO))

# Run PCA
pool_list <- c("Pool1HTO", "Pool2HTO", "Pool3HTO", "Pool4HTO", "Pool5HTO", "Pool6HTO", "Pool7HTO", "Pool8HTO", "Pool9HTO", "Pool11HTO", "Pool12HTO", "Pool13HTO")

for (pool_name in pool_list) {
  # Access the Seurat object by its name
  pool_object <- get(pool_name)
  
  # Run PCA
  pool_object <- RunPCA(pool_object, features = VariableFeatures(pool_object))
  
  # Run UMAP
  pool_object <- RunUMAP(pool_object, dims = 1:30, reduction.key = "UMAP_")
  
  # Assign the modified object back to its name
  assign(pool_name, pool_object)
}




# Integration with harmony with SCT transform
merged_seurat = merge(Pool1HTO, y=c(Pool2HTO, Pool3HTO, Pool4HTO, Pool5HTO, Pool6HTO, Pool7HTO, Pool8HTO, Pool9HTO, Pool11HTO, Pool12HTO, Pool13HTO))
obj.list <- SplitObject(merged_seurat, split.by = 'orig.ident')
for(i in 1:length(obj.list)){
  obj.list[[i]] <- SCTransform(obj.list[[i]])
}

for(i in 1:length(obj.list)){
  DefaultAssay(obj.list[[i]]) <- "SCT"
}

# Example list of Seurat objects
seurat_list <- list(Pool2HTO, Pool3HTO, Pool4HTO, Pool5HTO, Pool6HTO, Pool7HTO, Pool8HTO, Pool9HTO, Pool11HTO, Pool12HTO, Pool13HTO)

# Get the intersected row names (genes) across all Seurat objects
common_genes <- Reduce(intersect, lapply(seurat_list, rownames))

# Print the common genes
print(common_genes)

seurat_list_filtered <- lapply(seurat_list, function(seurat_obj) {
  subset(seurat_obj, features = common_genes)
})

Pool2HTO <- subset(Pool2HTO, features = common_genes)
Pool3HTO <- subset(Pool3HTO, features = common_genes)
Pool4HTO <- subset(Pool4HTO, features = common_genes)
Pool5HTO <- subset(Pool5HTO, features = common_genes)
Pool6HTO <- subset(Pool6HTO, features = common_genes)
Pool7HTO <- subset(Pool7HTO, features = common_genes)
Pool8HTO <- subset(Pool8HTO, features = common_genes)
Pool9HTO <- subset(Pool9HTO, features = common_genes)
Pool11HTO <- subset(Pool11HTO, features = common_genes)
Pool12HTO <- subset(Pool12HTO, features = common_genes)
Pool13HTO <- subset(Pool13HTO, features = common_genes)

# Find most variable features across samples to integrate
integ_features <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = 3000) 

# Merge normalized samples
merged_seurat <- merge(x = obj.list[[1]],
		       y = obj.list[2:length(obj.list)],
		       merge.data = TRUE)
DefaultAssay(merged_seurat) <- "SCT"

# Manually set variable features of merged Seurat object
VariableFeatures(merged_seurat) <- integ_features

# Calculate PCs using manually set variable features
merged_seurat <- RunPCA(merged_seurat, assay = "SCT", npcs = 50)

harmonized_seurat <- RunHarmony(merged_seurat, 
				group.by.vars = 'orig.ident', 
				reduction = "pca", assay.use = "SCT", reduction.save = "harmony")
harmonized_seurat <- RunUMAP(harmonized_seurat, reduction = "harmony", assay = "SCT", dims = 1:40)
seurat.integrated = harmonized_seurat 


saveRDS(harmonized_seurat,"/ictstr01/home/icb/bhavishya.nelakuditi/pre_processing/seurat_processig/integrated/new/harmonized_seurat.Rds")
##############################
p=DimPlot(seurat.integrated, group.by="orig.ident", split.by="orig.ident", reduction="umap", shuffle = TRUE, ncol=2, raster = FALSE)
ggsave(
  filename = "/ictstr01/home/icb/bhavishya.nelakuditi/pre_processing/seurat_processig/integrated/batch_effect_corrected_wnn1.png", 
  plot = p, 
  width = 20, 
  height = 30, 
  units = "cm"
)
p=DimPlot(seurat.integrated, group.by="orig.ident", reduction="umap", shuffle = TRUE, seed = 1, ncol=3, raster=FALSE)
ggsave(
  filename = "/ictstr01/home/icb/bhavishya.nelakuditi/pre_processing/seurat_processig/integrated/batch_effect_corrected_wnn2.png", 
  plot = p, 
  width = 20, 
  height = 15, 
  units = "cm"
)
p=DimPlot(seurat.integrated, group.by="CSclassification", reduction="umap", shuffle = TRUE, seed = 1, ncol=3, raster=FALSE)
ggsave(
  filename = "/ictstr01/home/icb/bhavishya.nelakuditi/pre_processing/seurat_processig/integrated/ig_plot_classification.png", 
  plot = p, 
  width = 20, 
  height = 15, 
  units = "cm"
)
seurat.integrated <- FindNeighbors(object = seurat.integrated, reduction = "harmony")
seurat.integrated <- FindClusters(seurat.integrated, resolution = 0.2)


# Projecting singlet identities on UMAP  visualization
p=DimPlot(seurat.integrated, group.by = "CSclassification", raster = FALSE)
ggsave("/ictstr01/home/icb/bhavishya.nelakuditi/pre_processing/seurat_processig/integrated/new/dimplotx.png",p,  width=10, height=6)
p=DimPlot(seurat.integrated, shuffle = T, seed = 1, group.by= "CSclassification", split.by= "CSclassification", ncol=3, raster=FALSE)
ggsave("/ictstr01/home/icb/bhavishya.nelakuditi/pre_processing/seurat_processig/integrated/new/dimplot2.png", p,width=20, height=10)
p=DimPlot(seurat.integrated, group.by= "seurat_clusters", raster=FALSE, label = TRUE)
ggsave("/ictstr01/home/icb/bhavishya.nelakuditi/pre_processing/seurat_processig/integrated/new/clusters.png", p, width=10, height=6)

Idents(seurat.integrated) = seurat.integrated$seurat_clusters
seurat.integrated <- PrepSCTFindMarkers(object = seurat.integrated)
'''