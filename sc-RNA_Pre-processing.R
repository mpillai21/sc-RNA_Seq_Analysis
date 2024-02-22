#Load Seurat
library(Seurat)

#Nomralize and Scale the data
data <- NormalizeData(data)
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(data), 10)
plot1 <- VariableFeaturePlot(data)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
p = plot1 + plot2
ggsave("./Output/Variable_feature.pdf",plot2 , units = "in", height = 5, width = 8, dpi = 72)

data <- ScaleData(data, features = VariableFeatures(data))

#save the normalized and scaled objects
saveRDS(data, "./Objects/data_normalized_scaled_date.rds")

#Linear Dimensionality Reduction
data <- RunPCA(data, features = VariableFeatures(object = data))
DimHeatmap(data, dims = 1:15, cells = 500, balanced = TRUE)
ElbowPlot(data)
ggsave("./Output/PCA_data.pdf", DimHeatmap(data, dims = 1:15, cells = 500, balanced = TRUE) , units = "in", height = 5, width = 8, dpi = 72)

#Perform clustering
data <- FindNeighbors(data, dims = 1:15)
data <- FindClusters(data, resolution = c(0.1, 0.2, 0.3, 0.5, 0.8))
Idents(data) <- "RNA_snn_res.0.3" #based on the communities 

#Non-linear reduction
data <- RunUMAP(data, dims = 1:10)
data <- RunTSNE(data, dims = 1:10)
saveRDS(data, "data_allreductions.rds")

#Plot the reductions
ggsave("./Output/UMAP_date.pdf", DimPlot(data, reduction = "umap", label = TRUE) , units = "in", height = 5, width = 8, dpi = 72)
ggsave("./Output/TSNE_date.pdf", DimPlot(data, reduction = "tsne") , units = "in", height = 5, width = 8, dpi = 72)
