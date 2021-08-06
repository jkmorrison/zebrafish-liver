library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(cowplot)

#load WT and MPI MT zebrafish liver scRNA-seq data 
WTfish.liver <- readRDS("WTzebrafishliver_03-18-21.rds")
MPIMTfish.liver <- readRDS("MPIMTzebrafishliver_03-24-21.rds")

#switch default assay to RNA 
DefaultAssay(WTfish.liver) <- "RNA"
DefaultAssay(MPIMTfish.liver) <- "RNA"

#combine fish and human liver scRNA-seq datasets 
combined.liver <- c(WTfish.liver, MPIMTfish.liver)

#normalize data and determine variable features 
combined.liver <- lapply(X = combined.liver, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = combined.liver)

allfish.anchors <- FindIntegrationAnchors(object.list = combined.liver, 
                                             anchor.features = features)
allfish.liver <- IntegrateData(anchors = allfish.anchors)

#perform integrated analysis
DefaultAssay(allfish.liver) <- "integrated"
zebrafish.liver.all <- ScaleData(allfish.liver, verbose = FALSE)
zebrafish.liver.all <- RunPCA(zebrafish.liver.all, npcs = 30, verbose = FALSE)
ElbowPlot(object = zebrafish.liver.all)
DimHeatmap(zebrafish.liver.all, dims = 1:15, cells = 500, balanced = TRUE)
DimHeatmap(zebrafish.liver.all, dims = 16:30, cells = 500, balanced = TRUE)

zebrafish.liver.all <- RunUMAP(zebrafish.liver.all, reduction = "pca", dims = 1:15)
zebrafish.liver.all <- FindNeighbors(zebrafish.liver.all, reduction = "pca", dims = 1:15)
zebrafish.liver.all <- FindClusters(zebrafish.liver.all, resolution = 0.8)

#visualization of clustering - Figure S2A
p2 <- DimPlot(zebrafish.liver.all, reduction = "umap", label = TRUE, repel = TRUE, label.size = 7, pt.size = 0.75) +
  theme(legend.position = "bottom", legend.text = element_text(size = 22))
p2
tiff("plot.tiff", units="in", width = 11, height = 9, res = 225)
p2
dev.off()

#save .rds file
saveRDS(zebrafish.liver.all, "Zf_all_liver.rds")

#read in .rds file
zebrafish.liver.all <- readRDS("Zf_all_liver.rds")

#identify cell type markers for each cluster and store in .csv
DefaultAssay(zebrafish.liver.all) <- "RNA"
recluster.markers <- FindAllMarkers(zebrafish.liver.all, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
recluster.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
write.csv(recluster.markers, "Zfall_liver_cluster_markers.csv")


