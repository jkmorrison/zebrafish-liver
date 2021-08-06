library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(cowplot)

Zf.all.liver <- readRDS("AllZf_liver_clustering.rds")

HSCs.and.partners <- WhichCells(Zf.all.liver, idents = c("zEC/HSC","zInf mac", "zNK/T", "zAgr2+", "zApln+", "zNon-inf mac"))

HSC.subset <- subset(Zf.all.liver, cells = HSCs.and.partners, do.clean = TRUE, do.scale = TRUE)

DefaultAssay(HSC.subset) <- "RNA"
HSC.subset <- NormalizeData(object = HSC.subset, normalization.method = "LogNormalize", 
                           scale.factor = 10000)
HSC.subset <- FindVariableFeatures(object = HSC.subset, mean.function = ExpMean, dispersion.function = LogVMR)

DefaultAssay(HSC.subset) <- "integrated"
HSC.subset.clustering <- ScaleData(object = HSC.subset)
HSC.subset.clustering <- RunPCA(object = HSC.subset.clustering, pc.genes = HSC.subset.clustering@var.genes, do.print = TRUE, pcs.print = 1:10, 
                    genes.print = 5)
ElbowPlot(object = HSC.subset.clustering)
DimHeatmap(HSC.subset.clustering, dims = 1:15, cells = 500, balanced = TRUE)

HSC.subset.clustering <- RunUMAP(HSC.subset.clustering, reduction = "pca", dims = 1:15)
HSC.subset.clustering <- FindNeighbors(HSC.subset.clustering, reduction = "pca", dims = 1:15)
HSC.subset.clustering <- FindClusters(HSC.subset.clustering, resolution = 0.8)


#visualization of clustering - Figure 3A
p2 <- DimPlot(HSC.subset.clustering, reduction = "umap", label = TRUE, repel = TRUE, label.size = 7, pt.size = 0.75) +
  theme(legend.text = element_text(size = 22))
p2
tiff("plot.tiff", units="in", width = 10, height = 7, res = 225)
p2
dev.off()

#visualization of clustering - split by genotype - Figure S2B
test <- HSC.subset.clustering
test$group <- factor(x = HSC.subset.clustering$group, levels = c("WT", "MPI_MT"))
p1 <- DimPlot(test, reduction = "umap", split.by = "group", label = FALSE, pt.size = 0.75) +
  theme(legend.position = "bottom", legend.text = element_text(size = 18))
p1
tiff("plot.tiff", units = "in", width = 12, height = 9, res = 225)
p1
dev.off()

saveRDS(HSC.subset.clustering, "AllZf_HSC_subset.rds")

HSC.subset.clustering <- readRDS("AllZf_HSC_subset.rds")

#number of cells per cluster by sample
write.csv(table(Idents(HSC.subset.clustering), HSC.subset.clustering$orig.ident), "AllZf_HSC_subset_cell_numbers_bysample.csv")

#identify cell type markers for each cluster and store in .csv
DefaultAssay(HSC.subset.clustering) <- "RNA"
recluster.markers <- FindAllMarkers(HSC.subset.clustering, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
recluster.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
write.csv(recluster.markers, "AllZf_HSC_subset_cluster_markers.csv")

#Feature Plot for EC and HSC marker gene expression - Fig S2 D-E
FeaturePlot(HSC.subset.clustering, features = c("hand2"), min.cutoff = "q9")

#Heatmap for nCount_RNA - Fig S2F
tiff("plot.tiff", units = "in", width = 9, height = 6, res = 250)
VlnPlot(HSC.subset.clustering, features = "nCount_RNA", pt.size = 0.05) +
  theme(axis.text = element_text(size = 16, color = "black")) +
  NoLegend()
dev.off()

#heatmap of marker genes for each cell type - Fig S2C
#cluster markers
markers.to.plot <- c("cxcl19", "il1b", "mfap4", 
                     "aoc2.1", "igfbp7", "ponzr1", 
                     "traf1","ccr9a", "tnfrsf9b", 
                     "ier2a", "socs3a", "cldn5b",
                     "nppb", "serpine1", "sele", 
                     "ncf1", "cxcr4b", "cfp", 
                     "ccl38.6", "trac", "dusp2", 
                     "sparc", "colec11", "ccl25b", 
                     "lygl1", "cd74a", "fcer1gl", 
                     "lgals2b", "agr2", "icn", 
                     "apln", "gp1bb", "hbegfb", 
                     "onecutl", "aqp10a", "abcb11b")

cell.limit <- WhichCells(object = HSC.subset.clustering, downsample = 100)
zf.HSC.sub.100cells <- subset(HSC.subset.clustering, cells = cell.limit, do.clean = TRUE, do.scale = TRUE)

DefaultAssay(zf.HSC.sub.100cells) <- "RNA"
all.set <- rownames(zf.HSC.sub.100cells)
zf.HSC.sub.100cells <- ScaleData(zf.HSC.sub.100cells, display.progress = F,features = all.set)

library(viridis)

tiff("plot.tiff", units = "in", width = 10, height = 14, res = 225)
DoHeatmap(zf.HSC.sub.100cells, features = c(markers.to.plot), size = 6.5, angle = 90) +
  NoLegend() +
  scale_fill_viridis() +
  theme(axis.text.y = element_text(size = 18, color = "black" , face = "italic"))
dev.off()

#HSC and EC cluster marker genes - Fig 3B
markers.to.plot <- c("kdrl", "fli1a", "oit3", "ptprb", "ramp2", "cdh5", "fabp11a",
                     "cldn5b", "sele",
                     "nppb", "rgs3a", "sox7", "etv2",
                     "hand2", "col1a1b","ifitm1", "ccl25b", "sparc", "steap4", "colec11")

tiff("plot.tiff", units = "in", width = 11, height = 6, res = 225)
DotPlot(HSC.subset.clustering, features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8) + 
  RotatedAxis() +
  theme(axis.text = element_text(size = 14, color = "black"),
        axis.text.x = element_text(face = "italic", size = 16))
dev.off()

#HSC differential gene expression mpi MT vs WT
fish.new2 <- HSC.subset.clustering
fish.new2$celltype.group <- paste(Idents(fish.new2), fish.new2$group, sep = "_")
fish.new2$celltype <- Idents(fish.new2)
Idents(fish.new2) <- "celltype.group"
cluster.DGE <- FindMarkers(fish.new2, ident.1 = "zHSC_MPI_MT", ident.2 = "zHSC_WT", verbose = FALSE)
head(cluster.DGE, n = 15)
write.csv(cluster.DGE, "Zf_HSC_DGE_MPIMTvWT.csv")



