library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(cowplot)


#Load datasets
WT2.data <- Read10X(data.dir = "Liver scRNA-seq/WT2/")
WT4.data <- Read10X(data.dir = "Liver scRNA-seq/WT4/")
WT5.data <- Read10X(data.dir = "Liver scRNA-seq/WT5/")

#Initialize Seurat object with raw data for each sample
WT1 <- CreateSeuratObject(counts = WT2.data, project = "WT1", min.cells = 5)
WT1$group <- "WT"

WT2 <- CreateSeuratObject(counts = WT4.data, project = "WT2", min.cells = 5)
WT2$group <- "WT"

WT3 <- CreateSeuratObject(counts = WT5.data, project = "WT3", min.cells = 5)
WT3$group <- "WT"

#QC visualization of features and counts with Violin plot
tiff("violin.tiff", units = "in", width = 6, heigh = 4, res = 200)
VlnPlot(WT1, features = "nCount_RNA", pt.size = 0.05) + 
  geom_hline(yintercept = 150, color = "red") + 
  geom_hline(yintercept = 200, color = "blue")
dev.off()

#FeatureScatter is typically used to visualize feature-feature relationships - 
#used here for the remove of cells with too low of features that they are debris or have 
#too many features and are likely doublets
plot1 <- FeatureScatter(WT2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1

#subset Seurat objects with nFeature or nCount limits to get rid of low quality cells and doublets
WT1 <- subset(WT1, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & nCount_RNA > 200)
WT2 <- subset(WT2, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & nCount_RNA > 200)
WT3 <- subset(WT3, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & nCount_RNA > 200)

#check for proper subsetting
VlnPlot(WT1, features = "nFeature_RNA", pt.size = 0.005)

#calculate mitochondrial RNA % for each cell in each Seurat object and plot 
#cells with high mitochondrial % indicate dead/dying/broken cells
WT1[["percent.mt"]] <- PercentageFeatureSet(WT1, pattern = "^mt-")
WT2[["percent.mt"]] <- PercentageFeatureSet(WT2, pattern = "^mt-")
WT3[["percent.mt"]] <- PercentageFeatureSet(WT3, pattern = "^mt-")

#check for mitochondrial transcript percentage distribution
#VlnPlot(WT1, features = "percent.mt", pt.size = 0.005)

#plot1 <- FeatureScatter(WT2, feature1 = "nCount_RNA", feature2 = "percent.mt")

#Enforce mitochondrial cutoffs - this is tissue dependent
WT1 <- subset(WT1, subset = percent.mt < 50)
WT2 <- subset(WT2, subset = percent.mt < 50)
WT3 <- subset(WT3, subset = percent.mt < 50)

#Double check Seurat object for proper subsetting
#VlnPlot(WT2, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), pt.size = 0.005)

#Normalize data in each Seurat object and  identify variable features
WT1 <- NormalizeData(WT1, verbose = FALSE)
WT1 <- FindVariableFeatures(WT1, selection.method = "vst", nfeatures = 2000)

WT2 <- NormalizeData(WT2, verbose = FALSE)
WT2 <- FindVariableFeatures(WT2, selection.method = "vst", nfeatures = 2000)

WT3 <- NormalizeData(WT3, verbose = FALSE)
WT3 <- FindVariableFeatures(WT3, selection.method = "vst", nfeatures = 2000)

#Test samples to determine PC number to be used 
test.sample <- WT2
test.sample <- ScaleData(test.sample, verbose = FALSE)
test.sample <- RunPCA(test.sample, npcs = 50, verbose = FALSE)
ElbowPlot(object = test.sample)

#perform integration with all samples for integrated clustering
anchors <- FindIntegrationAnchors(object.list = list(WT1, WT2, WT3), dims = 1:15)
fish.all <- IntegrateData(anchorset = anchors, dims = 1:15)

#perform integrated analysis
DefaultAssay(fish.all) <- "integrated"
fish.new <- ScaleData(fish.all, verbose = FALSE)
fish.new <- RunPCA(fish.new, npcs = 30, verbose = FALSE)
ElbowPlot(object = fish.new)
DimHeatmap(fish.new, dims = 1:15, cells = 500, balanced = TRUE)
DimHeatmap(fish.new, dims = 16:30, cells = 500, balanced = TRUE)

#tSNE and clustering
fish.new <- RunUMAP(fish.new, reduction = "pca", dims = 1:15)
fish.new <- FindNeighbors(fish.new, reduction = "pca", dims = 1:15)
fish.new <- FindClusters(fish.new, resolution = .8)

#save .rds file
saveRDS(fish.new, "WT_Zebrafish_Liver_Atlas.rds")

#Read in saved RDS
fish.new <- readRDS("WT_Zebrafish_Liver_Atlas.rds")

#visualization of clustering - Figure 1A
p2 <- DimPlot(fish.new, reduction = "umap", label = TRUE, repel = TRUE, label.size = 7, pt.size = 0.75) +
  theme(legend.position = "bottom", legend.text = element_text(size = 22))
p2
tiff("plot.tiff", units="in", width = 11, height = 9, res = 225)
p2
dev.off()

#visualization of clustering - Figure S1A
p1 <- DimPlot(fish.new, reduction = "umap", group.by = "orig.ident", shuffle = TRUE,
              label = FALSE, pt.size = 0.75) +
  theme(legend.text = element_text(size = 20), legend.position = "bottom", legend.box.just = "center")
p1
tiff("plot.tiff", units="in", width = 9.5, height = 8, res = 225)
p1
dev.off()

#identify conserved cell type markers for each cluster and store in .csv - Dataset S1
DefaultAssay(fish.new) <- "RNA"
liver.markers <- FindAllMarkers(fish.new, only.pos = FALSE, min.pct = 0.25, log2fc.threshold = 0.25)
liver.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
write.csv(liver.markers, "WT_Zf_liver_atlas_cluster_markers.csv")

#export .csv file containing number of cells per cluster by sample - Dataset S2
write.csv(table(Idents(fish.new), fish.new$orig.ident), "WT_Zf_liver_atlas_cell_numbers_bysample.csv")

#rescale so that gene labels are found in RNA scaled.data slot
DefaultAssay(fish.new) <- "RNA"
all.set <- rownames(fish.new)
fish.new <- ScaleData(fish.new, display.progress = F,features = all.set)

#Figure 1B - heatmap of marker genes for each cell type
#cluster markers
markers.to.plot <- c("fabp10a", "c3a.1", "apobb.1", "cp", "dbi", "aldob", "ahcy", 
                     "ldhbb", "diabloa", "pnp4b", "rnasel2", "crp3", "cbln8", 
                     "apoa1a", "eno3", "serpina1l", "p4ha3", "ldlra", "ctage5", 
                     "itln3", "crp4", "hamp", "anxa4", "lgals2b", "tm4sf4", 
                     "socs3a", "igfbp7", "colec11", "il1b", "cxcl19", "mfap4", 
                     "sla2", "dusp2", "trac",  
                     "hbaa2", "cahz", "nt5c2l1", "gramd1bb", "gadd45bb", "klf6a",
                     "marco", "cd74a", "grn1",
                     "agr2", "tcnl", "icn", "apln", "hbegfb", "gp1bb", 
                     "mmp13a", "mpx", "cdaa")

#marker expression - Heatmap
cell.limit <- WhichCells(object = fish.new, downsample = 100)
WTzf.liver.100cells <- subset(fish.new, cells = cell.limit, do.clean = TRUE, do.scale = TRUE)

DefaultAssay(WTzf.liver.100cells) <- "RNA"
all.set <- rownames(WTzf.liver.100cells)
WTzf.liver.100cells <- ScaleData(WTzf.liver.100cells, display.progress = F,features = all.set)

library(viridis)

tiff("plot.tiff", units = "in", width = 10, height = 14, res = 225)
DoHeatmap(WTzf.liver.100cells, features = c(markers.to.plot), size = 6.5, angle = 90) +
  NoLegend() +
  scale_fill_viridis() +
  theme(axis.text.y = element_text(size = 18, color = "black" , face = "italic"))
dev.off()

#feature plot for cell type marker gene expression - Figure S1C-D, S1F-K
FeaturePlot(fish.new, features = c("mpx"), min.cutoff = "q9")

#Heatmap for nCount_RNA - Figure S1E
tiff("plot.tiff", units = "in", width = 9, height = 6, res = 250)
VlnPlot(fish.new, features = "nCount_RNA", pt.size = 0.05) +
  theme(axis.text = element_text(size = 16, color = "black")) +
  NoLegend()
dev.off()



