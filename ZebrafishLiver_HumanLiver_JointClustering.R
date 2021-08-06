library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(cowplot)

#load our zebrafish liver scRNA-seq data 
fish.liver <- readRDS("WTFishLiver_HumanOrthologs_04-07-21.rds")

fish.liver$group <- "Fish"

#load Bader lab human liver scRNA-seq data
load("HumanLiver.Rdata")
human.liver <- HumanLiverSeurat
human.liver$group <- "Human"

#switch default assay to RNA 
DefaultAssay(fish.liver) <- "RNA"
DefaultAssay(human.liver) <- "RNA"

#combine fish and human liver scRNA-seq datasets 
combined.liver <- c(fish.liver, human.liver)

#normalize data and determine variable features 
combined.liver <- lapply(X = combined.liver, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = combined.liver)

fish.human.anchors <- FindIntegrationAnchors(object.list = combined.liver, 
                                             anchor.features = features)
fish.human.joint.liver <- IntegrateData(anchors = fish.human.anchors)

#perform integrated analysis
DefaultAssay(fish.human.joint.liver) <- "integrated"
joint.liver <- ScaleData(fish.human.joint.liver, verbose = FALSE)
joint.liver <- RunPCA(joint.liver, npcs = 30, verbose = FALSE)
ElbowPlot(object = joint.liver)
DimHeatmap(joint.liver, dims = 1:15, cells = 500, balanced = TRUE)
DimHeatmap(joint.liver, dims = 16:30, cells = 500, balanced = TRUE)

joint.liver <- RunUMAP(joint.liver, reduction = "pca", dims = 1:15)
joint.liver <- FindNeighbors(joint.liver, reduction = "pca", dims = 1:15)
joint.liver <- FindClusters(joint.liver, resolution = 0.8)

#visualization of clustering - Fig 2A
p2.joint <- DimPlot(joint.liver, reduction = "umap", label = TRUE, repel = TRUE, label.size = 7, pt.size = 0.75) +
  theme(legend.position = "bottom", legend.text = element_text(size = 22))
p2.joint
tiff("plot.tiff", units="in", width = 11, height = 9, res = 225)
p2.joint
dev.off()

#total joint clustering figure - Fig 2C
p1.joint <- DimPlot(joint.liver, reduction = "umap", split.by = "group", label = FALSE, pt.size = 0.75) +
  theme(text = element_text(size = 24, face = "plain"), legend.position = "bottom", legend.text = element_text(size = 22),
        axis.text = element_text(size = 16), axis.title = element_text(size = 16))
p1.joint
tiff("plot.tiff", units = "in", width = 12, height = 9, res = 225)
p1.joint
dev.off()

saveRDS(joint.liver, "Zf_Human_joint_liver_atlas.rds")
joint.liver <- readRDS("Zf_Human_joint_liver_atlas.rds")

#identify cell type markers for each cluster and store in .csv - Dataset S4
DefaultAssay(joint.liver) <- "RNA"
cluster.markers <- FindAllMarkers(joint.liver, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
cluster.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
write.csv(cluster.markers, "Zf_Human_joint_liver_cluster_markers.csv")

#number of cells per cluster per group - Dataset S5
write.csv(table(Idents(joint.liver), joint.liver$group), "Zf_Human_joint_liver_cell_numbers_bygroup.csv")

#heatmap of marker gene expression - Fig 2B
#cluster markers
markers.to.plot <- c("RNASE4", "CP", "C9", "FGA", "FGG", "SPINK1", 
                     "ALDOB", "AHCY", "FABP3", "TAT", "TF", "ADH1B", 
                     "TTR", "DCXR", "ANG", "DBI", "GAMT", "APOC2", 
                     "CCL5", "GNLY", "CD3D", "MARCO", "HMOX1", "MAFB", 
                     "SLC38A4", "F5", "ITIH2", "NKG7", "CMC1", "XCL2", 
                     "ANXA4", "ELF3", "ATF3", "MGP", "FCN2", "DNASE1L3", 
                     "SOCS3", "IGFBP7", "CLDN5", "IGKC", "IGLC2", "IGHG1",
                     "S100A8", "S100A9", "LYZ", "HBE1", "CREG1", "NT5DC4", 
                     "APLN", "PFN1", "TLN1", "STEAP4", "COLEC11", "SPARC", 
                     "AGR2", "TM4SF4", "PDZK1IP1", "STMN1", "TYMS","HMGB2")

#marker expression - Heatmap
cell.limit <- WhichCells(object = joint.liver, downsample = 150)
joint.liver.150cells <- subset(joint.liver, cells = cell.limit, do.clean = TRUE, do.scale = TRUE)

DefaultAssay(joint.liver.150cells) <- "RNA"
all.set <- rownames(joint.liver.150cells)
joint.liver.150cells <- ScaleData(joint.liver.150cells, display.progress = F,features = all.set)

library(viridis)

tiff("plot.tiff", units = "in", width = 10, height = 16, res = 225)
DoHeatmap(joint.liver.150cells, features = c(markers.to.plot), size = 6.5, angle = 90) + 
  NoLegend() +
  scale_fill_viridis() +
  theme(axis.text.y = element_text(size = 18, color = "black", face = "italic"))
dev.off()

#Dot plot for conserved markers for cell types of interest - Fig 2E
HSC.and.partners <- WhichCells(joint.liver, idents = c("zhNon-inf mac",
                                                       "zhNK",
                                                       "zhEC1",
                                                       "zhEC2",
                                                       "zhInf mac",
                                                       "zhEC/HSC",
                                                       "zhChol2",
                                                       zhgdTs))

HSC.and.partners.subset <- subset(joint.liver, cells = HSC.and.partners, do.clean = TRUE, do.scale = TRUE)

markers.to.plot <- c("HMOX1", "MAFB", "CTSB",
                     "XCL2", "DUSP2", "RUNX3",
                     "ID1", "IGFBP7", "CLDN5", "OIT3", "EGFL7", "GPM6A", "PTPRB",
                     "CSTA", "SRGN", "ARPC1B",
                     "COLEC11", "SPARC", "TAGLN", "COL1A1", "COL1A2", 
                     "ELF3", "ANXA4", "TM4SF4",
                     "STMN1", "TOP2A", "AURKB", "BIRC5", "NUSAP1", "CENPF")

tiff("plot.tiff", units = "in", width = 14, height = 6.5, res = 225)
DotPlot(HSC.and.partners.subset, features = markers.to.plot, cols = c("orangered2", "mediumturquoise", "white", "white"), dot.scale = 8, 
        split.by = "group") + 
  RotatedAxis() +
  theme(axis.text = element_text(size = 16, color = "black"),
        axis.text.x = element_text(face = "italic"))
dev.off()



