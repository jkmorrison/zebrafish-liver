library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(cowplot)

#load our zebrafish liver scRNA-seq data 
fish.HSC.subset <- readRDS("Zebrafish_HSC_subset_clustering_HumanOrthologs.rds")
fish.HSC.subset$group <- "Fish"

fish.HSC.subset <- RenameIdents(fish.HSC.subset, `0` = "zInf mac",
                                `1` = "zEC1", `2` = "NzK/T1", 
                                `3` = "zEC2", `4` = "zEC-HSC", 
                                `5` = "zEC3", `6` = "zNeu", 
                                `7` =  "zNK/T2", `8` = "zHSC",
                                `9` = "zNon-inf mac", `10` = "zAgr2+",
                                `11` = "zApln+", `12` = "zEC/Hep") 

#subset fish ECs and HSCs
fish.cells <- WhichCells(fish.HSC.subset, idents = c("zEC1", "zEC2",
                                                  "zEC-HSC", "zEC3",
                                                  "zHSC"))
fish.HSCs.ECs <- subset(fish.HSC.subset, cells = fish.cells, do.clean = TRUE, do.scale = TRUE)

#load Bader lab human liver scRNA-seq data
load("HumanLiver.Rdata")
human.liver <- HumanLiverSeurat
human.liver$group <- "Human"

#subset human ECs and HSCs
human.cells <- WhichCells(human.liver, idents = c("11", "12", "13", "20"))
human.HSCs.ECs <- subset(human.liver, cells = human.cells, do.clean = TRUE, do.scale = TRUE)

#switch default assay to RNA 
DefaultAssay(fish.HSCs.ECs) <- "RNA"
DefaultAssay(human.HSCs.ECs) <- "RNA"

#combine fish and human liver scRNA-seq datasets 
combined.HSCs.ECs <- c(fish.HSCs.ECs, human.HSCs.ECs)

#normalize data and determine variable features 
combined.HSCs.ECs <- lapply(X = combined.HSCs.ECs, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = combined.HSCs.ECs)

fish.human.anchors <- FindIntegrationAnchors(object.list = combined.HSCs.ECs, 
                                             anchor.features = features)
fish.human.joint.HSCs.ECs <- IntegrateData(anchors = fish.human.anchors)

#perform integrated analysis
DefaultAssay(fish.human.joint.HSCs.ECs) <- "integrated"
joint.HSCs.ECs <- ScaleData(fish.human.joint.HSCs.ECs, verbose = FALSE)
joint.HSCs.ECs <- RunPCA(joint.HSCs.ECs, npcs = 30, verbose = FALSE)
ElbowPlot(object = joint.HSCs.ECs)
DimHeatmap(joint.HSCs.ECs, dims = 1:15, cells = 500, balanced = TRUE)
DimHeatmap(joint.HSCs.ECs, dims = 16:30, cells = 500, balanced = TRUE)

joint.HSCs.ECs <- RunUMAP(joint.HSCs.ECs, reduction = "pca", dims = 1:9)
joint.HSCs.ECs <- FindNeighbors(joint.HSCs.ECs, reduction = "pca", dims = 1:9)
joint.HSCs.ECs <- FindClusters(joint.HSCs.ECs, resolution = 0.8)

#visualization of clustering - Fig 3E
p2.joint <- DimPlot(joint.HSCs.ECs, reduction = "umap", label = TRUE, repel = TRUE, label.size = 7, pt.size = 2) +
  theme(legend.text = element_text(size = 20))
p2.joint
tiff("plot.tiff", units="in", width = 10, height = 7, res = 225)
p2.joint
dev.off()

#visualization of clustering - Fig S4A
p1.joint <- DimPlot(joint.HSCs.ECs, reduction = "umap", split.by = "group", label = FALSE, pt.size = 1.5) +
  theme(text = element_text(size = 22, face = "plain"), axis.text = element_text(size = 15), 
        axis.title = element_text(size = 15),legend.position = "bottom", 
        legend.text = element_text(size = 22))
p1.joint
tiff("plot.tiff", units = "in", width = 12, height = 9, res = 225)
p1.joint
dev.off()

saveRDS(joint.HSCs.ECs, "Zf_Human_HSCs_ECs_joint_atlas.rds")
joint.HSCs.ECs <- readRDS("ZF_Human_HSCs_ECs_joint_atlas.rds")

#identify cell type markers for each cluster and store in .csv - Dataset S11
DefaultAssay(joint.HSCs.ECs) <- "RNA"
recluster.markers <- FindAllMarkers(joint.HSCs.ECs, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
recluster.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
write.csv(recluster.markers, "Zf_Human_HSCs_ECs_joint_atlas_cluster_markers.csv")

#number of cells per cluster per group - Dataset S12
write.csv(table(Idents(joint.HSCs.ECs), joint.HSCs.ECs$group), "Zf_Human_HSCs_ECs_joint_atlas_cell_numbers_bygroup.csv")

#Dotplot of conserved markers for cell types of interest - Fig 3G
HSCs <- WhichCells(joint.HSCs.ECs, idents = c("zhEC-HSC","zhHSC"))

HSC.only <- subset(joint.HSCs.ECs, cells = HSCs, do.clean = TRUE, do.scale = TRUE)

markers.to.plot <- c("ACTA2","GFAP", "PDGFRB", "LRAT", "DES", "COL1A1",
                     "COL1A2", "SPARC", "ANGPTL6", "STEAP4", 
                     "HAND2", "COLEC11")

tiff("cluster_dotplot.tiff", units = "in", width = 7, height = 4.5, res = 225)
DotPlot(HSC.only, features = markers.to.plot, cols = c("orangered2", "mediumturquoise"), dot.scale = 8, 
        split.by = "group") + RotatedAxis() +
  theme(axis.text = element_text(size = 15, color = "black"),
        axis.text.x = element_text(face = "italic"), legend.position = "bottom")
dev.off()

#heatmap of expression of human endothelial cell markers - Fig S4B
markers.to.plot <- c("RAMP2", "CALCRL",
                     "F8", "PECAM1", "MGP", "SPARCL1", "TM4SF1", "CLEC14A",
                     "ID1", "IGFBP7", "ADIRF", "CTGF", "VWF", "CD9", "C7",
                     "SRPX", "ID3", "CAV1", "GNG11", "AQP1", "HSPG2", "EMP1",
                     "SOX18", "CLDN5",
                     "FCGR2B", "LYVE1", "STAB2", "CCL14", "CLEC1B", "FCN2",
                     "S100A13", "FCN3", "CRHBP", "STAB1", "IFI27",
                     "CLEC4G", "CCL23", "OIT3", "RAMP3", "SGK1", 
                     "DNASE1L3", "LIFR", "SPARC", "ADGRL4", "EGFL7", "PCAT19",
                     "CDKN1C",
                     "ENG", "INMT", "PTGDS", "TIMP3", "RNASE1", "GPM6A", 
                     "PTPRB", "FAM167B", "LTC4S")

cell.limit <- WhichCells(object = joint.HSCs.ECs, downsample = 150)
joint.HSCs.ECs.150cells <- subset(joint.HSCs.ECs, cells = cell.limit, do.clean = TRUE, do.scale = TRUE)

DefaultAssay(joint.HSCs.ECs.150cells) <- "RNA"
all.set <- rownames(joint.HSCs.ECs.150cells)
joint.HSCs.ECs.150cells <- ScaleData(joint.HSCs.ECs.150cells, display.progress = F,features = all.set)

library(viridis)

tiff("plot.tiff", units = "in", width = 10, height = 13, res = 225)
DoHeatmap(joint.HSCs.ECs.150cells, features = c(markers.to.plot), size = 6.5, angle = 90) + 
  NoLegend() + 
  scale_fill_viridis() +
  theme(axis.text.y = element_text(size = 16, color = "black", face = "italic"))
dev.off()


