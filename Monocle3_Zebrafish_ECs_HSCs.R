library(monocle3)
library(ggplot2)
library(dplyr)
library(SeuratWrappers)
library(Seurat)
library(patchwork)
library(magrittr)
library(tidyverse)
library(tibble)
library(pheatmap)
library(viridis)

#load in zebrafish subset RDS and select hepatocyte populations
zebrafish.HSC.subset <- readRDS("AllZf_HSC_subset.rds")

#subset HSCs and ECs
HSC.EC.cells <- WhichCells(zebrafish.HSC.subset, idents = c("zEC1", "zEC2", "zEC-HSC", "zEC3", "zHSC"))
HSC.EC <- subset(zebrafish.HSC.subset, cells = HSC.EC.cells, do.clean = TRUE, do.scale = TRUE)

VlnPlot(HSC.EC, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"))
DimPlot(object = HSC.EC, reduction = "umap", group.by = "group", label = TRUE, repel = TRUE, label.size = 3.9, pt.size = 0.5)

#split Seurat object by group into list of WT and MPI MT Seurat objects
split <- SplitObject(HSC.EC, split.by = "group")

#pull cell barcodes for both group
WT.cells <- split$WT@assays$RNA@counts@Dimnames[[2]]

#subset HSC subset Seurat object by group (genotype) by using pulled cell barcodes
WT.HSC.EC <- subset(HSC.EC, cells = WT.cells, do.clean = TRUE, do.scale = TRUE)

#create cds object for Monocle analysis
data <- as(as.matrix(WT.HSC.EC@assays$RNA@data), 'sparseMatrix')

pd <- data.frame(WT.HSC.EC@meta.data)

fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))

cds <- new_cell_data_set(data,
                         cell_metadata = pd,
                         gene_metadata = fData)

cds <- preprocess_cds(cds, num_dim = 100)
plot_pc_variance_explained(cds)
cds <- reduce_dimension(cds)

cds <- cluster_cells(cds)

#assign cluster info from Seurat
list_cluster <- WT.HSC.EC@active.ident
names(list_cluster) <- WT.HSC.EC@assays[["RNA"]]@data@Dimnames[[2]]

cds@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster

### Assign UMAP coordinates from Seurat
cds@int_colData@listData$reducedDims@listData$UMAP <-WT.HSC.EC@reductions$umap@cell.embeddings

#confirm proper clustering and label - Fig S3A
plot <- plot_cells(cds, color_cells_by = "cluster", label_branch_points = FALSE, label_leaves = FALSE,
           group_label_size = 8, show_trajectory_graph = FALSE, cell_size = 2) +
  xlim(-13, -2) +
  ylim(-6, 4)
plot

tiff("plot.tiff", units = "in", height = 5.5, width = 7, res = 230)
plot
dev.off()

#Pseudotime on all partitions maintaing partition separations
cds <- learn_graph(cds)

#Moran's I differential gene expression analysis - Fig 3D
###graph autocorrelation analysis for comparing clusters
pr_graph_test_res <- graph_test(cds, neighbor_graph = "knn")
pr_deg_ids <- row.names(subset(pr_graph_test_res, q_value < 0.05))
sig_deg_graph_test <- pr_graph_test_res %>% filter (q_value < 0.05)
write.csv(sig_deg_graph_test, "sig_deg_bycluster_graph_test.csv")

###find modules of co-regulated genes by cluster
cds2 <- cds[pr_deg_ids,]
gene_module_df <- find_gene_modules(cds2, resolution = 1e-3)
write.csv(gene_module_df, "cluster_gene_modules.csv")

cell_group_df <- tibble::tibble(cell=row.names(colData(cds2)), 
                                cell_group=clusters(cds2)[colnames(cds2)])
agg_mat <- aggregate_gene_expression(cds2, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))

#creates heatmap of gene modules for clusters
plot1 <- pheatmap::pheatmap(agg_mat, cluster_rows=TRUE,
                            scale="column", clustering_method="ward.D2",
                            fontsize=13, number_color = "black")

plot1

tiff("plot.tiff", units = "in", height = 6, width = 8.5, res = 225)
plot1
dev.off()

#plots gene module onto clustering - Fig S3B
plot1 <- plot_cells(cds2, 
                    genes = gene_module_df %>% filter(module %in% 10),
                    group_cells_by = "cluster",
                    color_cells_by = "cluster",
                    show_trajectory_graph = FALSE,
                    label_cell_groups = FALSE,
                    cell_size = 1.5) +
  theme(text = element_text(size = 18),legend.text = element_text(size = 14),
        legend.title = element_text(size = 16), axis.text = element_text(size = 14,
                                                                         color = "black"),
        axis.title = element_text(size = 14)) +
  xlim(-13, -2) +
  ylim(-6, 4)

plot1

tiff("plot.tiff", units = "in", height = 5.5, width = 7.5, res = 230)
plot1
dev.off()

##heatmaps of co-expression module genes - Fig S3C
#get avg expression of genes for clusters
avgexp = AverageExpression(WT.HSC.EC, return.seurat = T)

#heatmap of gene module genes
module_genes <- gene_module_df %>% filter(module %in% 5)

plot1 <- DoHeatmap(avgexp, features = module_genes$id, size = 6, angle = 90,
                   draw.lines = FALSE) + 
  NoLegend() +
  scale_fill_viridis() +
  theme(axis.text.y = element_blank())
plot1

tiff("plot.tiff", units = "in", height = 9, width = 6.5, res = 225)
plot1
dev.off()

saveRDS(cds, "WT_Zf_HSCs_ECs_monocle3_analysis.rds")

cds <- readRDS("WT_Zf_HSCs_ECs_monocle3_analysis.rds")
write.table(get_citations(cds), "monocle_citations.txt")
