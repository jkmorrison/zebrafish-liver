#read in RDS
HSC.subset.humanorthologs <- readRDS("HSC_subset_clustering_HumanOrthologs.rds")

HSC.subset.humanorthologs <- RenameIdents(HSC.subset.humanorthologs, `0` = "Inflammatory Macrophages",
                                      `1` = "Endothelial Cells 1", `2` = "NK/T Cells 1", 
                                      `3` = "Endothelial Cells 2", `4` = "EC-HSC Doublets", 
                                      `5` = "Endothelial Cells 3", `6` = "Neutrophils", 
                                      `7` =  "NK/T Cells 2", `8` = "Hepatic Stellate Cells",
                                      `9` = "Non-inflammatory Macrophages", `10` = "agr2+ Cholangiocytes",
                                      `11` = "apln+ Cells", `12` = "EC-Hepatocyte Doublets") 

#split Seurat object by group into list of WT and MPI MT Seurat objects
split <- SplitObject(HSC.subset.humanorthologs, split.by = "group")

#pull cell barcodes for both group
WT.cells <- split$WT@assays$RNA@counts@Dimnames[[2]]
MPIMT.cells <- split$MPI_MT@assays$RNA@counts@Dimnames[[2]]

#subset HSC subset Seurat object by group (genotype) by using pulled cell barcodes
WT.HSC.subset <- subset(HSC.subset.humanorthologs, cells = WT.cells, do.clean = TRUE, do.scale = TRUE)
MPIMT.HSC.subset <- subset(HSC.subset.humanorthologs, cells = MPIMT.cells, do.clean = TRUE, do.scale = TRUE)

#save each as R data file
saveRDS(WT.HSC.subset, "WT_HSC_subset_HumanOrthologs.rds")
saveRDS(MPIMT.HSC.subset, "MPIMT_HSC_subset_HumanOrthologs.rds")



