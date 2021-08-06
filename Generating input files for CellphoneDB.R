library(Seurat)

#load rds of clustered scRNA-seq data
scrnaseq.data <- readRDS("liverclustering_allwithMTcutoff_dim17_res0-8.rds")

#generate metadata file with columns "Cell" and "cell_type"
forCellphoneDB.metadata <- data.frame(scrnaseq.data@active.ident)
forCellphoneDB.metadata

names(forCellphoneDB.metadata)[1] <- "cell_type"

head(forCellphoneDB.metadata)

write.csv(forCellphoneDB.metadata, "metadata_forCellphoneDB.csv")

#generate counts file for CellpphoneDB
forCellphoneDB.counts <- scrnaseq.data@assays$RNA@data
forCellphoneDB.counts

head(forCellphoneDB.counts)

write.csv(forCellphoneDB.counts, "counts_forCellphoneDB.csv")

counts <- read.csv("counts_forCellphoneDB.csv")
colnames(counts) <- gsub("\\.+", "-", colnames(counts))

colnames(counts)
write.csv(counts, "counts_forCellphoneDB.csv")
