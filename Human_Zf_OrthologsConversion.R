#load in zebrafish and human Seurat objects
fish <- readRDS("WTzebrafishliver_03-18-21.rds")

load("HumanLiver.Rdata")
human <- HumanLiverSeurat

#read the orthology file and prepare 
mart_export <- read.csv("HumanFish_Orthologs_v1.csv", row.names=NULL)
mart1 <- mart_export[,c(6,8)]
mart1$id <- paste0(mart1$Human.gene.name, "-", mart1$Gene.name)
mart2 <- subset(mart1, !duplicated(mart1$id))

write.csv(mart2, "human_zf_table.csv")

#import table with manual edits for genes of importance
mart2 <- read.csv("human_zf_table.csv")

#filter the duplicates
mart2 <- filter(mart2, !duplicated(Gene.name))
mart2 <- filter(mart2, !duplicated(Human.gene.name))

#subset data to the zebrafish genes with human orthologs
humangenes <- as.vector(unlist(mart2$Human.gene.name))
fishgenes <- as.vector(unlist(mart2$Gene.name))
fish1 <- subset(fish, features = fishgenes)
human1 <- subset(human, features = humangenes)

#convert gene names in RNA@counts and RNA@data slots to human orthologs
for (i in fishgenes){
  pn <- paste0("^", i, "$")
  p = which(grepl(pn, mart2$Gene.name))
  p1 <- p[1]
  f <- mart2$Human.gene.name[p1]
  pn2 <- i
  pn3 <- paste0("^", pn2, "$")
  k = grep(pn3, fish1@assays$RNA@counts@Dimnames[[1]])
  k1 <- k[1]
  fish1@assays$RNA@counts@Dimnames[[1]][k1] <- as.character(f)}

for (i in fishgenes){
  pn <- paste0("^", i, "$")
  p = which(grepl(pn, mart2$Gene.name))
  p1 <- p[1]
  f <- mart2$Human.gene.name[p1]
  pn2 <- i
  pn3 <- paste0("^", pn2, "$")
  k = grep(pn3, fish1@assays$RNA@data@Dimnames[[1]])
  k1 <- k[1]
  fish1@assays$RNA@data@Dimnames[[1]][k1]<-as.character(f)}

#convert meta.features gene names to human orthologs
fish1[["RNA"]]@meta.features <- data.frame(row.names = rownames(fish1[["RNA"]]))

#save the rds with gene names changed
saveRDS(fish1, file="WTFishLiverClustering_HumanOrthologs.rds")
