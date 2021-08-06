setwd("~/Documents/NYUAD/MPI/HumanCell_MPI")
directory <- getwd()
sampleFiles <- grep('.txt', list.files(directory), value=T)
sampleCondition <- c('ctrl', 'siMPI', 'ctrl', 'siMPI')
sampleNames <- sub('_rawCounts.txt', '' , sampleFiles)
sampleTable <- data.frame(sampleName = sampleNames,
                          fileName = sampleFiles,
                          condition = sampleCondition)
library(DESeq2)
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design= ~ condition)
### pre-filtering: remove rows whose counts <= 10
ddsHTSeq <- ddsHTSeq[rowSums(counts(ddsHTSeq)) > 10, ]
### set untreated or control first so that the direction of the logs fold changes compared to control
colData(ddsHTSeq)$condition<-factor(colData(ddsHTSeq)$condition, 
                                    levels=c('ctrl','siMPI'))

dds <- estimateSizeFactors(ddsHTSeq)
dds <- estimateDispersions(dds)
dds <- DESeq(ddsHTSeq)
### transformation
rld <- rlog(dds, blind=FALSE)
head(assay(rld), 3)
library("pheatmap")
sampleDists <- dist(t(assay(rld)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- sampleNames
colnames(sampleDistMatrix) <- sampleCondition
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

Counts_dds <- counts(dds, normalized=TRUE)
colnames(Counts_dds) <- rownames(sampleDistMatrix)
library("IHW")
resIHW <- results(dds, filterFun=ihw)
summary(resIHW)
plotMA(resIHW, ylim=c(-8,8),alpha = 0.05, main='siMPI vs. ctrl - Human Cell', cex=1.0)


library('biomaRt')
listMarts()    # to see which database options are present
ensembl <- useMart('ensembl')  # using ensembl database data
listDatasets(ensembl)     # function to see which datasets are present in ensembl
mart_hg <- useDataset('hsapiens_gene_ensembl', useMart('ensembl'))
# human: 'hgnc_symbol', mouse: 'mgi_symbol', zebrafish: 'zfin_id_symbol' 
Ensembl_ids <- as.character(row.names(Counts_dds))
Gene_symbols <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"),
                      filters='ensembl_gene_id',values=Ensembl_ids, mart=mart_hg) 

#Gene_GO <- getBM(attributes=c("ensembl_gene_id", "mgi_symbol",'go_id',"interpro_description"),
#                filter='ensemble_gene_id', values = Ensembl_ids, mart=mart_mm)

### Summary
Dat <- cbind(as.data.frame(resIHW), Counts_dds)
Dat$ensembl_id <- row.names(Counts_dds)
Dat <- merge(Dat, Gene_symbols, by.x='ensembl_id', by.y='ensembl_gene_id')
### output
write.table(Dat, 'HumanCell_siMPI_vs_ctrl_DEG.csv', quote=F, sep='\t')
