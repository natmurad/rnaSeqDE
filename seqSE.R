## Import libraries
library(tximport)
library(DESeq2)
library(tidyverse)
library(data.table)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(ggplot2)
library(gplots)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(PoiClaClu)
library(genefilter)
library(sva)
library(glmpca)

# set directory
work_dir="/work/directory/counts/"
setwd("/work/directory/counts/")

# reading count files
counts_name <- list.files(file.path(work_dir), "*.genes.results")
txi.rsem <- tximport(counts_name, type = "rsem", txIn = FALSE, txOut = FALSE)
txi.rsem$length[txi.rsem$length == 0] <- 1
head(txi.rsem$counts)
dim(txi.rsem$counts)

# prepare metadata
metadata <- read.csv("metadata.tsv", h = T, sep=",")
sampleTable <- data.frame(condition = factor(rep(c("Control", "Control_Treat1",
                                                   "Treat2", "Treat2_Treat1"),
                                                 each = 3)))
colData <- metadata[,c(1,3,4)]
colData <- cbind(sampleTable, colData)

# load rsem object
rownames(sampleTable) <- colnames(txi.rsem$counts)
dds <- DESeqDataSetFromTximport(txi.rsem,
                                colData = colData,
                                design = ~ treat1 + treat2 + treat1*treat2)
# differential expression
dds <- DESeq(dds)
res <- results(dds)

# filter
keep <- rowSums(counts(dds)) > 1
ddsfiltered <- dds[keep,]
nrow(ddsfiltered)

## sample dists

# variance stabilizing transformation (VST)
vsd <- vst(dds, blind = FALSE)
head(assay(vsd), 3)

sampleDists <- dist(t(assay(vsd)))
sampleDists

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$treat1, vsd$treat2, sep = " - ")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")))(255)

pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

# poisson dist
poisd <- PoissonDistance(t(counts(dds)))
samplePoisDistMatrix <- as.matrix(poisd$dd)
rownames(samplePoisDistMatrix) <- paste( dds$treat1, dds$treat2, sep=" - ")
colnames(samplePoisDistMatrix) <- NULL
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors)

# PCA plot using Generalized PCA
gpca <- glmpca(counts(ddsfiltered), L=2)
gpca.dat <- gpca$factors
gpca.dat$gravity <- ddsfiltered$treat1
gpca.dat$stim <- ddsfiltered$treat2

ggplot(gpca.dat, aes(x = dim1, y = dim2, color = treat1, shape = treat2)) +
  geom_point(size =3) + coord_fixed() + ggtitle("glmpca - Generalized PCA")

# pca log data
rld <- rlog(dds)
rldfilter <- rlog(ddsfiltered)
plotPCA(rld, intgroup = c("treat1","treat2"))
plotPCA(rldfilter, intgroup = c("treat1","treat2"))

## look at results
resNames <- results(dds)
head(res)
summary(res)

res <- results(ddsfiltered, lfcThreshold = 0.01)

# which genes are different? blue points are significant and gray points are not logfoldchange
plotMA(res, ylim=c(-5,5))

# add significance col
res1 <- as.data.frame(res)
res1 <- mutate(res1, sig=ifelse(res1$padj<0.1, "FDR<0.1", "Not sig"))

res1[which(abs(res1$log2FoldChange)<1.0), "sig"] = "Not sig"

# changing refseq to gene symbols
res$refseq <- rownames(res1)

#columns(org.Mm.eg.db)
#org.Mm.eg.db

res$SYMBOL = mapIds(org.Mm.eg.db,
                    key = res1$refseq,
                    column = "SYMBOL",
                    keytype = "ENSEMBL",
                    multivals = "first")

res$ENTREZ = mapIds(org.Mm.eg.db,
                    key = res1$refseq,
                    column = "ENTREZID",
                    keytype = "ENSEMBL",
                    multivals = "first")

# contrasts
resTreat1 <- results(ddsfiltered, contrast=c("treat1","control", "treat"))
summary(resTreat1)
resStim$SYMBOL = mapIds(org.Mm.eg.db,
                        key = res1$refseq,
                        column = "SYMBOL",
                        keytype = "ENSEMBL",
                        multivals = "first")

resStim$ENTREZ = mapIds(org.Mm.eg.db,
                        key = res1$refseq,
                        column = "ENTREZID",
                        keytype = "ENSEMBL",
                        multivals = "first")

treat1Sig <- subset(resTreat1, padj < 0.05)

resTreat2 <- results(ddsfiltered, contrast=c("treat2","control", "treat"))
summary(resTreat2)
resTreat2$SYMBOL = mapIds(org.Mm.eg.db,
                           key = res1$refseq,
                           column = "SYMBOL",
                           keytype = "ENSEMBL",
                           multivals = "first")

resTreat2$ENTREZ = mapIds(org.Mm.eg.db,
                           key = res1$refseq,
                           column = "ENTREZID",
                           keytype = "ENSEMBL",
                           multivals = "first")

treat2Sig <- subset(resTreat2, padj < 0.05)

# multiple testing
sum(res$pvalue < 0.05, na.rm=TRUE)
sum(!is.na(res$pvalue))

sum(resTreat1$pvalue < 0.01, na.rm=TRUE)
sum(!is.na(resTreat1$pvalue))

sum(resTreat2$pvalue < 0.01, na.rm=TRUE)
sum(!is.na(resTreat2$pvalue))

# strongest down-reg
resSig <- subset(res, padj < 0.1)
head(resSig[ order(resSig$log2FoldChange), ])

# strongest up-reg
head(resSig[ order(resSig$log2FoldChange, decreasing = TRUE), ])

# gene clustering
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 20)

mat  <- assay(vsd)[topVarGenes,]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("gravity","stim")])
pheatmap(mat, annotation_col = anno)

resSig <- subset(res, padj < 0.05)


write.table(resSig, "DEall", sep = "\t", col.names = TRUE, row.names = TRUE, dec = ".")
write.table(treat1Sig, "DEtreat1", sep = "\t", col.names = TRUE,
            row.names = TRUE, dec = ".")
write.table(treat2Sig, "DEtreat2", sep = "\t", col.names = TRUE, row.names = TRUE, dec = ".")
