## Import libraries
library(tximport)
library(DESeq2)
library(tidyverse)
library(data.table)
library(AnnotationDbi)
#library(org.Mm.eg.db)
library(ggplot2)
library(gplots)
library(dplyr)
library(genefilter)
library(sva)
library(org.Hs.eg.db)
library(EnhancedVolcano)
library(magrittr)

work_dir = "/work/directory/"
setwd(work_dir)

# reading count files
counts_name <- list.files(file.path(work_dir), "*.genes.results")

txi.rsem <- tximport(counts_name, type = "rsem", txIn = FALSE, txOut = FALSE)
txi.rsem$length[txi.rsem$length == 0] <- 1
head(txi.rsem$counts)
dim(txi.rsem$counts)

# prepare metadata
metadata <- read.csv("metadata.csv", h = T, sep=",")

# gene names
genes <- rownames(txi.rsem$counts)
#genes <- gsub("gene-", "", genes)

# DE

dds <-  DESeqDataSetFromTximport(txi.rsem, colData = sampleTable,
                                 design = ~ condition)
dds <- DESeq(dds)

# filter
keep <- rowSums(counts(dds)) > 10
ddsfiltered <- dds[keep,]
nrow(ddsfiltered)
resultsNames(ddsfiltered)

# filtered

res <- results(ddsfiltered, contrast = c("condition", "1", "2"))
res <- lfcShrink(dds_1,
                 contrast = c("condition", "no", "BHB"),
                 res=res, type = 'ashr')

res$SYMBOL = mapIds(org.Hs.eg.db,
                    key = rownames(res),
                    column = "SYMBOL",
                    keytype = "ENSEMBL",
                    multivals = "first")

#res$ENTREZ = mapIds(org.Hs.eg.db,
#                    key = rownames(res),
#                    column = "ENTREZID",
#                    keytype = "ENSEMBL",
#                    multivals = "first")

library(EnhancedVolcano)
EnhancedVolcano(res,
                lab = res$SYMBOL,
                x = 'log2FoldChange',
                y = 'pvalue',
                xlim = c(-7,7),
                title = "Neurons - no LPS vs LPS",
                subtitle = "LogFC vs -Log10Pvalue"
)

# cluster
vsd <- vst(dds, blind = FALSE)
head(assay(vsd), 3)

library("genefilter")
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 200)
mat  <- assay(vsd)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[c("condition")])
pheatmap(mat, annotation_col = anno, show_rownames = F, border_color = FALSE)

## Gene Ontology

library(clusterProfiler)
library(enrichplot)
library(ggplot2)

# organism
organism = "org.Hs.eg.db"
#BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

original_gene_list <- res$log2FoldChange

# name the vector
names(original_gene_list) <- res$SYMBOL

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)


gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "SYMBOL", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "none")

# Dotplot
require(DOSE)
dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)

# network
x2 <- pairwise_termsim(gse)
emapplot(x2, showCategory = 10)


# categorySize can be either 'pvalue' or 'geneNum'
cnetplot(x2, categorySize="pvalue", foldChange=gene_list, showCategory = 5,
         cex_label_category = 0.5, cex_label_gene=0.5, node_label='all')