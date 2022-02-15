library(DESeq2)
library(org.Hs.eg.db)
library(tibble)
library(dplyr)
library(tidyr)
library(fgsea)
library(ggplot2)
library(reshape2)
library(ComplexHeatmap)
library(circlize)


# already we performed DESeq2 analysis and have statitics for workig on it
res$row <- rownames(res)

# important notice: if you have not such stats in your result (say comming from edgeR),
# you may need to create a rank metric for your genes. To do this:
res$stat = -log10(res$pvalue)/sign(res$log2FoldChange)

# Map Ensembl gene IDs to symbol. First create a mapping table.
ens2symbol <- AnnotationDbi::select(org.Hs.eg.db,
                                    key=rownames(res), 
                                    # key=res$SYMBOL,
                                    columns="SYMBOL",
                                    keytype="ENSEMBL")
names(ens2symbol)[1] <- "row"

ens2symbol <- as_tibble(ens2symbol)
ens2symbol

# joining
res <- merge(data.frame(res), ens2symbol, by=c("row"))
#res$SYMBOL <- res$SYMBOL.x
res2 <- res %>% 
  dplyr::select(SYMBOL, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(SYMBOL) %>% 
  summarize(stat=mean(stat))
#res2
res2$stat[which(res2$stat==-Inf)] <- 0
res2$stat[which(res2$stat==Inf)] <- 0

# creating  a named vector [ranked genes]
ranks <- res2$stat
names(ranks) <- res2$SYMBOL
ranks <- na.omit(ranks)
pathways.hallmark <- gmtPathways("/pathwayslist/h.all.v7.4.symbols.gmt")

head(pathways.hallmark)

####### IF NOT HUMAN NEED TO MAP, IF HUMAN SKIP
# creating  a named vector [ranked genes]
library(biomaRt)
mart <- useDataset("mmusculus_gene_ensembl", mart=useMart("ensembl"))
bm <- getBM(attributes=c("ensembl_gene_id", "hsapiens_homolog_associated_gene_name"), mart=mart) %>%
  distinct() %>%
  as_tibble() %>%
  na_if("") %>% 
  na.omit()

res3 <- merge(res, bm, by.x="row", by.y="ensembl_gene_id")

ranks <- res3$stat
names(ranks) <- res3$hsapiens_homolog_associated_gene_name
ranks <- na.omit(ranks)
pathways.hallmark <- gmtPathways("/aln/counts/h.all.v7.4.symbols.gmt")

######################

#Running fgsea algorithm:
fgseaRes <- fgseaMultilevel(pathways=pathways.hallmark, stats=ranks)

# Tidy the results:
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES)) # order by normalized enrichment score (NES)

# To see what genes are in each of these pathways:
gene.in.pathway <- pathways.hallmark %>% 
  enframe("pathway", "SYMBOL") %>% 
  unnest(cols = c(SYMBOL)) %>% 
  inner_join(res, by="SYMBOL")

# Plot the normalized enrichment scores - BarPlot
# Color the bar indicating whether or not the pathway was significant:
fgseaResTidy$adjPvalue <- ifelse(fgseaResTidy$padj <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")
ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways Enrichment Score from GSEA")
