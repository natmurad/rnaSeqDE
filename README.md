# :dna: Sequencing Data Analysis

## Docker Environment

**Environment**

 - Download [Snakemake Docker Image](https://hub.docker.com/r/snakemake/snakemake)

 - Start container

 ```docker run --name SnakemakeV1 -it <image_id> bash```
 
 - Copy your files to container

  ```docker cp <directory/files> <id_container>:<directory>```
 
 - Create environment from file *(in the container)*:
 
```conda env create -f environment.yaml```

## seqSE workflow

**Steps**:

1 - Quality control with *trimmomatic*

2 - Quality report with *fasqtc* and *multiqc*

3 - Align to reference with *STAR*

4 - Getting expression counts & TPM with *RSEM* (option: *featureCounts*)

**Run pipeline:**

```snakemake -s seqSE -c 50 --use-conda```

## Differential Expressed Genes with DESeq2

 - seqSE.R - R code to perform DE analysis using *DEseq2*, map ID and Gene information from Ensembl. Exploratory analysis of samples with heatmap, PCA Analysis, hierarchical clustering, MA plot and Volcano plot of DE genes. Creating tables with different filters, comparisons, UP and DOWN regulated genes.

## GO Enrichment Analysis


## Hallmarks

## References and Useful Links



