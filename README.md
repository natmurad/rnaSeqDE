# :dna: Sequencing Data Analysis

## Docker Environment

 - Download the Dockerfile and inside same directory run:

 ```docker build -t <IMAGE_TAG> .```

 - Check the IMAGE ID (look for the name of your image tag):

 ```docker images -a```
 
 - Start container

 ```docker run --name <CONTAINER_NAME> -it <IMAGE_ID> bash```
 
 - Copy your files to container

  ```docker cp <DIRECTORY/FILES> <CONTAINER_ID>:<DIRECTORY>```
 
## seqSE Workflow

**Steps**:

1 - Quality control with *trimmomatic*

2 - Quality report with *fasqtc* and *multiqc*

3 - Align to reference with *STAR*

4 - Getting expression counts & TPM with *RSEM* (option: *featureCounts*)


## Prepare code before running:

### Change in the file *"seqSE"* the path to the *config.yaml* file.

### In *config.yaml* file change the paths for all folder.


## Run pipeline:

```snakemake -s seqSE -c 50```


## Differential Expressed Genes with DESeq2

 - seqSE.R - R code to perform DE analysis using *DEseq2*, map ID and Gene information from Ensembl. Exploratory analysis of samples with heatmap, PCA Analysis, hierarchical clustering, MA plot and Volcano plot of DE genes.
 
 - nasaCode.R - [NASA Gene Lab code](https://github.com/nasa/GeneLab_Data_Processing/tree/master/RNAseq/GLDS_Processing_Scripts/GLDS-168/04-05-DESeq2_NormCounts_DGE) that I like to use.

## GO Analysis

 - GOanalysis.R - GO analysis using *clusterprofiler* & *enrichplot*

## Pathways

- fgsea.R - enriched pathways using *fgsea*

## References and Useful Links

[NASA Gene Lab script](https://github.com/nasa/GeneLab_Data_Processing/tree/master/RNAseq/GLDS_Processing_Scripts/GLDS-168/04-05-DESeq2_NormCounts_DGE)

[DESeq2 experimental design and interpretation](https://rstudio-pubs-static.s3.amazonaws.com/329027_593046fb6d7a427da6b2c538caf601e1.html#example-1-two-group-comparison) *Steven Xijin Ge*

[Gene Set Enrichment Analysis with ClusterProfiler](https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/) *Mohammed Khalfan*

[Analyzing RNA-seq data with DESeq2](http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#ma-plot) *Michael I. Love, Simon Anders, Wolfgang Huber*

[DE Analysis script](https://github.com/ACSoupir/Bioinformatics_RNASeq/blob/master/Mouse_RNA_Seq_p53_genotoxic.Rmd) *Alex Soupir*

[RNA-seq Analysis with R](https://sbc.shef.ac.uk/workshops/2019-01-14-rna-seq-r/rna-seq-gene-set-testing.nb.html) *Stephane Ballereau, Mark Dunning, Oscar Rueda, Ashley Sawle*
