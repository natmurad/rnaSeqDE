# :dna: Sequencing Data Analysis

## seqSE workflow

**Steps**:

1 - Quality control with *trimmomatic*

2 - Quality report with *fasqtc* and *multiqc*

3 - Mapping to reference with *Star*

4 - Getting expression counts with *RSEM*

```snakemake -s seqSE -c 50 --use-conda```

**Environment**

 - Download [Snakemake Docker Image](https://hub.docker.com/r/snakemake/snakemake)

 - Start container

 ```docker run --name SnakemakeV1 -it <image_id> bash```
 
 - Copy your files to container

  ```docker cp <directory/files> <id_container>:<directory>```
 
 - Create environment from file *(in the container)*:
 
```conda env create -f environment.yaml```
