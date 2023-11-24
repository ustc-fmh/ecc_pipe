# Cross-Species User Add

## 1.download genome fa file
```
## ucsc
https://hgdownload.soe.ucsc.edu/downloads.html
## ncbi
https://www.ncbi.nlm.nih.gov/datasets/genome/
```
## 2.upload fa file to acc_pipe
```
## upload fa file in resource
$ecc_pipe/resource/user_add/reference
## prepare files
cd $ecc_pipe/resource/user_add/reference
bwa index $genome.fa
minimap2 -d $genome.mmi $genome.fa
samtools faidx $genome.fa
```
## 3.choose params in Detect mode
```
## user_ref params in config.yaml
user_ref = $genome
```
## Here is a example for TAIR10 genome
* [Example TAIR10: Shell Code](TAIR10.pbs)
* [Example TAIR10: config](TAIR10.yaml)

## NOTE
-   AA is not support this function, if you have more question about AA, you can see issues below.
-   [issue](https://github.com/AmpliconSuite/AmpliconSuite-pipeline/issues/45)