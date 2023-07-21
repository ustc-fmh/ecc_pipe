# eccDNA_pipe: A pipeline for detecting eccDNA from NGS/TGS  
A pipeline for detecting Extrachromosomal Circular DNA from NGS/TGS in Circle-Map, AmpliconArchitect, CReSIL

## Installation instructions
### Install Conda environment
```
cd ecc_pipe
conda env create -f ./install/env.yml
```
## resource and AA sample data share
Please install resource.zip and unzip in the ./ecc_pipe before set environment
- [resource](https://rec.ustc.edu.cn/share/1db60bb0-2705-11ee-b797-2951fa8ccb8e)
- [AA sample data](https://rec.ustc.edu.cn/share/01b661c0-ceeb-11ed-a387-976f74b3f711)
- passwd:ustc

### Activate environment and set environment for AA
```
conda activate ecc_pipe
## Set environment for AA
master=$PWD
#### data_repo
cd ./resource/AA/AmpliconArchitect/data_repo/
echo export AA_DATA_REPO=$PWD/ >> ~/.bashrc
touch coverage.stats && chmod a+rw coverage.stats
cd $master
#### AC
echo export AC_SRC=$PWD/resource/AA/AmpliconClassifier >> ~/.bashrc
#### AA
echo export AA_SRC=$PWD/resource/AA/AmpliconArchitect/src >> ~/.bashrc
#### mosek
cd ./resource/AA/AmpliconArchitect/
echo export MOSEKPLATFORM=linux64x86 >> ~/.bashrc
export MOSEKPLATFORM=linux64x86
echo export PATH=$PATH:$PWD/mosek/8/tools/platform/$MOSEKPLATFORM/bin >> ~/.bashrc
echo export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD/mosek/8/tools/platform/$MOSEKPLATFORM/bin >> ~/.bashrc
echo export MOSEKLM_LICENSE_FILE=$PWD/mosek/8/licenses >> ~/.bashrc
cd $PWD/mosek/8/tools/platform/linux64x86/python/3/
python3 setup.py install
```
### Install CReSIL and install reference for homer
```
cd $master
cd ./resource/cresil-master
pip install .

## Install homer reference
#### $conda_path = ~/miniconda3/envs/ecc_pipe
cd $conda_path/share/homer
perl configureHomer.pl -install mm10
perl configureHomer.pl -install hg38
```

## 01.QC
First set all params in ***QC_config.yaml*** and run python in shell for QC by fastqc and trim_galore
```
python3 ecc_pipe_master.py --QC -n 4 --config configfile_path
## --QC: set the function; no input
## -n: set max threads for all data
## --config: set QC_config.yaml file path 
```
all params in the ***QC_config.yaml*** file show:
-   **ecc_master_path** - The full path for ecc_pipe; eg: /home/user/project/ecc_pipe
-   **output_path** - output path, eg: ./example/01.upstream/result/01.QC
-   **input_file** - fastq file path, eg: file_name: ['r1.fq.gz', 'r2.fq.gz']
-   **threads** - set max threads for one data, eg:1
-   **paired** - single or paired sequence; eg: True or False

## 02.Detect
First set all params in ***circlemap_config.yaml, AA_config.yaml, cresil_config.yaml***
```
python3 ecc_pipe_master.py --Detect --tool circlemap -n 24 --config configfile_path
python3 ecc_pipe_master.py --Detect --tool AA -n 24 --config configfile_path
python3 ecc_pipe_master.py --Detect --tool cresil -n 24 --config configfile_path
## --Detect: set the function; no input
## --tool: set detect tools for NGS/TGS
## -n: set max threads for all data
## --config: set tool_config.yaml file path
```
all params in the ***circlemap_config.yaml*** file show:
-   **ecc_master_path** - The full path for ecc_pipe; eg: /home/user/project/ecc_pipe
-   **output_path** - output path, eg: ./example/01.upstream/result/02.circlemap
-   **input_file** - fastq file path, eg: file_name: ['r1.fq.gz', 'r2.fq.gz']
-   **threads** - set max threads for one data, eg:12
-   **reference** - hg38 or mm10


all params in the ***AA_config.yaml*** file show:
-   **ecc_master_path** - The full path for ecc_pipe; eg: /home/user/project/ecc_pipe
-   **output_path** - output path, eg: ./example/01.upstream/result/02.AA
-   **input_file** - fastq file path, eg: file_name: ['r1.fq.gz', 'r2.fq.gz']
-   **threads** - set max threads for one data, eg:12
-   **reference** - hg38 or mm10
-   **cnvkit_dir** - the path for cnvkit eg: $conda_path/bin/cnvkit.py
-   **rscript** - Rscript path eg: $conda_path/bin/Rscript
-   **python_path** - Python path eg: $conda_path/bin/python
-   **cngain** - Set a custom threshold for the CN gain considered by AA, default:4
-   **cnsize** - Set a custom threshold for CN interval size considered by AA, default:10000


all params in the ***cresil_config.yaml*** file show:
-   **ecc_master_path** - The full path for ecc_pipe; eg: /home/user/project/ecc_pipe
-   **output_path** - output path, eg: ./example/01.upstream/result/02.cresil
-   **input_file** - fastq file path, eg: file_name: ['r1.fq.gz', 'r2.fq.gz']
-   **threads** - set max threads for one data, eg:12
-   **reference** - hg38 or mm10
-   **wgs** - if use WGS mode, default: 0

Detect tools raw website:
- [circlemap](https://github.com/iprada/Circle-Map)
- [AA](https://github.com/jluebeck/AmpliconSuite-pipeline)
- [cresil](https://github.com/visanuwan/cresil)

## 03.Analysis
we apply three mode to analysis the eccDNA result file by python: Distribution, DEG, Visualize
### Distribution
This mode calculates the length and chromatin distribution and 
        performs annotation of RepeatMasker with Homer, 
        as well as enhancers, super-enhancers, SNPs, and eQTLs with [eccAtlas](http://lcbb.swjtu.edu.cn/eccDNAatlas/)  
```
python3 ecc_pipe_master.py --Analysis --mode Distribution --tool circlemap \
            --file_path example/02.downstream/circlemap/DNARCAT1.circle.bed \
            --geno hg38 --trim 0.5
```
-   **Analysis** - set the function; no input
-   **mode** - set the mode, eg: Distribution
-   **tool** - str in ['circlemap', 'AA', 'cresil', 'other']
-   **file_path** - eccDNA result file from 02.Detect; if set tool=='other',please set the file_path contain 6 columns:['Chr', 'Start', 'End', 'Count', 'eccID', 'Length']; eg: ./example/02.downstream/other/result_ecc.txt
-   **peak_path** - default: no input; if tool=='other', please set the peak bed file contain 5 columns: ['chr', 'start', 'end', 'count', 'id'] eg: ./example/02.downstream/other/result_peak.bed
-   **geno** - hg38 or mm10
-   **trim** - overlap ratio for annotate enhancer/super_enhancer/snp/eQTL ; float in [0,1]; eg:0.5
-   **circlemap_qc** - if set tool == 'circlemap', on/off qc the eccDNA[note: if tool != circlemap not need this params] result eg:1
    
### DEG
Prior to running this mode, the Distribution mode must be executed first. 
        Afterwards, this mode calculates the gene matrix by intersecting the eccDNA matrix, 
        and performs differential expression analysis (DEG) and Gene Ontology (GO) annotation using Deseq2/clusterprofile.
```
python3 ecc_pipe_master.py --Analysis --mode DEG \
        --path_share example/02.downstream/deg_test/ \
        --group_file example/02.downstream/deg_test/group.txt \
        --count_type gene --geno hg38 --trim 1 --log2fc 2 --pvalue 0.01
```
-   **Analysis** - set the function; no input
-   **mode** - set the mode, eg: DEG
-   **path_share** - path contain all sample data, eg: example/02.downstream/deg_test/
-   **group_file** - txt file contain group info, eg: example/02.downstream/deg_test/group.txt
-   **count_type** - gene or region, gene: compute all gene in ecc region; region: compute ecc region in gene
-   **trim** - overlap ratio for compute gene matrix ; float in [0,1]; eg: 1    
-   **log2fc** - log2foldchange by deseq2 result; eg: 1 
-   **pvalue** - pvalue cut by deseq2 result; eg: 0.05 
        
### Visualize
this mode Visualize the eccDNA by Circlize
```
python3 ecc_pipe_master.py --Analysis --mode Visualize \
        --bed_file example/02.downstream/cresil/ecc_pipe_result/cresil_result.analysis.bed \
        --geno hg38 --ecc_id 1
```
-   **Analysis** - set the function; no input
-   **mode** - set the mode, eg: Visualize
-   **bed_file** - peak_list bed file;
-   **geno** - hg38 or mm10
-   **ecc_id** - ecc id in the peak_list

## single sample test run
Please read pbs:
    [upstream](https://github.com/ustc-fmh/ecc_pipe/tree/main/example/01.upstream/script)
    [downstream](https://github.com/ustc-fmh/ecc_pipe/tree/main/example/02.downstream/script)

## Multi sample analysis 
Please read ipynb:
    [Multi sample](https://github.com/ustc-fmh/ecc_pipe/blob/main/example/github_online/Multi_sample_analysis.ipynb)

## NOTE
-   **1** - The YAML file delimiter should not be a tab; it should consist of four spaces.
-   **2** - The reference supports "hg38" and "mm10".
-   **3** - The input file for Cresil is fq but not fq.gz.

## Citation
Please cite the following article if you use eccDNA_pipe in your research
> xxxx

## License and Copyright
eccDNA_pipe is distributed under the terms of the USTC