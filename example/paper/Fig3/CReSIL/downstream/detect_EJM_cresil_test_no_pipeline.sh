#!/bin/sh
#An example for serial job.
#SBATCH -J Cresil_EJM_no_pipeline
#SBATCH -o /home/lifesci/luosongwen/eccDNA/stimulate_data/test_data/logs/Cresil_EJM_no_pipeline_%j.log
#SBATCH -e /home/lifesci/luosongwen/eccDNA/stimulate_data/test_data/logs/Cresil_EJM_no_pipeline_%j.err
#SBATCH --qos=aepmemqos
#SBATCH --exclusive
#SBATCH --time=240:00:00
#SBATCH -N 1 -n 40 -p 2TB-AEP-Mem
echo Running on $SBATCH_PARTITION paratation
echo Start time is `date`
source /home/lifesci/luosongwen/miniconda3/etc/profile.d/conda.sh
conda activate ecc_pipe_old
#echo Directory is $PWD
echo This job runs on the following nodes:
echo $SLURM_JOB_NODELIST
echo This job has allocated $SLURM_JOB_CPUS_PER_NODE cpu core


cd /gpfs/home/lifesci/luosongwen/eccDNA/ecc_pipe
R1fastq=/home/lifesci/luosongwen/eccDNA/stimulate_data/test_data/Cresil_real_data/cresil_fastq/SRR18143375_1.fastq

## Run trim 
cresil trim -t 40 -fq $R1fastq -r /gpfs/home/lifesci/luosongwen/eccDNA/ecc_pipe/resource/cresil-master/reference/hg38/hg38.mmi -o Cresil_EJM_test_result

## Run eccDNA identification for enriched data
cresil identify -t 40 -fa /gpfs/home/lifesci/luosongwen/eccDNA/ecc_pipe/resource/cresil-master/reference/hg38/hg38.fa -fai /gpfs/home/lifesci/luosongwen/eccDNA/ecc_pipe/resource/cresil-master/reference/hg38/hg38.fa.fai -fq $R1fastq -trim Cresil_EJM_test_result/trim.txt

## Run eccDNA annotation
cresil annotate -t 40 -rp reference.rmsk.bed -cg reference.cpg.bed -gb reference.gene.bed -identify Cresil_EJM_test_result/eccDNA_final.txt

echo End time: `date`


